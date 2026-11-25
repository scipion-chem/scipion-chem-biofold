# **************************************************************************
# *
# * Authors:   Joaquin Algorta (joaquin.algorta@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os
import shutil

import pyworkflow.protocol.params as params

from pwchem import Plugin

from biofold.constants import CHAI_DIC
from .protocol_cofolding_base import COFOLD_ATTRS_HELP, ProtCofoldingBase


class ProtChaiCofolding(ProtCofoldingBase):
    # assigned, not a bare literal: a """...""" % VAR expression is not a docstring
    __doc__ = """Cofold a receptor with each ligand of a set using Chai-1, as a docking protocol.

    Every ligand is predicted in complex with the receptor chains. The receptor MSA is
    computed once, on the first ligand, and reused by all the others.

    Each predicted complex is superposed onto the input receptor (CA atoms), so all poses
    share its frame. Every output pose carries its own predicted receptor as proteinFile,
    and the set keeps the input receptor.

    Chai-1 does not predict binding affinity, so this protocol gives poses and confidence
    scores only. For ranking a library by predicted potency use the Boltz-2 cofolding
    protocol, which has an affinity model; Chai-1 is best used alongside it, with the
    docking consensus protocol, to see which poses two methods agree on.

    Unlike the other cofolding backends, Chai-1's command line predicts a single complex
    per call, so the model is loaded again for every ligand. That overhead is paid once
    per ligand whatever the prediction settings, which makes this protocol better suited
    to a shortlist than to a large screen.
    %s
    Reference: Chai Discovery, Chai-1 (https://github.com/chaidiscovery/chai-lab)
    """ % COFOLD_ATTRS_HELP
    _label = 'chai-1 cofolding'

    BACKEND, MOL_CLASS = 'chai', 'Chai'

    # -------------------------- DEFINE param functions ----------------------
    def _definePredictionParams(self, form):
        form.addParam('diffusionSamples', params.IntParam, default=4,
                      label='Poses per ligand: ',
                      help='Number of diffusion samples (predicted complexes) per ligand '
                           '(--num-diffn-samples).')
        form.addParam('trunkRecycles', params.IntParam, default=3, expertLevel=params.LEVEL_ADVANCED,
                      label='Trunk recycles: ',
                      help='Times the trunk feeds its output back as input to refine the pair '
                           'representation. Chai-1 default: 3.')
        form.addParam('samplingSteps', params.IntParam, default=200, expertLevel=params.LEVEL_ADVANCED,
                      label='Diffusion timesteps: ',
                      help='Denoising steps used to build each pose from noise '
                           '(--num-diffn-timesteps). Chai-1 default: 200.')
        form.addParam('trunkSamples', params.IntParam, default=1, expertLevel=params.LEVEL_ADVANCED,
                      label='Trunk samples: ',
                      help='Independent trunk passes per ligand. Above 1, Chai nests its output '
                           'in trunk_<i> subdirectories, which this protocol does not read, so '
                           'leave it at 1.')
        form.addParam('seeds', params.StringParam, default='42', expertLevel=params.LEVEL_ADVANCED,
                      label='Seeds: ',
                      help='Random seed. Several comma-separated seeds predict each ligand that many '
                           'times, multiplying the number of output poses. Chai-1 takes one seed per '
                           'call, so every extra seed is a separate run that loads the model again: '
                           'it costs a full prediction, not just another sample.')

    # --------------------------- STEPS functions ------------------------------
    def writeJob(self, ligId, smiles, rec, msa):
        """One FASTA per ligand. Chai computes its own MSA while predicting (see
        predictStep), so there is no MSA step and `msa` is empty."""
        # One record per receptor chain, in order, then the ligand: Chai names the chains
        # A, B, C... by record order, which is what the extraction step assumes.
        records = ''.join(f">protein|name=rec_{i}\n{ch['seq']}\n" for i, ch in enumerate(rec['chains']))
        with open(os.path.join(self.getInputsDir(), f'{ligId}.fasta'), 'w') as f:
            f.write(f'{records}>ligand|name={ligId}\n{smiles}\n')

    def predictStep(self):
        """One Chai call per ligand and seed: its CLI predicts a single complex at a time.

        The MSA is computed by the first call that needs it and then reused, so the
        server is queried once however many ligands and seeds there are.
        """
        rec = self.getReceptorInfo()
        gpus = self.getGpuIds()
        msaDir = os.path.join(self.getWorkDir(), 'msa')

        for ligId, seed in ((l, sd) for l in rec['ligands'] for sd in self.getSeeds()):
            outDir = os.path.join(self.getPredictionsDir(), ligId, f'seed_{seed}')
            # Chai refuses a directory that already exists, which also makes this the
            # resume check: a ligand and seed already predicted are skipped on a continued run.
            if os.path.exists(outDir):
                continue

            args = [f'"{os.path.join(self.getInputsDir(), f"{ligId}.fasta")}"', f'"{outDir}"',
                    f'--num-trunk-recycles {self.trunkRecycles.get()}',
                    f'--num-diffn-timesteps {self.samplingSteps.get()}',
                    f'--num-diffn-samples {self.diffusionSamples.get()}',
                    f'--num-trunk-samples {self.trunkSamples.get()}',
                    f'--seed {seed}']
            if gpus:
                args.append('--device cuda:0')

            if not self.useMsaServer.get():
                pass
            elif os.path.isdir(msaDir):
                args.append(f'--msa-directory "{msaDir}"')
            else:
                args.append('--use-msa-server')
                args.append(f'--msa-server-url {self.msaServerUrl.get()}')

            Plugin.runCondaCommand(self, ' '.join(args), CHAI_DIC, 'chai-lab fold',
                                   gpuIdx=','.join(gpus) if gpus else None)

            # Keep the MSA the first run generated, so the rest reuse it
            generated = os.path.join(outDir, 'msas')
            if self.useMsaServer.get() and not os.path.isdir(msaDir) and os.path.isdir(generated):
                shutil.copytree(generated, msaDir)

    # --------------------------- UTILS functions -----------------------------------
    def getSeeds(self):
        return [int(s) for s in self.seeds.get().split(',') if s.strip()]

    def getPredictionsDir(self):
        """Chai writes <this>/<ligand>/seed_<s>/, one call per ligand and seed"""
        return os.path.join(self.getWorkDir(), 'chai_results')

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = super()._validate()
        try:
            if not self.getSeeds():
                errors.append('At least one seed is required.')
        except ValueError:
            errors.append(f'Seeds must be whole numbers separated by commas: {self.seeds.get()}')
        return errors

    def _methods(self):
        return ['Protein-ligand complexes were cofolded with Chai-1.']
