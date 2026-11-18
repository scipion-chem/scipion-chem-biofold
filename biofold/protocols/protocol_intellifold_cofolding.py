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
import json
import os

import pyworkflow.protocol.params as params

from pwchem import Plugin

from biofold.constants import INTELLIFOLD_DIC
from .protocol_cofolding_base import COFOLD_ATTRS_HELP, ProtCofoldingBase


class ProtIntelliFoldCofolding(ProtCofoldingBase):
    # assigned, not a bare literal: a """...""" %% VAR expression is not a docstring
    __doc__ = """Cofold a receptor with each ligand of a set using IntelliFold, as a docking protocol.

    Every ligand is predicted in complex with the receptor chains. The receptor MSA is
    computed only once and shared by all ligand jobs, and all jobs run in a single
    IntelliFold call, so the model is loaded once.

    Each predicted complex is superposed onto the input receptor (CA atoms), so all poses
    share its frame. Every output pose carries its own predicted receptor as proteinFile,
    and the set keeps the input receptor.

    IntelliFold does not predict binding affinity, so this protocol gives poses and
    confidence scores only. For ranking a library by predicted potency use the Boltz-2
    cofolding protocol, which has an affinity model; IntelliFold is best used alongside
    it, with the docking consensus protocol, to see which poses two methods agree on.
    %s
    Note on the pLDDT columns: IntelliFold computes pLDDT once per ligand and writes the
    same values to every diffusion sample, so _cofoldLigPLDDT and _cofoldPLDDT are
    identical across the poses of one ligand and cannot be used to tell them apart. Rank
    the poses of a ligand by _cofoldConfidence or _cofoldLigIpTM, which do vary. (Boltz-2
    does report a per-pose pLDDT.)

    Reference: IntelliGen AI, IntelliFold (https://github.com/IntelliGen-AI/IntelliFold)
    """ % COFOLD_ATTRS_HELP
    _label = 'intellifold cofolding'

    BACKEND, MOL_CLASS = 'intellifold', 'IntelliFold'
    MSA_ENV = INTELLIFOLD_DIC

    # -------------------------- DEFINE param functions ----------------------
    def _definePredictionParams(self, form):
        form.addParam('diffusionSamples', params.IntParam, default=4,
                      label='Poses per ligand: ',
                      help='Number of diffusion samples (predicted complexes) per ligand.')
        form.addParam('model', params.EnumParam, default=2, choices=['v1', 'v2', 'v2-flash'],
                      label='Model: ', display=params.EnumParam.DISPLAY_HLIST,
                      help='IntelliFold weights to use. v2-flash is the fastest and the '
                           'IntelliFold default.')
        form.addParam('recyclingSteps', params.IntParam, default=10, expertLevel=params.LEVEL_ADVANCED,
                      label='Recycling steps: ',
                      help='Times the trunk feeds its output back as input to refine the pair '
                           'representation before the structure is generated. IntelliFold '
                           'default: 10.')
        form.addParam('samplingSteps', params.IntParam, default=200, expertLevel=params.LEVEL_ADVANCED,
                      label='Sampling steps: ',
                      help='Denoising steps of the diffusion module used to build each pose from '
                           'noise. Fewer steps are faster but give rougher structures. '
                           'IntelliFold default: 200.')
        form.addParam('seed', params.StringParam, default='42', expertLevel=params.LEVEL_ADVANCED,
                      label='Seed: ',
                      help='Random seed. Several comma-separated seeds predict each ligand that '
                           'many times, multiplying the number of output poses.')

    # --------------------------- STEPS functions ------------------------------
    def writeJob(self, ligId, smiles, rec, msa):
        """One job per ligand, in the Boltz-1 schema. It has no affinity block because
        IntelliFold does not predict an affinity."""
        proteins = [{'protein': {'id': ids if len(ids) > 1 else ids[0], 'sequence': seq,
                                 'msa': msa.get(seq, {}).get('msa', 'empty')}}
                    for seq, ids in self.getChainsBySequence(rec).items()]
        job = {'version': 1,
               'sequences': proteins + [{'ligand': {'id': rec['ligandChain'], 'smiles': smiles}}]}
        # IntelliFold reads the file with a YAML parser, and JSON is valid YAML: writing it
        # this way needs no yaml module, which scipion3 does not ship.
        with open(os.path.join(self.getInputsDir(), f'{ligId}.yaml'), 'w') as f:
            json.dump(job, f, indent=2)

    def predictStep(self):
        """All ligand jobs in one IntelliFold call, so the model is loaded once"""
        gpus = self.getGpuIds()
        args = [f'"{self.getInputsDir()}"',
                f'--out_dir "{self.getWorkDir()}"',
                f'--recycling_iters {self.recyclingSteps.get()}',
                f'--sampling_steps {self.samplingSteps.get()}',
                f'--num_diffusion_samples {self.diffusionSamples.get()}',
                f'--model {self.getEnumText("model")}',
                f'--seed {self.seed.get()}',
                f'--num_workers {self.numberOfThreads.get()}',
                '--output_format mmcif']
        # The MSA is already written into every YAML, so the server must not be called again
        Plugin.runCondaCommand(self, ' '.join(args), INTELLIFOLD_DIC, 'intellifold predict',
                               gpuIdx=','.join(gpus) if gpus else None)

    # --------------------------- UTILS functions -----------------------------------
    def getPredictionsDir(self):
        """IntelliFold nests its output under the stem of the input directory"""
        return os.path.join(self.getWorkDir(), os.path.basename(self.getInputsDir()), 'predictions')

    # --------------------------- INFO functions -----------------------------------
    def _methods(self):
        return ['Protein-ligand complexes were cofolded with IntelliFold.']
