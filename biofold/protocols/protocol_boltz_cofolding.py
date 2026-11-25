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

from biofold.constants import BOLTZ_DIC
from .protocol_cofolding_base import COFOLD_ATTRS, COFOLD_ATTRS_HELP, ProtCofoldingBase


class ProtBoltzCofolding(ProtCofoldingBase):
    # assigned, not a bare literal: a """...""" %% VAR expression is not a docstring
    __doc__ = """Cofold a receptor with each ligand of a set using Boltz-2, as a docking protocol.

    Every ligand is predicted in complex with the receptor chains, and Boltz-2's affinity
    module scores the binding. The receptor MSA is computed only once and shared by all
    ligand jobs, and all jobs run in a single Boltz call, so the model is loaded once.

    Each predicted complex is superposed onto the input receptor (CA atoms), so all poses
    share its frame. Every output pose carries its own predicted receptor as proteinFile,
    and the set keeps the input receptor.

    Boltz-2 is the only cofolding backend that predicts binding affinity, which is what
    makes it usable for ranking a library rather than only for posing a known binder.
    %s
    Boltz-2 adds two columns of its own:

    - *_boltzAffinity*: predicted log10(IC50) in uM, lower binds stronger. The affinity
      model takes the top-ranked pose, crops the pocket around the ligand, re-samples it
      and predicts from its best sample, so there is one value per ligand, shared by all
      of its poses. Boltz's own docs advise using it to compare active molecules with one
      another, not to tell binders from inactives.
    - *_boltzIntProbability*: probability that the ligand is a binder. The value to rank
      by when screening.

    Affinity is only supported for ligands up to 128 heavy atoms; it was trained on
    ligands up to 56.

    Reference: Passaro et al., bioRxiv 2025 (Boltz-2)
    """ % COFOLD_ATTRS_HELP
    _label = 'boltz-2 cofolding'

    BACKEND, MOL_CLASS = 'boltz', 'Boltz'
    MSA_ENV = BOLTZ_DIC
    SCORE_ATTRS = dict(COFOLD_ATTRS,
                       affinity='_boltzAffinity',
                       bindingProbability='_boltzIntProbability')

    # -------------------------- DEFINE param functions ----------------------
    def _defineExtraMsaParams(self, group):
        group.addParam('maxMsaSeqs', params.IntParam, default=8192,
                       label='Max MSA sequences: ', expertLevel=params.LEVEL_ADVANCED,
                       help='Lower it (e.g. 512-1024) to reduce GPU memory use.')

    def _definePredictionParams(self, form):
        form.addParam('diffusionSamples', params.IntParam, default=4,
                      label='Poses per ligand: ',
                      help='Number of diffusion samples (predicted complexes) per ligand.')
        form.addParam('infPot', params.BooleanParam, default=True,
                      label="Inference potentials: ",
                      help='Use steering potentials to improve the physical plausibility of the '
                           'poses. This runs several steered trajectories per pose, so it is one '
                           'of the most expensive options: turn it off for screening.')
        form.addParam('recyclingSteps', params.IntParam, default=3, expertLevel=params.LEVEL_ADVANCED,
                      label='Recycling steps: ',
                      help='Times the trunk (pairformer) feeds its output back as input to refine '
                           'the pair representation before the structure is generated. More '
                           'steps can improve hard cases; each extra step adds roughly one '
                           'trunk pass of compute.')
        form.addParam('samplingSteps', params.IntParam, default=200, expertLevel=params.LEVEL_ADVANCED,
                      label='Sampling steps: ',
                      help='Denoising steps of the diffusion module used to build each pose from '
                           'noise. Fewer steps are faster but give rougher, less accurate '
                           'structures. Boltz default: 200.')
        form.addParam('seed', params.IntParam, default=42, expertLevel=params.LEVEL_ADVANCED,
                      label='Seed: ',
                      help='Random seed, so a run can be reproduced. Boltz takes a single seed and '
                           'cannot run several in one call.')
        form.addParam('stepScale', params.FloatParam, default=1.638, expertLevel=params.LEVEL_ADVANCED,
                      label='Step scale: ',
                      help='Diffusion temperature: lower values give more diverse samples.')
        form.addParam('maxParallelSamples', params.IntParam, default=5, expertLevel=params.LEVEL_ADVANCED,
                      label='Max parallel samples: ',
                      help='Samples predicted at once. Lower it to reduce GPU memory use.')
        form.addParam('affinityMWcorr', params.BooleanParam, default=False,
                      label="Molecular weight correction: ", expertLevel=params.LEVEL_ADVANCED,
                      help='Apply the molecular weight correction to the affinity prediction.')
        form.addParam('diffusionSamplesAff', params.IntParam, default=5, expertLevel=params.LEVEL_ADVANCED,
                      label='Diffusion samples for affinity: ',
                      help='The affinity is predicted by a separate Boltz-2 affinity model. It '
                           'takes the top-ranked pose, crops it to the pocket around the ligand '
                           '(up to 256 residues/ligand tokens), runs its trunk once and then '
                           're-samples the cropped complex this many times. The sample with the '
                           'best ipTM is used to predict the affinity; the samples are not written '
                           'out. One affinity value per ligand, whatever the number of output '
                           'poses. More samples cost time linearly (they run one by one).')
        form.addParam('samplingStepsAff', params.IntParam, default=200, expertLevel=params.LEVEL_ADVANCED,
                      label='Sampling steps for affinity: ',
                      help='Denoising steps the affinity model uses for each of its own samples '
                           '(see "Diffusion samples for affinity"). Lowering it speeds up the '
                           'affinity stage, at the cost of a rougher pose to predict the affinity '
                           'from.')

    # --------------------------- STEPS functions ------------------------------
    def writeJob(self, ligId, smiles, rec, msa):
        """One Boltz job per ligand, with the affinity of that ligand requested"""
        proteins = [{'protein': {'id': ids if len(ids) > 1 else ids[0], 'sequence': seq,
                                 'msa': msa.get(seq, {}).get('msa', 'empty')}}
                    for seq, ids in self.getChainsBySequence(rec).items()]
        job = {'version': 1,
               'sequences': proteins + [{'ligand': {'id': rec['ligandChain'], 'smiles': smiles}}],
               'properties': [{'affinity': {'binder': rec['ligandChain']}}]}
        # Boltz reads the file with a YAML parser, and JSON is valid YAML: writing it this
        # way needs no yaml module, which scipion3 does not ship.
        with open(os.path.join(self.getInputsDir(), f'{ligId}.yaml'), 'w') as f:
            json.dump(job, f, indent=2)

    def predictStep(self):
        """All ligand jobs in one Boltz call. Finished jobs are skipped on continue."""
        gpus = self.getGpuIds()
        args = [f'"{self.getInputsDir()}"',
                f'--out_dir "{self.getWorkDir()}"',
                f'--cache "{os.path.abspath(os.path.join(Plugin.getVar(BOLTZ_DIC["home"]), "mol"))}"',
                f'--recycling_steps {self.recyclingSteps.get()}',
                f'--sampling_steps {self.samplingSteps.get()}',
                f'--diffusion_samples {self.diffusionSamples.get()}',
                f'--max_parallel_samples {self.maxParallelSamples.get()}',
                f'--step_scale {self.stepScale.get()}',
                f'--seed {self.seed.get()}',
                f'--max_msa_seqs {self.maxMsaSeqs.get()}',
                f'--diffusion_samples_affinity {self.diffusionSamplesAff.get()}',
                f'--sampling_steps_affinity {self.samplingStepsAff.get()}',
                f'--num_workers {self.numberOfThreads.get()}']
        if self.infPot.get():
            args.append('--use_potentials')
        if self.affinityMWcorr.get():
            args.append('--affinity_mw_correction')
        args.append(f'--accelerator gpu --devices {len(gpus)}' if gpus else '--accelerator cpu')

        Plugin.runCondaCommand(self, ' '.join(args), BOLTZ_DIC, 'boltz predict',
                               gpuIdx=','.join(gpus) if gpus else None)
        # Boltz names its output boltz_results_<input dir>. Renamed only once it succeeded,
        # so a continued run still finds (and skips) the predictions already done
        os.replace(os.path.join(self.getWorkDir(), 'boltz_results_inputs'), self.getResultsDir())

    # --------------------------- UTILS functions -----------------------------------
    def getResultsDir(self):
        return os.path.join(self.getWorkDir(), 'boltz_results')

    def getPredictionsDir(self):
        return os.path.join(self.getResultsDir(), 'predictions')

    # --------------------------- INFO functions -----------------------------------
    def _methods(self):
        return ['Protein-ligand complexes were cofolded and their binding affinity predicted '
                'with Boltz-2 (Passaro et al. 2025).']
