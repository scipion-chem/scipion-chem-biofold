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
import glob
import json
import os
import shutil

import pyworkflow.protocol.params as params

from pwchem import Plugin

from biofold.constants import ALPHAFOLD3_DIC
from .protocol_cofolding_base import COFOLD_ATTRS_HELP, SCRIPTS_DIR, ProtCofoldingBase


class ProtAlphaFold3Cofolding(ProtCofoldingBase):
    # assigned, not a bare literal: a """...""" % VAR expression is not a docstring
    __doc__ = """Cofold a receptor with each ligand of a set using AlphaFold 3, as a docking protocol.

    Every ligand is predicted in complex with the receptor chains. The receptor MSA is
    computed only once, on the ColabFold server, and shared by all ligand jobs, and all jobs
    run in a single AlphaFold 3 call, so the model is loaded and compiled once per input size.

    AlphaFold 3's own database search is not used: it needs about 630 GB of genetic
    databases. The MSA comes from the server instead, and no structural templates are used.
    Without the server the run is single-sequence, which is faster and usually less accurate.

    Each predicted complex is superposed onto the input receptor (CA atoms), so all poses
    share its frame. Every output pose carries its own predicted receptor as proteinFile,
    and the set keeps the input receptor.

    AlphaFold 3 does not predict binding affinity, so this protocol gives poses and
    confidence scores only. For ranking a library by predicted potency use the Boltz-2
    cofolding protocol, which has an affinity model; AlphaFold 3 is best used alongside it,
    with the docking consensus protocol, to see which poses two methods agree on.
    %s
    AlphaFold 3 runs on a single GPU: with several selected, the first one is used.

    This protocol does not install AlphaFold 3. Point ALPHAFOLD3_HOME at a checkout of
    its repository (the folder that holds run_alphafold.py), and keep its conda environment
    named alphafold3. The model parameters are obtained from Google DeepMind and their
    use, like that of the predictions, is subject to its terms of use.

    Reference: Abramson et al., Nature 2024 (AlphaFold 3)
    """ % COFOLD_ATTRS_HELP
    _label = 'alphafold3 cofolding'

    BACKEND, MOL_CLASS = 'alphafold3', 'AlphaFold3'
    MSA_ENV = ALPHAFOLD3_DIC

    # -------------------------- DEFINE param functions ----------------------
    def _definePredictionParams(self, form):
        form.addParam('diffusionSamples', params.IntParam, default=4,
                      label='Poses per ligand: ',
                      help='Number of diffusion samples (predicted complexes) per ligand and seed '
                           '(--num_diffusion_samples). AlphaFold 3 default: 5.')
        form.addParam('recyclingSteps', params.IntParam, default=10, expertLevel=params.LEVEL_ADVANCED,
                      label='Recycling steps: ',
                      help='Times the trunk feeds its output back as input to refine the pair '
                           'representation before the structure is generated (--num_recycles). '
                           'AlphaFold 3 default: 10.')
        form.addParam('seeds', params.StringParam, default='42', expertLevel=params.LEVEL_ADVANCED,
                      label='Seeds: ',
                      help='Random seed. Several comma-separated seeds predict each ligand that '
                           'many times, multiplying the number of output poses.')
        form.addParam('unifiedMemory', params.BooleanParam, default=True, expertLevel=params.LEVEL_ADVANCED,
                      label='Unified memory: ',
                      help='Let AlphaFold 3 spill GPU memory to host RAM instead of failing with an '
                           'out-of-memory error. Recommended for GPUs below 80 GB (AlphaFold 3 was '
                           'tuned on A100/H100 80 GB); it only slows down runs that would not fit '
                           'otherwise.')
        form.addParam('modelDir', params.StringParam, default='', expertLevel=params.LEVEL_ADVANCED,
                      label='Model parameters directory: ',
                      help='Folder with the AlphaFold 3 parameters (af3.bin or af3.bin.zst, '
                           '--model_dir). Leave empty to look for a "models" folder inside '
                           'ALPHAFOLD3_HOME, then beside it.')

    # --------------------------- STEPS functions ------------------------------
    def writeJob(self, ligId, smiles, rec, msa):
        """One AlphaFold 3 input JSON per ligand. The job name is the ligand id, which is also
        the directory AlphaFold 3 writes that ligand's predictions into."""
        proteins = [{'protein': dict({'id': ids if len(ids) > 1 else ids[0], 'sequence': seq,
                                      'templates': []},
                                     # no MSA means an explicit single-sequence run: left unset,
                                     # AlphaFold 3 would try to build one from the databases
                                     **msa.get(seq, {'unpairedMsa': '', 'pairedMsa': ''}))}
                    for seq, ids in self.getChainsBySequence(rec).items()]
        job = {'name': ligId, 'modelSeeds': self.getSeeds(),
               'sequences': proteins + [{'ligand': {'id': rec['ligandChain'], 'smiles': smiles}}],
               'dialect': 'alphafold3',
               # version 2 is the first that takes MSAs from a path, kept out of the JSON
               'version': 2}
        with open(os.path.join(self.getInputsDir(), f'{ligId}.json'), 'w') as f:
            json.dump(job, f, indent=1)

    def predictStep(self):
        """All the ligands still to do in one AlphaFold 3 call, so the model is compiled once.

        A continued run only predicts the ligands that have no final summary yet.
        """
        rec = self.getReceptorInfo()
        todo = [ligId for ligId in rec['ligands'] if not os.path.exists(self.getSummaryFile(ligId))]
        if not todo:
            return

        # --input_dir takes every JSON of a folder, so the pending ones are gathered in one
        todoDir = os.path.join(self.getWorkDir(), 'inputs_todo')
        shutil.rmtree(todoDir, ignore_errors=True)
        os.makedirs(todoDir)
        for ligId in todo:
            shutil.copy(os.path.join(self.getInputsDir(), f'{ligId}.json'), todoDir)

        args = [f'"{self.getRunScript()}"',
                f'--input_dir "{todoDir}"',
                f'--output_dir "{self.getPredictionsDir()}"',
                f'--model_dir "{self.getModelDir()}"',
                # the MSA is already in every JSON, so no database search is run
                '--norun_data_pipeline',
                f'--num_recycles {self.recyclingSteps.get()}',
                f'--num_diffusion_samples {self.diffusionSamples.get()}',
                # an interrupted ligand leaves a non-empty folder, which would otherwise be
                # answered with a new time-stamped one
                '--force_output_dir',
                f'--jax_compilation_cache_dir "{os.path.join(self.getWorkDir(), "jax_cache")}"']

        # Through env, not export: CUDA_VISIBLE_DEVICES goes in front of the program, and a
        # variable assignment before "export" never reaches what runs after it
        envVars = ['XLA_FLAGS=--xla_gpu_enable_triton_gemm=false']
        if self.unifiedMemory.get():
            envVars += ['XLA_PYTHON_CLIENT_PREALLOCATE=false', 'TF_FORCE_UNIFIED_MEMORY=true',
                        'XLA_CLIENT_MEM_FRACTION=3.2']

        gpus = self.getGpuIds()
        # one process drives one GPU
        gpuStr = f'CUDA_VISIBLE_DEVICES={gpus[0]}' if gpus else ''
        self.runJob(f'{self.getEnvActivation()} && {gpuStr} env {" ".join(envVars)} python ',
                    ' '.join(args), env=Plugin.getEnviron(), numberOfThreads=1)

    def runMsaScript(self, args):
        """cofoldMsa.py needs only the standard library, but runs in the AlphaFold 3 env"""
        self.runJob(f'{self.getEnvActivation()} && python {os.path.join(SCRIPTS_DIR, "cofoldMsa.py")}',
                    args, env=Plugin.getEnviron())

    # --------------------------- UTILS functions -----------------------------------
    def getEnvActivation(self):
        """Activation line for the environment, by its exact name. pwchem's helpers build it as
        name-version, which is not how an AlphaFold 3 checkout is usually set up."""
        return f"{Plugin.getCondaActivationCmd()}conda activate {ALPHAFOLD3_DIC['env']}"

    def getAlphaFold3Home(self):
        home = Plugin.getVar(ALPHAFOLD3_DIC['home'])
        return os.path.abspath(home) if home else ''

    def getRunScript(self):
        return os.path.join(self.getAlphaFold3Home(), 'run_alphafold.py')

    def getModelDir(self):
        """The parameters folder: the given one, else a "models" folder in or beside the checkout"""
        given = (self.modelDir.get() or '').strip()
        if given:
            return os.path.abspath(given)
        home = self.getAlphaFold3Home()
        for cand in (os.path.join(home, 'models'), os.path.join(os.path.dirname(home), 'models')):
            if home and glob.glob(os.path.join(cand, 'af3.bin*')):
                return cand
        return ''

    def getSeeds(self):
        return [int(s) for s in self.seeds.get().split(',') if s.strip()]

    def getPredictionsDir(self):
        """AlphaFold 3 writes <this>/<job name>/seed-<s>_sample-<k>/, one job per ligand"""
        return os.path.join(self.getWorkDir(), 'alphafold3_results')

    def getSummaryFile(self, ligId):
        """Written once a job is complete: the marker of a ligand already predicted"""
        return os.path.join(self.getPredictionsDir(), ligId, f'{ligId}_summary_confidences.json')

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = super()._validate()
        if not os.path.isfile(self.getRunScript()):
            errors.append('run_alphafold.py not found: set ALPHAFOLD3_HOME to the AlphaFold 3 repository.')
        if not self.getModelDir():
            errors.append('AlphaFold 3 model parameters (af3.bin or af3.bin.zst) not found: set '
                          '"Model parameters directory".')
        try:
            self.getSeeds()
        except ValueError:
            errors.append(f'Seeds must be whole numbers separated by commas: {self.seeds.get()}')
        return errors

    def _methods(self):
        return ['Protein-ligand complexes were cofolded with AlphaFold 3 (Abramson et al. 2024).']
