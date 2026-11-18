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

from biofold.constants import PROTENIX_DIC
from .protocol_cofolding_base import COFOLD_ATTRS_HELP, ProtCofoldingBase
from .protocol_protenix import ProtProtenix

MODELS = ProtProtenix.models


class ProtProtenixCofolding(ProtCofoldingBase):
    # assigned, not a bare literal: a """...""" % VAR expression is not a docstring
    __doc__ = """Cofold a receptor with each ligand of a set using Protenix, as a docking protocol.

    Every ligand is predicted in complex with the receptor chains. The receptor MSA is
    computed only once and shared by all ligand jobs, and all jobs run in a single
    Protenix call, so the model is loaded once.

    Each predicted complex is superposed onto the input receptor (CA atoms), so all poses
    share its frame. Every output pose carries its own predicted receptor as proteinFile,
    and the set keeps the input receptor.

    Protenix does not predict binding affinity, so this protocol gives poses and
    confidence scores only. For ranking a library by predicted potency use the Boltz-2
    cofolding protocol, which has an affinity model; Protenix is best used alongside it,
    with the docking consensus protocol, to see which poses two methods agree on.
    %s
    Reference: ByteDance, Protenix (https://github.com/bytedance/Protenix)
    """ % COFOLD_ATTRS_HELP
    _label = 'protenix cofolding'

    BACKEND, MOL_CLASS = 'protenix', 'Protenix'
    MSA_ENV = PROTENIX_DIC

    # -------------------------- DEFINE param functions ----------------------
    def _definePredictionParams(self, form):
        form.addParam('diffusionSamples', params.IntParam, default=4,
                      label='Poses per ligand: ',
                      help='Number of diffusion samples (predicted complexes) per ligand.')
        form.addParam('model', params.EnumParam, default=0, choices=MODELS,
                      label='Model: ',
                      help='Protenix weights to use:\n\n'
                           '- *protenix-v2*: the newest model (464M parameters), and the default.\n'
                           '- *protenix_base_default_v1.0.0*: trained with a data cutoff aligned '
                           'with AlphaFold3 (2021-09-30). Recommended for fair benchmarking and '
                           'comparison with other methods.\n'
                           '- *protenix_base_20250630_v1.0.0*: trained with a more recent data '
                           'cutoff (2025-06-30).\n'
                           '- *protenix_base_default_v0.5.0*: previous version, kept for '
                           'backward compatibility.\n\n'
                           'The weights of a model are downloaded the first time it is used.')
        form.addParam('recyclingSteps', params.IntParam, default=10, expertLevel=params.LEVEL_ADVANCED,
                      label='Recycling steps: ',
                      help='Times the trunk feeds its output back as input to refine the pair '
                           'representation before the structure is generated. Protenix '
                           'default: 10.')
        form.addParam('samplingSteps', params.IntParam, default=200, expertLevel=params.LEVEL_ADVANCED,
                      label='Sampling steps: ',
                      help='Denoising steps of the diffusion module used to build each pose from '
                           'noise. Fewer steps are faster but give rougher structures. '
                           'Protenix default: 200.')
        form.addParam('seeds', params.StringParam, default='101', expertLevel=params.LEVEL_ADVANCED,
                      label='Seeds: ',
                      help='Random seed. Several comma-separated seeds predict each ligand that '
                           'many times, multiplying the number of output poses.')

    def _defineExtraMsaParams(self, group):
        group.addParam('msaServerMode', params.EnumParam, default=0,
                       choices=['protenix', 'colabfold'], condition='useMsaServer',
                       label='MSA server: ', display=params.EnumParam.DISPLAY_HLIST,
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Which MMseqs2 server Protenix queries. The "MSA server URL" above '
                            'is not used by Protenix, which picks the endpoint from this choice.')

    # --------------------------- STEPS functions ------------------------------
    def getMsaArgs(self):
        """Protenix picks its MSA endpoint from a server mode rather than a URL"""
        return f'{self.getEnumText("msaServerMode")} "{self.getProtenixRoot()}"'

    def writeJob(self, ligId, smiles, rec, msa):
        """One Protenix job per ligand. The job name is the ligand id, which is also the
        directory Protenix writes that ligand's predictions into."""
        proteins = [{'proteinChain': dict({'sequence': seq, 'count': len(ids), 'id': ids},
                                          **msa.get(seq, {}))}
                    for seq, ids in self.getChainsBySequence(rec).items()]
        job = {'name': ligId,
               'sequences': proteins + [{'ligand': {'ligand': smiles, 'count': 1,
                                                    'id': [rec['ligandChain']]}}]}
        # Protenix takes a list of jobs per file
        with open(os.path.join(self.getInputsDir(), f'{ligId}.json'), 'w') as f:
            json.dump([job], f, indent=1)

    def predictStep(self):
        """All ligand jobs in one Protenix call, so the model is loaded once"""
        args = [f'-i "{self.getInputsDir()}"',
                f'-o "{self.getPredictionsDir()}"',
                # the MSA is already in every job JSON, so no further search is run
                '--use_msa TRUE',
                f'--cycle {self.recyclingSteps.get()}',
                f'--step {self.samplingSteps.get()}',
                f'--sample {self.diffusionSamples.get()}',
                f'--seeds {self.seeds.get()}',
                f'--model_name {self.getEnumText("model")}']

        # Protenix JIT-compiles a fused LayerNorm at runtime, so it needs a CUDA toolchain
        prefix = ("export CUDA_HOME=$CONDA_PREFIX && "
                  "export CPATH=$CONDA_PREFIX/targets/x86_64-linux/include:$CPATH && "
                  "export LD_LIBRARY_PATH=$CONDA_PREFIX/targets/x86_64-linux/lib:$LD_LIBRARY_PATH && "
                  f"export PROTENIX_ROOT_DIR={self.getProtenixRoot()} && protenix pred")

        gpus = self.getGpuIds()
        Plugin.runCondaCommand(self, ' '.join(args), PROTENIX_DIC, prefix,
                               gpuIdx=','.join(gpus) if gpus else None,
                               cwd=self.getProtenixRoot())

    # --------------------------- UTILS functions -----------------------------------
    def getProtenixRoot(self):
        return os.path.abspath(Plugin.getVar(PROTENIX_DIC['home']))

    def getPredictionsDir(self):
        """Protenix writes <this>/<job name>/seed_<s>/predictions/, one job per ligand"""
        return os.path.join(self.getWorkDir(), 'protenix_results')

    # --------------------------- INFO functions -----------------------------------
    def _methods(self):
        return ['Protein-ligand complexes were cofolded with Protenix.']
