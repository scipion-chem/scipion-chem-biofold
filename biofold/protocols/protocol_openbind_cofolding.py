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

import pyworkflow.protocol.params as params

from pwchem import Plugin

from biofold.constants import OPENBIND_DIC
from .protocol_cofolding_base import COFOLD_ATTRS_HELP, ProtCofoldingBase


class ProtOpenBindCofolding(ProtCofoldingBase):
    # assigned, not a bare literal: a """...""" % VAR expression is not a docstring
    __doc__ = """Cofold a receptor with each ligand of a set using OpenBind, as a docking protocol.

    OpenBind is the OpenFold3 model run with the OpenBind weights, through its
    `run_openfold predict` command.

    Every ligand is predicted in complex with the receptor chains. The receptor MSA is
    computed only once and shared by all ligand queries, and all queries run in a single
    OpenFold3 call, so the model is loaded once.

    Each predicted complex is superposed onto the input receptor (CA atoms), so all poses
    share its frame. Every output pose carries its own predicted receptor as proteinFile,
    and the set keeps the input receptor.

    OpenBind does not predict binding affinity, so this protocol gives poses and
    confidence scores only. For ranking a library by predicted potency use the Boltz-2
    cofolding protocol, which has an affinity model; OpenBind is best used alongside it,
    with the docking consensus protocol, to see which poses two methods agree on.
    %s
    Templates are turned off, so a run depends only on the sequences, the ligand and the
    MSA.

    This protocol does not install OpenBind. Point OPENBIND_HOME at an existing install
    and keep its conda environment named as the plugin expects (openbind-0).

    Reference: AlQuraishi Laboratory, OpenFold3 (https://github.com/aqlaboratory/openfold3)
    """ % COFOLD_ATTRS_HELP
    _label = 'openbind cofolding'

    BACKEND, MOL_CLASS = 'openbind', 'OpenBind'
    # The MSA client is a command of its own (run_openfold align-msa-server), not a python
    # API, so the shared cofoldMsa.py script is not used: see runPrepare below.
    MSA_ENV = None

    # -------------------------- DEFINE param functions ----------------------
    def _definePredictionParams(self, form):
        form.addParam('diffusionSamples', params.IntParam, default=4,
                      label='Poses per ligand: ',
                      help='Number of diffusion samples (predicted complexes) per ligand '
                           '(--num-diffusion-samples). OpenFold3 default: 5.')
        form.addParam('seeds', params.StringParam, default='42', expertLevel=params.LEVEL_ADVANCED,
                      label='Seeds: ',
                      help='Random seed. Several comma-separated seeds predict each ligand that '
                           'many times, multiplying the number of output poses.')
        form.addParam('ckptFile', params.StringParam, default='', expertLevel=params.LEVEL_ADVANCED,
                      label='Checkpoint file: ',
                      help='Model weights (--inference-ckpt-path). Leave empty to use the only '
                           '.pt file found in OPENBIND_HOME.')

    # --------------------------- STEPS functions ------------------------------
    def runPrepare(self):
        """Receptor MSA once, then every ligand query in a single file.

        OpenFold3 takes all of its queries in one JSON file, so there is no per-ligand job
        file to write and the writeJob hook of the base class is not used.
        """
        rec = self.getReceptorInfo()
        msa = self.computeMsa(rec) if self.useMsaServer.get() else {}

        queries = {}
        for ligId, smiles in rec['ligands'].items():
            chains = [dict({'molecule_type': 'protein', 'chain_ids': ids, 'sequence': seq},
                           **msa.get(seq, {}))
                      for seq, ids in self.getChainsBySequence(rec).items()]
            chains.append({'molecule_type': 'ligand', 'chain_ids': [rec['ligandChain']],
                           'smiles': smiles})
            # Without an MSA, OpenFold3 falls back to a single-sequence one
            queries[ligId] = {'chains': chains, 'use_msas': bool(msa)}

        os.makedirs(self.getInputsDir(), exist_ok=True)
        with open(self.getQueryFile(), 'w') as f:
            json.dump({'queries': queries}, f, indent=1)

    def computeMsa(self, rec):
        """Query the MSA server once, for the receptor chains alone.

        Returns {receptor sequence: {the query keys that point at its alignments}}, ready
        to paste into every ligand query. Identical chains share one sequence, so one MSA.
        """
        # the command writes query_msa.json straight into this directory, without creating it
        msaDir = os.path.join(self.getWorkDir(), 'msa')
        os.makedirs(msaDir, exist_ok=True)
        recQuery = os.path.join(self.getWorkDir(), 'receptor_query.json')
        with open(recQuery, 'w') as f:
            json.dump({'queries': {'receptor': {'chains': [
                {'molecule_type': 'protein', 'chain_ids': ids, 'sequence': seq}
                for seq, ids in self.getChainsBySequence(rec).items()]}}}, f, indent=1)

        # JSON is valid YAML, so writing the settings this way needs no yaml module
        settings = os.path.join(self.getWorkDir(), 'msaSettings.yml')
        with open(settings, 'w') as f:
            json.dump({'server_url': self.msaServerUrl.get()}, f)

        args = f'--query-json "{recQuery}" --output-dir "{msaDir}" ' \
               f'--msa-computation-settings-yaml "{settings}"'
        Plugin.runCondaCommand(self, args, OPENBIND_DIC,
                               f'{self.getCacheExport()}run_openfold align-msa-server')

        # The command writes back the same query with the alignment paths filled in
        with open(os.path.join(msaDir, 'query_msa.json')) as f:
            updated = json.load(f)
        msaKeys = ('main_msa_file_paths', 'paired_msa_file_paths')
        return {ch['sequence']: {k: ch[k] for k in msaKeys if ch.get(k)}
                for query in updated['queries'].values() for ch in query['chains']
                if ch.get('sequence')}

    def predictStep(self):
        """Every ligand query in one OpenFold3 call, so the model is loaded once"""
        gpus = self.getGpuIds()
        args = [f'--query-json "{self.getQueryFile()}"',
                f'--output-dir "{self.getPredictionsDir()}"',
                f'--inference-ckpt-path "{self.getCheckpoint()}"',
                f'--num-diffusion-samples {self.diffusionSamples.get()}',
                # the alignments are already in the query file, so the server is not called
                '--use-msa-server False',
                '--use-templates False',
                f'--runner-yaml "{self.writeRunnerYaml(gpus)}"']

        Plugin.runCondaCommand(self, ' '.join(args), OPENBIND_DIC,
                               f'{self.getCacheExport()}run_openfold predict',
                               gpuIdx=','.join(gpus) if gpus else None)

    def getCacheExport(self):
        """Tie the run to this install: OpenFold3 reads its cache for the default
        checkpoint and for a runner.yml of its own, which ours then overrides. Without it
        the cache is the user's ~/.openfold3, whose settings would leak into the run.

        Through env, not export: runCondaCommand puts CUDA_VISIBLE_DEVICES in front of the
        program, and a variable assignment before "export" never reaches what runs after it."""
        home = Plugin.getVar(OPENBIND_DIC['home'])
        return f'env OPENFOLD_CACHE="{home}" ' if home else ''

    def writeRunnerYaml(self, gpus):
        """The settings OpenFold3 takes only from its runner YAML, not from the command line.

        Written as JSON, which is valid YAML, so no yaml module is needed.
        """
        config = {
            'experiment_settings': {'seeds': self.getSeeds(),
                                    # a continued run keeps the predictions already written
                                    'skip_existing': True},
            'pl_trainer_args': {'accelerator': 'gpu' if gpus else 'cpu',
                                'devices': max(len(gpus), 1)},
            # cif keeps the per-atom pLDDT in the B-factor column, which is where the
            # ligand pLDDT is read from; the per-atom PAE/PDE files are not used here
            'output_writer_settings': {'structure_format': 'cif',
                                       'write_full_confidence_scores': False}}

        runnerYaml = os.path.join(self.getWorkDir(), 'runner.yml')
        with open(runnerYaml, 'w') as f:
            json.dump(config, f, indent=1)
        return runnerYaml

    # --------------------------- UTILS functions -----------------------------------
    def getQueryFile(self):
        return os.path.join(self.getInputsDir(), 'queries.json')

    def getSeeds(self):
        return [int(s) for s in self.seeds.get().split(',') if s.strip()]

    def getCheckpoint(self):
        """The weights to predict with: the given file, or the only .pt in OPENBIND_HOME"""
        ckpt = (self.ckptFile.get() or '').strip()
        if ckpt:
            return os.path.abspath(ckpt)
        home = Plugin.getVar(OPENBIND_DIC['home'])
        found = sorted(glob.glob(os.path.join(home, '*.pt'))) if home else []
        return found[0] if len(found) == 1 else ''

    def getPredictionsDir(self):
        """OpenFold3 writes <this>/<query name>/seed_<s>/, one query per ligand"""
        return os.path.join(self.getWorkDir(), 'openbind_results')

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = super()._validate()
        if not self.getCheckpoint():
            errors.append('No model weights found: set "Checkpoint file", or leave exactly one '
                          '.pt file in OPENBIND_HOME.')
        try:
            self.getSeeds()
        except ValueError:
            errors.append(f'Seeds must be whole numbers separated by commas: {self.seeds.get()}')
        return errors

    def _methods(self):
        return ['Protein-ligand complexes were cofolded with OpenBind (OpenFold3).']
