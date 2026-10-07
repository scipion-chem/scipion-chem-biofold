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
import re

import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
from pyworkflow.protocol.params import USE_GPU, GPU_LIST
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import convertToSdf, getBaseFileName, getChainIds

from biofold.constants import BOLTZ_DIC

FROM_STRUCTURE, FROM_SEQUENCE = 0, 1
AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWYX')

SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'scripts')

# results.json key -> output attribute
SCORE_ATTRS = {'affinity': '_boltzAffinity',
               'bindingProbability': '_boltzIntProbability',
               'confidence': '_boltzConfidence',
               'ligandIpTM': '_boltzIpTM',
               'plddt': '_boltzPLDDT',
               'iplddt': '_boltzIPLDDT',
               'receptorRMSD': '_boltzRecRMSD'}


class ProtBoltzCofolding(EMProtocol):
    """Cofold a receptor with each ligand of a set using Boltz-2, as a docking protocol.

    Every ligand is predicted in complex with the receptor chains, and Boltz-2's
    affinity module scores the binding. The receptor MSA is computed only once
    and shared by all ligand jobs, and all jobs run in a single Boltz call, so
    the model is loaded once.

    The receptor is given as a structure (AtomStruct) or as protein sequences.

    Each predicted complex is superposed onto the input receptor (CA atoms), so
    all poses share its frame. Every output pose carries its own predicted
    receptor as proteinFile, and the set keeps the input receptor. With sequence
    input there is no input structure: the first prediction is the reference
    frame and the set's receptor.

    Output columns:
      - _boltzAffinity: predicted log10(IC50) in uM (lower binds stronger).
        Computed on the top-ranked model and copied to every pose of the ligand.
      - _boltzIntProbability: probability that the ligand is a binder (0-1).
      - _boltzConfidence: Boltz confidence score used to rank the models.
      - _boltzIpTM: ligand interface pTM.
      - _boltzPLDDT / _boltzIPLDDT: mean pLDDT of the complex / of the interface.
      - _boltzRecRMSD: CA RMSD (A) of the predicted receptor to the input one
        (to the reference prediction with sequence input). Large values mean
        Boltz did not reproduce that conformation.

    Affinity is only supported for ligands up to 128 heavy atoms; it was trained
    on ligands up to 56.

    Reference: Passaro et al., bioRxiv 2025 (Boltz-2)
    """
    _label = 'boltz-2 cofolding'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addHidden(USE_GPU, params.BooleanParam, default=True,
                       label="Use GPU for execution: ",
                       help="Boltz is very slow on CPU.")
        form.addHidden(GPU_LIST, params.StringParam, default='0', label="Choose GPU IDs",
                       help="Comma-separated GPU devices. Ligand jobs are distributed among them.")

        form.addSection(label='Input')
        group = form.addGroup('Input specifications')
        group.addParam('receptorOrigin', params.EnumParam, default=FROM_STRUCTURE,
                       label='Receptor from: ', choices=['Structure', 'Sequence'],
                       display=params.EnumParam.DISPLAY_HLIST)
        group.addParam('inputAtomStruct', params.PointerParam, pointerClass='AtomStruct',
                       label='Receptor structure: ', condition=f'receptorOrigin == {FROM_STRUCTURE}',
                       help='Its protein chains are cofolded with every ligand.')
        group.addParam('inputSequence', params.TextParam, width=80,
                       label='Receptor sequence: ', condition=f'receptorOrigin == {FROM_SEQUENCE}',
                       help='Protein sequence in one-letter code. Separate the chains of a '
                            'complex with ":" (e.g. MKV...:MKV... for a homodimer); they are '
                            'named A, B, ... Spaces and line breaks are ignored.')
        group.addParam('chains', params.StringParam, default='',
                       label='Receptor chains: ', condition=f'receptorOrigin == {FROM_STRUCTURE}',
                       help='Chains to include: use the wizard or type comma-separated ids '
                            '(e.g. A,B). Empty uses every protein chain.')
        group.addParam('inputSmallMolecules', params.PointerParam, pointerClass='SetOfSmallMolecules',
                       label='Ligand set: ')

        group = form.addGroup('MSA')
        group.addParam('useMsaServer', params.BooleanParam, default=True,
                       label='Compute MSA on server: ',
                       help='Compute the receptor MSA once with the MMseqs2 server. If disabled, '
                            'Boltz runs in single-sequence mode (faster, usually less accurate).')
        group.addParam('msaServerUrl', params.StringParam, default='https://api.colabfold.com',
                       label='MSA server URL: ', condition='useMsaServer',
                       expertLevel=params.LEVEL_ADVANCED)
        group.addParam('maxMsaSeqs', params.IntParam, default=8192,
                       label='Max MSA sequences: ', expertLevel=params.LEVEL_ADVANCED,
                       help='Lower it (e.g. 512-1024) to reduce GPU memory use.')

        form.addSection(label='Parameters')
        form.addParam('diffusionSamples', params.IntParam, default=4,
                      label='Poses per ligand: ',
                      help='Number of diffusion samples (predicted complexes) per ligand.')
        form.addParam('infPot', params.BooleanParam, default=True,
                      label="Inference potentials: ",
                      help='Use steering potentials to improve the physical plausibility of the poses.')
        form.addParam('recyclingSteps', params.IntParam, default=3, expertLevel=params.LEVEL_ADVANCED,
                      label='Recycling steps: ')
        form.addParam('samplingSteps', params.IntParam, default=200, expertLevel=params.LEVEL_ADVANCED,
                      label='Sampling steps: ')
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
                      label='Diffusion samples for affinity: ')

        form.addParallelSection(threads=2, mpi=0)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.prepareStep)
        self._insertFunctionStep(self.predictStep, needsGPU=True)
        self._insertFunctionStep(self.extractStep)
        self._insertFunctionStep(self.createOutputStep)

    def prepareStep(self):
        """Ligands to SDF, then receptor MSA (once) + one Boltz YAML per ligand"""
        ligands = [{'id': ligId, 'file': os.path.abspath(convertToSdf(self, mol.getFileName()))}
                   for ligId, mol in self.getInputMolsDic().items()]
        ligandsJson = self._getExtraPath('ligands.json')
        with open(ligandsJson, 'w') as f:
            json.dump(ligands, f)

        chains = self.getChainsArg()
        msaServer = self.msaServerUrl.get() if self.useMsaServer.get() else 'none'
        receptor = self.getSequence() if self.isFromSequence() else self.getOriginalReceptorFile(getLink=False)
        args = f'prepare "{receptor}" "{os.path.abspath(ligandsJson)}" ' \
               f'"{self.getWorkDir()}" {chains} {msaServer}'
        Plugin.runScript(self, 'boltzCofold.py', args, env=BOLTZ_DIC, scriptDir=SCRIPTS_DIR)

    def predictStep(self):
        """All ligand jobs in one Boltz call. Finished jobs are skipped on continue."""
        gpus = self.getGpuIds()
        args = [f'"{os.path.join(self.getWorkDir(), "inputs")}"',
                f'--out_dir "{self.getWorkDir()}"',
                f'--cache "{os.path.abspath(os.path.join(Plugin.getVar(BOLTZ_DIC["home"]), "mol"))}"',
                f'--recycling_steps {self.recyclingSteps.get()}',
                f'--sampling_steps {self.samplingSteps.get()}',
                f'--diffusion_samples {self.diffusionSamples.get()}',
                f'--max_parallel_samples {self.maxParallelSamples.get()}',
                f'--step_scale {self.stepScale.get()}',
                f'--max_msa_seqs {self.maxMsaSeqs.get()}',
                f'--diffusion_samples_affinity {self.diffusionSamplesAff.get()}',
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

    def extractStep(self):
        """Superpose each prediction onto the receptor; write poses, receptors and scores"""
        predDir = os.path.join(self.getResultsDir(), 'predictions')
        args = f'extract "{self.getWorkDir()}" "{predDir}" "{self.getOutputDir()}"'
        Plugin.runScript(self, 'boltzCofold.py', args, env=BOLTZ_DIC, scriptDir=SCRIPTS_DIR)

    def createOutputStep(self):
        with open(os.path.join(self.getOutputDir(), 'results.json')) as f:
            results = json.load(f)
        if not results:
            raise Exception('Boltz produced no prediction. Check the run log '
                            '(e.g. GPU out of memory, unreadable ligands).')

        inputMols = self.getInputMolsDic()
        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
        for res in results:
            newMol = SmallMolecule()
            newMol.copy(inputMols[res['id']], copyId=False)
            newMol.setPoseFile(os.path.relpath(res['poseFile']))
            newMol.setProteinFile(os.path.relpath(res['receptorFile']))
            newMol.setPoseId(res['poseId'])
            newMol.setGridId(1)
            newMol.setMolClass('Boltz')
            newMol.setDockId(self.getObjId())
            # Every attribute on every item: a Set fixes its columns from the first one
            for key, attr in SCORE_ATTRS.items():
                setattr(newMol, attr, pwobj.Float(res[key]))
            outputSet.append(newMol)

        outputSet.updateMolClass()
        outputSet.setProteinFile(self.getOriginalReceptorFile())
        outputSet.setDocked(True)
        self._defineOutputs(outputSmallMolecules=outputSet)
        self._defineSourceRelation(self.inputSmallMolecules, outputSet)

    # --------------------------- UTILS functions -----------------------------------
    def getWorkDir(self):
        return os.path.abspath(self._getExtraPath('boltz'))

    def getResultsDir(self):
        return os.path.join(self.getWorkDir(), 'boltz_results')

    def getOutputDir(self):
        return os.path.abspath(self._getPath('outputLigands'))

    def getInputMolsDic(self):
        """Boltz job id -> input molecule. The id is also the YAML/record name, so it
        must be filesystem-safe and unique (molName alone may repeat across conformers)."""
        return {f"{re.sub(r'[^A-Za-z0-9_.-]', '_', mol.getMolName())}_{mol.getObjId()}": mol.clone()
                for mol in self.inputSmallMolecules.get()}

    def getGpuIds(self):
        if not getattr(self, USE_GPU).get():
            return []
        return [g.strip() for g in getattr(self, GPU_LIST).get().split(',') if g.strip()]

    def getChainsArg(self):
        """Chain ids for the script: 'all', or 'A,B' from typed ids or the wizard's JSON"""
        value = (self.chains.get() or '').strip()
        if not value:
            return 'all'
        try:
            chainIds = getChainIds(value)
        except ValueError:  # typed ids, not wizard JSON
            chainIds = value.split(',')
        return ','.join(c.strip() for c in chainIds)

    def isFromSequence(self):
        return self.receptorOrigin.get() == FROM_SEQUENCE

    def getSequence(self):
        return re.sub(r'\s', '', self.inputSequence.get() or '').upper()

    def getOriginalReceptorFile(self, getLink=True):
        """The input receptor, by default the hard link kept inside the protocol (used by viewers).
        With sequence input, the reference predicted receptor (exists after the extract step)."""
        if self.isFromSequence():
            return os.path.relpath(os.path.join(self.getOutputDir(), 'reference_receptor.pdb'))
        recFile = self.inputAtomStruct.get().getFileName()
        if not getLink:
            return os.path.abspath(recFile)
        recDir = self._getExtraPath('originalReceptor')
        os.makedirs(recDir, exist_ok=True)
        recLink = os.path.join(recDir, getBaseFileName(recFile))
        if not os.path.exists(recLink):
            os.link(recFile, recLink)
        return recLink

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if self.isFromSequence():
            chains = self.getSequence().split(':')
            summary.append(f'Receptor: {len(chains)} sequence(s), {sum(map(len, chains))} residues')
        elif self.inputAtomStruct.get():
            summary.append(f'Receptor: {os.path.basename(self.inputAtomStruct.get().getFileName())}')
        if self.inputSmallMolecules.get():
            summary.append(f'Ligands: {self.inputSmallMolecules.get().getSize()} molecule(s)')
        if self.hasAttribute('outputSmallMolecules'):
            summary.append(f'Output poses: {self.outputSmallMolecules.getSize()}')
        return summary

    def _validate(self):
        errors = []
        if self.isFromSequence():
            chains = self.getSequence().split(':')
            if not all(chains):
                errors.append('Empty receptor sequence (or empty chain between ":").')
            badChars = set(''.join(chains)) - AMINO_ACIDS
            if badChars:
                errors.append(f'Invalid characters in the receptor sequence: {"".join(sorted(badChars))}')
        elif not self.inputAtomStruct.get():
            errors.append('A receptor structure is required.')
        return errors

    def _methods(self):
        return ['Protein-ligand complexes were cofolded and their binding affinity predicted '
                'with Boltz-2 (Passaro et al. 2025).']
