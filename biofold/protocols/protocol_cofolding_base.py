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
"""Shared base for the cofolding protocols (Boltz-2, Chai-1, Protenix, IntelliFold).

Every backend predicts the complex of a receptor with each ligand of a set and returns
a docked SetOfSmallMolecules, so the form, the ligand and receptor handling, the
superposition and the output building are all the same. A subclass supplies its own
prediction parameters, its command line and where its predictions land.
"""
import json
import os
import re

import pyworkflow.object as pwobj
import pyworkflow.protocol.params as params
from pyworkflow.protocol.params import USE_GPU, GPU_LIST
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.constants import RDKIT_DIC, SCORCH2_DIC
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import convertToSdf, getBaseFileName, getChainIds

SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'scripts')

FROM_STRUCTURE, FROM_SEQUENCE, FROM_SET = 0, 1, 2
FROM_MOLS, FROM_LIBRARY = 0, 1
AMINO_ACIDS = set('ACDEFGHIKLMNPQRSTVWYX')

COFOLD_ATTRS = {'ligPLDDT': '_cofoldLigPLDDT',
                'ligIpTM': '_cofoldLigIpTM',
                'confidence': '_cofoldConfidence',
                'plddt': '_cofoldPLDDT',
                'receptorRMSD': '_cofoldRecRMSD'}

COFOLD_ATTRS_HELP = """
Output columns (all on a 0-1 scale except the RMSD):

- *_cofoldLigPLDDT*: mean pLDDT over the ligand atoms only. The best single measure of
  how confident the model is about the pose itself.
- *_cofoldLigIpTM*: ipTM of the protein-ligand interface. With several protein chains it
  is the best interface the ligand makes with any one of them.
- *_cofoldConfidence*: the backend's own ranking score, which it uses to order its models.
- *_cofoldPLDDT*: pLDDT over the whole complex. Dominated by the protein, so it says
  whether the fold worked at all, not whether the pose is good. Do not rank with it.
- *_cofoldRecRMSD*: CA RMSD (A) of the predicted receptor against the input one. Large
  values mean the backend did not reproduce the input conformation.
"""


def safeName(name):
    """Filesystem-safe version of a molecule name"""
    return re.sub(r'[^A-Za-z0-9_.-]', '_', name)


class ProtCofoldingBase(EMProtocol):
    """Cofold a receptor with each ligand of a set. Not used directly."""
    _label = None

    #: results.json key -> output attribute; subclasses extend COFOLD_ATTRS
    SCORE_ATTRS = COFOLD_ATTRS
    #: name passed to cofoldUtils.py extract, and the molClass of the output molecules
    BACKEND = None
    MOL_CLASS = None
    #: conda env of the tool whose MSA client computes the receptor MSA, or None when the
    #: tool computes it itself while predicting
    MSA_ENV = None

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addHidden(USE_GPU, params.BooleanParam, default=True,
                       label="Use GPU for execution: ",
                       help="Cofolding is very slow on CPU.")
        form.addHidden(GPU_LIST, params.StringParam, default='0', label="Choose GPU IDs",
                       help="Comma-separated GPU devices.")

        form.addSection(label='Input')
        group = form.addGroup('Input specifications')
        group.addParam('receptorOrigin', params.EnumParam, default=FROM_STRUCTURE,
                       label='Receptor from: ', choices=['Structure', 'Sequence', 'SetOfSequences'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       help='Where the receptor comes from:\n\n'
                            '- *Structure*: an AtomStruct object. Its protein chains (or the ones '
                            'selected) are cofolded, and every prediction is superposed onto it.\n'
                            '- *Sequence*: type or paste the protein sequence in one-letter code; '
                            'separate the chains of a complex with ":".\n'
                            '- *SetOfSequences*: pick one or several sequences from a '
                            'SetOfSequences object; each becomes one chain.\n\n'
                            'With Sequence or SetOfSequences there is no input structure: the '
                            'first prediction is used as the reference frame and receptor.')
        group.addParam('inputAtomStruct', params.PointerParam, pointerClass='AtomStruct',
                       label='Receptor structure: ', condition=f'receptorOrigin == {FROM_STRUCTURE}',
                       help='Its protein chains are cofolded with every ligand.')
        group.addParam('inputSequence', params.TextParam, width=80,
                       label='Receptor sequence: ', condition=f'receptorOrigin == {FROM_SEQUENCE}',
                       help='Protein sequence in one-letter code. Separate the chains of a '
                            'complex with ":" (e.g. MKV...:MKV... for a homodimer); they are '
                            'named A, B, ... Spaces and line breaks are ignored.')
        group.addParam('inputSetOfSequences', params.PointerParam, pointerClass='SetOfSequences',
                       allowsNull=True, label='Receptor sequences: ',
                       condition=f'receptorOrigin == {FROM_SET}')
        group.addParam('chooseSeqs', params.StringParam, default='All',
                       label='Sequences to use: ', condition=f'receptorOrigin == {FROM_SET}',
                       help='Sequences of the set that form the receptor, one chain each (named '
                            'A, B, ... in set order). Use the wizard, or "All".')
        group.addParam('chains', params.StringParam, default='',
                       label='Receptor chains: ', condition=f'receptorOrigin == {FROM_STRUCTURE}',
                       help='Chains to include: use the wizard or type comma-separated ids '
                            '(e.g. A,B). Empty uses every protein chain. Fewer chains make the '
                            'prediction much cheaper.')
        group.addParam('ligandOrigin', params.EnumParam, default=FROM_MOLS,
                       label='Ligands from: ', choices=['SetOfSmallMolecules', 'SmallMoleculesLibrary'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       help='- *SetOfSmallMolecules*: molecule files (converted to SMILES).\n'
                            '- *SmallMoleculesLibrary*: a SMILES table; the SMILES are used directly, '
                            'with no file conversion.')
        group.addParam('inputSmallMolecules', params.PointerParam, pointerClass='SetOfSmallMolecules',
                       label='Ligand set: ', condition=f'ligandOrigin == {FROM_MOLS}', allowsNull=True)
        group.addParam('inputLibrary', params.PointerParam, pointerClass='SmallMoleculesLibrary',
                       label='Ligand library: ', condition=f'ligandOrigin == {FROM_LIBRARY}', allowsNull=True)
        group.addParam('smilesHeader', params.StringParam, default='SMI',
                       label='SMILES column: ', condition=f'ligandOrigin == {FROM_LIBRARY}',
                       help='Library column holding the SMILES. Use the wizard to pick it. Ligands are '
                            'named from the "NAME" column when the library has one.')

        group = form.addGroup('MSA')
        group.addParam('useMsaServer', params.BooleanParam, default=True,
                       label='Compute MSA on server: ',
                       help='Compute the receptor MSA once with the MMseqs2 server and reuse it for '
                            'every ligand. If disabled, the backend runs in single-sequence mode '
                            '(faster, usually less accurate).')
        group.addParam('msaServerUrl', params.StringParam, default='https://api.colabfold.com',
                       label='MSA server URL: ', condition='useMsaServer',
                       expertLevel=params.LEVEL_ADVANCED)
        self._defineExtraMsaParams(group)

        group = form.addGroup('Stereochemistry')
        group.addParam('checkStereo', params.BooleanParam, default=True,
                       label='Discard wrong stereochemistry: ',
                       help='Cofolding models place the atoms by diffusion, so a stereocentre can '
                            'come out mirrored even though the input defines it.\n\n'
                            'With this on, every pose is checked against its input SMILES with '
                            "PoseBusters' identity test, and any pose whose tetrahedral "
                            'stereochemistry disagrees is left out of the output.')

        form.addSection(label='Parameters')
        self._definePredictionParams(form)
        form.addParallelSection(threads=2, mpi=0)

    def _defineExtraMsaParams(self, group):
        """Backend-specific MSA parameters. Override when needed."""
        pass

    def _definePredictionParams(self, form):
        """The backend's own prediction parameters. Subclasses must implement."""
        raise NotImplementedError

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.prepareStep)
        self._insertFunctionStep(self.predictStep, needsGPU=True)
        self._insertFunctionStep(self.extractStep)
        if self.checkStereo.get():
            self._insertFunctionStep(self.stereoStep)
        self._insertFunctionStep(self.createOutputStep)

    def prepareStep(self):
        """Ligands to SMILES/SDF, then the backend's receptor MSA (once) and input files"""
        if self.isFromLibrary():
            # One .smi per ligand, so each output molecule has its own input file
            os.makedirs(self._getExtraPath('libraryLigands'), exist_ok=True)
            ligands = []
            for ligId, (name, smi) in self.getLibraryLigands().items():
                with open(self.getLibraryMolFile(ligId), 'w') as f:
                    f.write(f'{smi}\t{name}\n')
                ligands.append({'id': ligId, 'smiles': smi})
        else:
            ligands = [{'id': ligId, 'file': os.path.abspath(convertToSdf(self, mol.getFileName()))}
                       for ligId, mol in self.getInputMolsDic().items()]
        ligandsJson = self._getExtraPath('ligands.json')
        with open(ligandsJson, 'w') as f:
            json.dump(ligands, f)

        # Receptor chains and ligand SMILES are resolved once, in the RDKit env
        args = f'receptor "{self.getReceptorArg()}" "{os.path.abspath(ligandsJson)}" ' \
               f'"{self.getWorkDir()}" {self.getChainsArg()}'
        Plugin.runScript(self, 'cofoldUtils.py', args, env=RDKIT_DIC, scriptDir=SCRIPTS_DIR)

        self.runPrepare()

    def runPrepare(self):
        """Compute the receptor MSA once, then write one job file per ligand.

        Only the MSA needs the prediction tool's env, since its client ships with the tool.
        The job files are plain text, so they are written here, in the protocol.
        """
        rec = self.getReceptorInfo()
        msa = {}
        if self.MSA_ENV is not None and self.useMsaServer.get():
            args = f'{self.BACKEND} "{self.getWorkDir()}" {self.getMsaArgs()}'
            Plugin.runScript(self, 'cofoldMsa.py', args, env=self.MSA_ENV, scriptDir=SCRIPTS_DIR)
            with open(os.path.join(self.getWorkDir(), 'msa.json')) as f:
                msa = json.load(f)

        os.makedirs(self.getInputsDir(), exist_ok=True)
        for ligId, smi in rec['ligands'].items():
            self.writeJob(ligId, smi, rec, msa)

    def getMsaArgs(self):
        """What follows the backend name and work dir on cofoldMsa.py's command line"""
        return self.msaServerUrl.get()

    def writeJob(self, ligId, smiles, rec, msa):
        """Write the backend's input file for one ligand into the inputs directory.

        `rec` is <workDir>/receptor.json and `msa` maps each receptor sequence to the
        job-file keys that point at its MSA (empty when no MSA was computed).
        Subclasses must implement.
        """
        raise NotImplementedError

    def getReceptorInfo(self):
        """<workDir>/receptor.json: receptor chains, ligand SMILES and the ligand chain id"""
        with open(os.path.join(self.getWorkDir(), 'receptor.json')) as f:
            return json.load(f)

    @staticmethod
    def getChainsBySequence(rec):
        """{receptor sequence: [chain ids]}. Identical chains share one entity, and so
        one MSA."""
        entities = {}
        for ch in rec['chains']:
            entities.setdefault(ch['seq'], []).append(ch['id'])
        return entities

    def predictStep(self):
        """Run the backend. Subclasses must implement."""
        raise NotImplementedError

    def extractStep(self):
        """Superpose each prediction onto the receptor; write poses, receptors and scores.

        Shared by every backend and needs no code from the prediction tool, so it runs in
        pwchem's RDKit env rather than the tool's own.
        """
        args = (f'extract {self.BACKEND} "{self.getWorkDir()}" '
                f'"{self.getPredictionsDir()}" "{self.getOutputDir()}"')
        Plugin.runScript(self, 'cofoldUtils.py', args, env=RDKIT_DIC, scriptDir=SCRIPTS_DIR)

    def stereoStep(self):
        """Runs in pwchem's scorch2 env, the one its PoseBusters package is installed into."""
        args = f'"{self.getWorkDir()}" "{self.getOutputDir()}"'
        Plugin.runScript(self, 'cofoldStereo.py', args, env=SCORCH2_DIC, scriptDir=SCRIPTS_DIR)

    def getStereoVerdicts(self):
        """{pose file: stereoOK}. A pose that could not be checked is kept."""
        stereoFile = os.path.join(self.getOutputDir(), 'stereo.json')
        if not self.checkStereo.get() or not os.path.exists(stereoFile):
            return {}
        with open(stereoFile) as f:
            return {pose: v['stereoOK'] for pose, v in json.load(f).items()}

    def getPredictionsDir(self):
        """Directory holding one subdirectory of predictions per ligand id."""
        raise NotImplementedError

    def createOutputStep(self):
        with open(os.path.join(self.getOutputDir(), 'results.json')) as f:
            results = json.load(f)
        if not results:
            raise Exception(f'{self.MOL_CLASS} produced no prediction. Check the run log '
                            f'(e.g. GPU out of memory, unreadable ligands).')

        verdicts = self.getStereoVerdicts()
        kept = [r for r in results if verdicts.get(r['poseFile'], True) is not False]
        if not kept:
            raise Exception(f'Every pose has the wrong stereochemistry, so {self.MOL_CLASS} '
                            f'produced no output.')
        dropped = len(results) - len(kept)
        if dropped:
            lost = sorted({r['id'] for r in results} - {r['id'] for r in kept})
            self.info(f'Stereochemistry: {dropped} of {len(results)} poses discarded')
            if lost:
                self.info(f'No pose left for: {", ".join(lost)}')
        results = kept

        scoreKeys = [k for k in self.SCORE_ATTRS if any(r.get(k) is not None for r in results)]

        inputMols = self.getInputMolsDic()
        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
        for res in results:
            newMol = SmallMolecule()
            newMol.copy(inputMols[res['id']], copyId=False)
            newMol.setPoseFile(os.path.relpath(res['poseFile']))
            newMol.setProteinFile(os.path.relpath(res['receptorFile']))
            newMol.setPoseId(res['poseId'])
            newMol.setGridId(1)
            newMol.setMolClass(self.MOL_CLASS)
            newMol.setDockId(self.getObjId())
            # Every attribute on every item: a Set fixes its columns from the first one
            for key in scoreKeys:
                setattr(newMol, self.SCORE_ATTRS[key], pwobj.Float(res.get(key) or 0.0))
            outputSet.append(newMol)

        outputSet.updateMolClass()
        outputSet.setProteinFile(self.getOriginalReceptorFile())
        outputSet.setDocked(True)
        self._defineOutputs(outputSmallMolecules=outputSet)
        self._defineSourceRelation(self.inputLibrary if self.isFromLibrary() else self.inputSmallMolecules,
                                   outputSet)

    # --------------------------- UTILS functions -----------------------------------
    def getWorkDir(self):
        return os.path.abspath(self._getExtraPath(self.BACKEND))

    def getInputsDir(self):
        return os.path.join(self.getWorkDir(), 'inputs')

    def getOutputDir(self):
        return os.path.abspath(self._getPath('outputLigands'))

    def getInputMolsDic(self):
        """Job id -> input molecule. The id also names the backend's input file and its
        prediction directory, so it must be filesystem-safe and unique (molName alone may
        repeat across conformers)."""
        if self.isFromLibrary():
            return {ligId: SmallMolecule(smallMolFilename=os.path.relpath(self.getLibraryMolFile(ligId)), molName=name)
                    for ligId, (name, _) in self.getLibraryLigands().items()}
        return {f"{safeName(mol.getMolName())}_{mol.getObjId()}": mol.clone()
                for mol in self.inputSmallMolecules.get()}

    def isFromLibrary(self):
        return self.ligandOrigin.get() == FROM_LIBRARY

    def getLibrary(self):
        """Used by the SMILES column wizard"""
        return self.inputLibrary.get()

    def getLibraryLigands(self):
        """Job id -> (name, SMILES) for every library row with a SMILES.
        The library file has no header line; its columns follow getHeaders()."""
        library = self.getLibrary()
        headers = library.getHeaders()
        smiIdx = headers.index(self.smilesHeader.get())
        nameIdx = headers.index('NAME') if 'NAME' in headers else None
        ligands = {}
        with open(library.getFileName()) as f:
            for i, line in enumerate(f, 1):
                row = line.rstrip('\n').split(library.getSeparator())
                smi = row[smiIdx].strip() if len(row) > smiIdx else ''
                if not smi:
                    continue
                name = row[nameIdx].strip() if nameIdx is not None and len(row) > nameIdx else ''
                name = name or f'mol{i}'
                ligands[f'{safeName(name)}_{i}'] = (name, smi)
        return ligands

    def getLibraryMolFile(self, ligId):
        return os.path.abspath(self._getExtraPath('libraryLigands', f'{ligId}.smi'))

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
        """Receptor given as sequences (typed or from a set), so no input structure"""
        return self.receptorOrigin.get() != FROM_STRUCTURE

    def getSequence(self):
        """Receptor chains as SEQ1:SEQ2:..."""
        if self.receptorOrigin.get() == FROM_SET:
            chosen = {n.strip() for n in (self.chooseSeqs.get() or '').split(',')}
            seqs = [seq.getSequence() for seq in self.inputSetOfSequences.get()
                    if 'All' in chosen or seq.getSeqName() in chosen]
            return ':'.join(seqs).upper()
        return re.sub(r'\s', '', self.inputSequence.get() or '').upper()

    def getReceptorArg(self):
        """What the prepare script takes: a structure file, or SEQ1:SEQ2:..."""
        return self.getSequence() if self.isFromSequence() else self.getOriginalReceptorFile(getLink=False)

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
        if self.isFromLibrary() and self.inputLibrary.get():
            summary.append(f'Ligands: library of {self.inputLibrary.get().getLength()} molecule(s)')
        elif not self.isFromLibrary() and self.inputSmallMolecules.get():
            summary.append(f'Ligands: {self.inputSmallMolecules.get().getSize()} molecule(s)')
        if self.hasAttribute('outputSmallMolecules'):
            summary.append(f'Output poses: {self.outputSmallMolecules.getSize()}')
        return summary

    def _validate(self):
        errors = []
        if self.receptorOrigin.get() == FROM_SET and not self.inputSetOfSequences.get():
            errors.append('A receptor SetOfSequences is required.')
        elif self.isFromSequence():
            chains = self.getSequence().split(':')
            if not all(chains):
                errors.append('Empty receptor sequence: no sequence selected, or an empty chain between ":".')
            badChars = set(''.join(chains)) - AMINO_ACIDS
            if badChars:
                errors.append(f'Invalid characters in the receptor sequence: {"".join(sorted(badChars))}')
        elif not self.inputAtomStruct.get():
            errors.append('A receptor structure is required.')

        if self.isFromLibrary():
            if not self.inputLibrary.get():
                errors.append('A ligand library is required.')
            elif self.smilesHeader.get() not in self.inputLibrary.get().getHeaders():
                errors.append(f'SMILES column "{self.smilesHeader.get()}" is not in the library headers: '
                              f'{", ".join(self.inputLibrary.get().getHeaders())}')
        elif not self.inputSmallMolecules.get():
            errors.append('A ligand set is required.')
        return errors
