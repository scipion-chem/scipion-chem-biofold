# **************************************************************************
# *
# * Authors:   Blanca Pueche (blanca.pueche@cnb.csis.es)
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
import string

import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pyworkflow.object import String

from pwchem import Plugin
from biofold.constants import PROTENIX_DIC

from pwem.objects import  AtomStruct, SetOfAtomStructs
from pwchem.utils.utilsFasta import parseFasta


class ProtProtenix(EMProtocol):
    """
    Protocol to use Protenix model.

    AI Generated:

       #todo
    """
    _label = 'protenix modelling'
    models = ['protenix_base_default_v1.0.0', 'protenix_base_20250630_v1.0.0', 'protenix_base_default_v0.5.0']

    # -------------------------- DEFINE param functions ----------------------
    def _addInputForm(self, form):
        form.addParam('inputSequence', params.PointerParam,
                      pointerClass='Sequence', allowsNull=True,
                      label="Input sequence: ", condition='inputOrigin==0',
                      help='Select the sequence object to add to the set')

        form.addParam('inputAtomStruct', params.PointerParam,
                      pointerClass='AtomStruct', allowsNull=True,
                      label="Input structure: ", condition='inputOrigin==1',
                      help='Select the AtomStruct object whose sequence to add to the set')

        form.addParam('entityType', params.EnumParam, default=0, condition='inputOrigin in [0,1]',
                      label='Input entity type: ', choices=['Protein', 'DNA', 'RNA'],
                      help='Input entity type to add to the set')

        form.addParam('inpChain', params.StringParam,
                      label='Input chain: ', condition='inputOrigin == 1',
                      help='Specify the protein chain to use as sequence.')

        form.addParam('inpPositions', params.StringParam,
                      label='Input positions: ', condition='inputOrigin in [0,1]',
                      help='Specify the positions of the sequence to add in the output.')

        form.addParam('addInput', params.LabelParam,
                      label='Add input: ', condition='inputOrigin in [0,1]',
                      help='Add sequence to the output set')

    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='Input')
        form.addParam('inputOrigin', params.EnumParam, default=0,
                      label='Input origin: ', choices=['Sequence', 'AtomStruct', 'fasta file'],
                      help='Input origin to add to the set')
        self._addInputForm(form)

        form.addParam('inputList', params.TextParam, width=100, condition='inputOrigin in [0,1]',
                      default='', label='List of inputs: ',
                      help='The list of input to use for the final output set.')

        form.addParam('file', params.FileParam, condition='inputOrigin == 2',
                      label='Sequence file: ',
                      help='Select the fasta file.')


        group = form.addGroup('Parameters')
        group.addParam('msa', params.BooleanParam, default=True,
                      label="Run with MSAs: ",
                      help='Choose whether to run with MSAs for improved performance.')
        group.addParam('sample', params.IntParam, default=5,
                      label='Diffusion samples for affinity: ', help="Number of diffusion samples for affinity.")
        group.addParam('steps', params.IntParam, default=200,
                        label='Sampling steps: ', help="Number of sampling steps for prediction.")
        group.addParam('model', params.EnumParam, default=0, choices=self.models,
                        label='Model: ', help="Model options:\n\n- protenix_base_default_v1.0.0: Default model trained with a data cutoff aligned with AlphaFold3 (2021-09-30). "
                                              "Recommended for fair benchmarking and comparison with other state-of-the-art methods.\n\n- "
                                              "protenix_base_20250630_v1.0.0: Updated model trained with a more recent data cutoff (2025-06-30). "
                                              "Recommended for practical applications and best performance.\n\n- protenix_base_default_v0.5.0: Previous model version. "
                                              "Maintained for backward compatibility with older workflows.")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        if self.inputOrigin.get() == 2:
            self._insertFunctionStep(self.createProtenixJsonFromFastaStep)
        else:
            self._insertFunctionStep(self.createProtenixInputFileStep)

        self._insertFunctionStep(self.runProtenixStep)
        #self._insertFunctionStep(self.extractScoreStep)
        #self._insertFunctionStep(self.createOutputStep)

    def createProtenixJsonFromFastaStep(self):
        fastaPath = os.path.abspath(self.file.get())
        seqDic = parseFasta(fastaPath)

        chainIdIter = iter(string.ascii_uppercase)
        sequences = []

        for _, sequence in seqDic.items():
            chainId = next(chainIdIter)

            entityType = self.guessEntityType(sequence)

            sequences.append(
                self._buildProtenixEntity(entityType, sequence, chainId)
            )

        jsonPath = os.path.abspath(self._getPath("protenix_input.json"))

        with open(jsonPath, "w") as f:
            json.dump([
                {
                    "name": "intellifold_job",
                    "sequences": sequences
                }
            ], f, indent=2)

    def createProtenixInputFileStep(self):
        chainIdIter = iter(string.ascii_uppercase)
        sequences = []

        for inputLine in self.inputList.get().split('\n'):
            if not inputLine.strip():
                continue

            inpJson = json.loads(inputLine.split(')')[1].strip())
            seqDic = parseFasta(os.path.abspath(inpJson['seqFile']))
            _, sequence = next(iter(seqDic.items()))

            entityType = inpJson.get('entity', 'protein')
            chainId = next(chainIdIter)

            sequences.append(
                self._buildProtenixEntity(entityType, sequence, chainId)
            )

        jsonPath = os.path.abspath(self._getPath("protenix_input.json"))

        with open(jsonPath, "w") as f:
            json.dump([
                {
                    "name": "intellifold_job",
                    "sequences": sequences
                }
            ], f, indent=2)

    def runProtenixStep(self):
        jsonPath = os.path.abspath(self._getPath("protenix_input.json"))
        args = ['-i', str(jsonPath)]

        args.append(f'-o {os.path.join(os.path.abspath(self._getPath()), "protenix_results")}')

        if self.msa.get():
            args.append("--use_msa TRUE")
        else:
            args.append("--use_msa FALSE")

        args.append(f" --sample {self.sample.get()}")
        args.append(f" --step {self.steps.get()}")
        args.append(f" --model_name {self.getEnumText('model')}")
        args.append(f" --nhmmer_n_cpu {self.numberOfThreads.get()}")

        program_prefix = (
            "export CUDA_HOME=$CONDA_PREFIX && "
            "export PROTENIX_DEVICE=cpu && "
            "protenix pred"
        )

        Plugin.runCondaCommand(
            self,
            args=" ".join(args),
            condaDic=PROTENIX_DIC,
            program=program_prefix,
            cwd=os.path.abspath(Plugin.getVar(PROTENIX_DIC['home']))
        )


    def extractScoreStep(self):
        """Extract per-residue score and compute mean score per model"""
        resultsPath = os.path.join(os.path.abspath(self._getPath()), "chai_results")
        self.meanScore = {}

        extraFiles = self.getExtraFiles()

        for cifPath in extraFiles:
            cifName = os.path.basename(cifPath)
            modelName = os.path.splitext(cifName)[0]

            headers = []
            scoreValues = []
            seenResidues = set()

            with open(cifPath) as f:
                for line in f:
                    if line.startswith('_atom_site.'):
                        headers.append(line.strip())
                    elif line.startswith('ATOM'):
                        break

            colIndex = {h.split('.')[-1]: i for i, h in enumerate(headers)}

            if 'B_iso_or_equiv' not in colIndex:
                raise Exception(f"No score field in {cifName}")

            with open(cifPath) as f:
                for line in f:
                    if not line.startswith('ATOM'):
                        continue
                    cols = re.sub(r'\s+', ' ', line.strip()).split()
                    resnum = int(cols[colIndex['auth_seq_id']])
                    score = float(cols[colIndex['B_iso_or_equiv']])
                    if resnum not in seenResidues:
                        seenResidues.add(resnum)
                        scoreValues.append(score)

            self.meanScore[modelName] = sum(scoreValues) / len(scoreValues)

        self.bestModel = max(self.meanScore, key=self.meanScore.get)

    def createOutputStep(self):
        extraFiles = self.getExtraFiles()
        outputSet = SetOfAtomStructs.create(self._getPath())

        bestSrc = None

        for cifPath in extraFiles:
            atomStruct = AtomStruct(filename=cifPath)
            atomStruct.origin = String()
            atomStruct.setAttributeValue('origin', 'Chai')
            outputSet.append(atomStruct)

            modelName = os.path.splitext(os.path.basename(cifPath))[0]
            if modelName == self.bestModel:
                bestSrc = cifPath

        if bestSrc is None:
            raise Exception(f"Best model {self.bestModel} not found among output files.")

        bestStruct = AtomStruct(filename=bestSrc)
        bestStruct.origin = String()
        bestStruct.setAttributeValue('origin', 'Chai')

        self._defineOutputs(
            outputBestAtomStruct=bestStruct,
            outputSetOfAtomStructs=outputSet
        )

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []

        resultsPath = os.path.join(os.path.abspath(self._getPath()), "chai_results")

        try:
            extraFiles = self.getExtraFiles()
        except Exception as e:
            summary.append(f"No CIF files found in {resultsPath}: {e}")
            return summary

        if not extraFiles:
            summary.append(f"No CIF files found in {resultsPath}.")
            return summary

        scores = {}
        logFile = os.path.join(self._getPath('logs'), "run.stdout")

        if os.path.exists(logFile):
            with open(logFile) as f:
                for line in f:
                    m = re.search(r"Score=([\d.]+), writing output to (.+\.cif)", line)
                    if m:
                        score = float(m.group(1))
                        cifPath = m.group(2)

                        relPath = os.path.relpath(cifPath, resultsPath)
                        modelName = os.path.splitext(relPath)[0]

                        scores[modelName] = score

        if not scores:
            summary.append("No scores detected in log file.")
            return summary

        summary.append("Scores per model:")

        for modelName in sorted(scores.keys()):
            summary.append(f"  {modelName}: {scores[modelName]:.4f}")

        bestModel = max(scores, key=scores.get)
        summary.append(f"\nBest structure (highest score): {bestModel}.cif")

        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def getExtraFiles(self):
        extraFiles = []
        resultsPath = os.path.join(os.path.abspath(self._getPath()), "chai_results")

        for root, dirs, files in os.walk(resultsPath):
            for name in files:
                if name.lower().endswith(".cif"):
                    extraFiles.append(os.path.join(root, name))

        if not extraFiles:
            raise Exception("No CIF files found in chai_results.")

        return extraFiles

    def guessEntityType(self, sequence):
        seq = sequence.upper()
        dnaLetters = set("ACGT")
        rnaLetters = set("ACGU")
        proteinLetters = set("ACDEFGHIKLMNPQRSTVWY")

        seqSet = set(seq)

        # check for RNA (U present, T absent)
        if "U" in seqSet and "T" not in seqSet:
            return "rna"
        # check for DNA (T present, U absent)
        elif "T" in seqSet and "U" not in seqSet and seqSet <= dnaLetters:
            return "dna"
        # check for protein: contains amino acid letters not in DNA/RNA
        elif seqSet <= proteinLetters:
            return "protein"
        else:
            return "protein"

    def ensureFastaHasNames(self):
        fastaPath = os.path.abspath(self.file.get())
        outputPath = os.path.join(self._getPath('input.fasta'))

        with open(fastaPath) as f:
            lines = f.readlines()

        fixedLines = []
        counter = 1
        sequenceBuffer = ""
        currentHeader = None


        for line in lines:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if currentHeader is not None:
                    entityType = self.guessEntityType(sequenceBuffer)
                    self.writeSequence(entityType, currentHeader, sequenceBuffer, counter, fixedLines)
                    counter += 1

                header = line[1:]
                if "|name=" in header:
                    header = header.split("|name=")[-1]
                currentHeader = header
                sequenceBuffer = ""
            else:
                sequenceBuffer += line.strip()

        if currentHeader is not None:
            entityType = self.guessEntityType(sequenceBuffer)
            self.writeSequence(entityType, currentHeader, sequenceBuffer, counter, fixedLines)
            counter += 1

        with open(outputPath, 'w') as f:
            f.writelines(fixedLines)

        self.NEWFILE = True

    def writeSequence(self, entityType, header, sequence, counter, fixedLines):
        if not header:
            header = f"seq{counter}"
        fixedLines.append(f">{entityType}|{header}\n")
        for i in range(0, len(sequence), 80):
            fixedLines.append(sequence[i:i + 80] + "\n")

    def _buildProtenixEntity(self, entityType, sequence, chainId):
        entityType = entityType.lower()
        if not sequence:
            sequence = "X"  # placeholder sequence to avoid empty
        if entityType == "protein":
            return {"proteinChain": {"id": chainId, "sequence": sequence}}
        elif entityType == "dna":
            return {"dnaSequence": {"id": chainId, "sequence": sequence}}
        elif entityType == "rna":
            return {"rnaSequence": {"id": chainId, "sequence": sequence}}
        else:
            return {"proteinChain": {"id": chainId, "sequence": sequence}}
