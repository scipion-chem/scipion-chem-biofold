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
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol
from pyworkflow.object import String

from pwchem import Plugin
from biofold.constants import CHAI_DIC

from pwem.objects import  AtomStruct, SetOfAtomStructs


class ProtChai(EMProtocol):
    """
    Protocol to use Chai-1 model.

    AI Generated:

        ProtChai - User Manual

        Overview
        --------
        The ProtChai protocol performs biomolecular structure prediction using the
        Chai-1 modelling framework. It supports proteins, DNA, and RNA inputs and
        generates 3D structural models in CIF format.

        This protocol integrates sequence preprocessing, FASTA preparation,
        Chai-1 execution, confidence score extraction, and structured output
        generation within the Scipion framework.

        Inputs
        ------
        Input origin can be selected from:

        1. **Sequence**:
           - A Sequence object defined within Scipion.
           - Supports protein, DNA, or RNA sequences.

        2. **AtomStruct**:
           - Extracts sequence from a provided structure.
           - Allows specifying chain ID and subsequence positions.

        3. **FASTA file**:
           - Direct input of a FASTA file containing one or more sequences.
           - Entity type (protein/DNA/RNA) is automatically inferred if not
             explicitly defined.
           - FASTA headers are reformatted to comply with Chai-1 requirements.

        Additional input options:
        - **entityType**: Define biological type (Protein, DNA, RNA).
        - **inpChain**: Specify chain when using AtomStruct input.
        - **inpPositions**: Select specific residue ranges.
        - **inputList**: Merge multiple inputs into a single modelling job.

        Prediction Parameters
        ---------------------
        - **msa**: Enable usage of the MSA server to improve prediction quality.
        - **trunkRecycles**: Number of trunk recycling refinement steps.
        - **timeSteps**: Number of diffusion sampling timesteps.
        - **trunkSamples**: Number of trunk samples per prediction.
        - **diffNsamples**: Number of diffusion samples for affinity estimation.

        Workflow
        --------
        1. **Input Processing**:
           - Parses sequences from Sequence objects, AtomStruct objects,
             or FASTA files.
           - Automatically formats sequences into a Chai-compatible FASTA file.
           - Infers entity type when required.

        2. **Model Execution**:
           - Executes `chai-lab fold` inside the configured Conda environment.
           - Applies selected recycling, sampling, and diffusion parameters.
           - Optionally connects to the MSA server.
           - Stores generated CIF files in the protocol working directory.

        3. **Score Extraction**:
           - Reads per-residue confidence scores from the `B_iso_or_equiv`
             field in the output CIF files.
           - Computes the mean score per predicted model.
           - Automatically selects the highest-scoring structure.

        4. **Output Creation**:
           - Wraps all predicted CIF models into a SetOfAtomStructs object.
           - Identifies and returns the best predicted structure separately.

        Outputs
        -------
        - **outputBestAtomStruct**: Best predicted 3D structure in CIF format
          (highest mean confidence score).
        - **outputSetOfAtomStructs**: Set containing all predicted structures.
        - Prediction results stored in the protocol output directory.

        Practical Recommendations
        -------------------------
        - Enable MSA usage for improved structural accuracy.
        - Increase recycling steps for enhanced refinement at the cost of runtime.
        - Adjust diffusion and sampling parameters depending on sequence length.
        - Use multiple trunk samples for more robust predictions in difficult cases.

        Interpretation
        --------------
        The resulting AtomStruct corresponds to the top-ranked Chai-1 prediction
        according to the mean per-residue confidence score.

        The output structures can be used for:

        - Structural visualization
        - Docking studies
        - Interaction analysis
        - Molecular dynamics simulations
        - Comparative structural modelling

        Warnings
        --------
        - Large sequences increase computational time and memory usage.
        - MSA server availability may affect runtime.
        - Ensure sequences are biologically valid and correctly formatted.

        Final Perspective
        -----------------
        ProtChai provides an integrated interface for Chai-1-based biomolecular
        structure prediction within Scipion. It streamlines sequence preparation,
        model execution, score evaluation, and structured result retrieval,
        enabling reproducible and automated structural modelling workflows.
    """
    _label = 'chai-1 modelling'
    NEWFILE = False

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


        form = form.addGroup('Parameters')
        form.addParam('msa', params.BooleanParam, default=True,
                      label="Run with MSAs: ",
                      help='Choose whether to run with MSAs for improved performance.')
        form.addParam('trunkRecycles', params.IntParam, default=3, expertLevel=params.LEVEL_ADVANCED,
                        label='Recycling steps: ', help="Number of recycling steps for prediction.")
        form.addParam('timeSteps', params.IntParam, default=200,
                        label='Sampling steps: ', help="Number of sampling steps for prediction.")
        form.addParam('trunkSamples', params.IntParam, default=1, expertLevel=params.LEVEL_ADVANCED,
                        label='Trunk samples: ', help="Number of trunk samples for prediction.")
        form.addParam('diffNsamples', params.IntParam, default=5, expertLevel=params.LEVEL_ADVANCED,
                       label='Difussion samples for affinity: ', help="Number of diffusion samples for affinity.")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        if (self.inputOrigin.get() != 2):
            self._insertFunctionStep(self.createInputFileStep)
        else:
            self._insertFunctionStep(self.ensureFastaHasNames)
        self._insertFunctionStep(self.runChaiStep)
        self._insertFunctionStep(self.extractScoreStep)
        self._insertFunctionStep(self.createOutputStep)

    def createInputFileStep(self):
        path = os.path.abspath(self._getPath())
        if not os.path.exists(path):
            os.makedirs(path)

        fastaPath = os.path.join(path, "input.fasta")
        with open(fastaPath, 'w') as f:
            counter = 1
            for line in self.inputList.get().splitlines():
                line = line.strip()
                if not line:
                    continue

                try:
                    if ')' in line:
                        json_part = line.split(')', 1)[1].strip()
                    else:
                        json_part = line
                    inpDict = json.loads(json_part)
                except json.JSONDecodeError:
                    continue

                base_name = inpDict.get("name", "unknown")
                seqFile = inpDict.get("seqFile")
                entity = inpDict.get("entity", "protein").lower()

                sequence = ""
                with open(seqFile) as sf:
                    for l in sf:
                        if l.startswith(">"):
                            continue
                        sequence += l.strip()

                uniqueName = f"{base_name}_{counter}"
                counter += 1
                f.write(f">{entity}|name={uniqueName}\n")
                f.write(f"{sequence}\n")

    def runChaiStep(self):
        if (self.inputOrigin.get() == 2 and not self.NEWFILE):
            filePath = os.path.abspath(self.file.get())
        else:
            filePath = os.path.abspath(self._getPath('input.fasta'))
        args = [str(filePath)]

        args.append(os.path.join(os.path.abspath(self._getPath()), "chai_results"))

        if self.msa.get():
            args.append("--use-msa-server")

        args.append(f" --num-trunk-recycles {self.trunkRecycles.get()}")
        args.append(f" --num-diffn-timesteps {self.timeSteps.get()}")
        args.append(f" --num-trunk-samples {self.trunkSamples.get()}")
        args.append(f" --num-diffn-samples {self.diffNsamples.get()}")

        Plugin.runCondaCommand(
            self,
            args=" ".join(args),
            condaDic=CHAI_DIC,
            program="chai-lab fold",
            cwd=os.path.abspath(Plugin.getVar(CHAI_DIC['home']))
        )

    def extractScoreStep(self):
        """Extract per-residue score and compute mean score per model"""
        resultsPath = os.path.join(os.path.abspath(self._getPath()), "chai_results")
        self.meanScore = {}

        extraFiles = self.getExtraFiles()

        for cifName in extraFiles:
            cifPath = os.path.join(resultsPath, cifName)
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
        resultsPath = os.path.join((self._getPath()), "chai_results")
        extraFiles = self.getExtraFiles()
        outputSet = SetOfAtomStructs.create(self._getPath())
        for cifName in extraFiles:
            atomStruct = AtomStruct(filename=os.path.join(resultsPath, cifName))
            atomStruct.origin = String()
            atomStruct.setAttributeValue('origin', 'Chai')
            outputSet.append(atomStruct)

        bestSrc = os.path.join(resultsPath, self.bestModel + '.cif')
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
                        modelName = os.path.splitext(os.path.basename(cifPath))[0]
                        scores[modelName] = score

        summary.append("Scores per model:")
        for modelName in sorted(scores.keys()):
            summary.append(f"  {modelName}: {scores[modelName]}")

        if scores:
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

        for name in sorted(os.listdir(resultsPath)):
            if name.lower().endswith('.cif'):
                extraFiles.append(name)

        if not extraFiles:
            raise Exception("No CIF files found in the selected folder.")

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