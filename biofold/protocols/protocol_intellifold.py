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
import shutil
import string
import re
import os
import pyworkflow.protocol.params as params
from biofold.objects import IntelliFoldEntity
from pwem.protocols import EMProtocol
from pyworkflow.object import String

from pwchem import Plugin
from biofold.constants import INTELLIFOLD_DIC

from pwem.objects import AtomStruct, SetOfAtomStructs
from pwchem.protocols.Sequences.protocol_define_sequences import ProtDefineSetOfSequences
from pwchem.utils.utilsFasta import parseFasta
from .protocol_chai import ProtChai



class ProtIntelliFold(EMProtocol):
    """
    Protocol to use IntelliFold model.

    AI Generated:

        ProtIntelliFold - User Manual

        Overview
        --------
        The ProtIntelliFold protocol performs biomolecular structure prediction
        using the IntelliFold diffusion-based modelling framework. It supports
        proteins, DNA, and RNA inputs and generates multiple 3D structural models
        in CIF format.

        This protocol integrates sequence preprocessing, configuration file
        generation, IntelliFold execution, confidence extraction, and structured
        output creation within the Scipion framework.

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
           - Entity type (protein/DNA/RNA) is automatically inferred.

        Additional input options:
        - **entityType**: Define biological type (Protein, DNA, RNA).
        - **inpChain**: Specify chain when using AtomStruct input.
        - **inpPositions**: Select specific residue ranges.
        - **inputList**: Merge multiple inputs into a single modelling job.

        Prediction Parameters
        ---------------------
        - **recyclingSteps**: Number of recycling refinement iterations.
        - **samplingSteps**: Number of diffusion sampling steps.
        - **diffusionSamples**: Number of diffusion samples for structure prediction.
        - **msa**: Enable usage of multiple sequence alignment server.
        - **model**: Select prediction model version:
            - v1: Original IntelliFold model.
            - v2: Improved accuracy, slower inference.
            - v2-flash: Recommended default; faster and more accurate.

        Workflow
        --------
        1. **Input Processing**:
           - Parses sequences from Sequence objects, AtomStruct objects,
             or FASTA files.
           - Automatically infers entity type when required.
           - Generates a structured `input.json` file describing all entities.

        2. **Configuration Generation**:
           - Converts `input.json` into `input.yaml` required by IntelliFold.

        3. **Model Execution**:
           - Executes `intellifold predict` inside the configured Conda environment.
           - Applies selected recycling, sampling, diffusion, and model parameters.
           - Optionally connects to the MSA server.
           - Stores predictions in the protocol working directory.

        4. **Confidence Extraction and Ranking**:
           - Reads ranking scores and confidence metrics from
             `_summary_confidences.json` files.
           - Extracts:
               - ranking_score
               - pLDDT
               - pTM
           - Computes model ranking and selects the best predicted structure.
           - Generates a `results.txt` summary file.

        5. **Output Creation**:
           - Copies all predicted CIF files into the protocol output folder.
           - Wraps them into a SetOfAtomStructs object.
           - Identifies and returns the best predicted structure separately.

        Outputs
        -------
        - **outputBestAtomStruct**: Best predicted 3D structure in CIF format
          (highest ranking score).
        - **outputSetOfAtomStructs**: Set containing all predicted structures.
        - **results.txt**: Tabulated summary including ranking score, pLDDT, and pTM.

        Practical Recommendations
        -------------------------
        - Use the 'v2-flash' model for optimal balance between speed and accuracy.
        - Increase recycling steps for improved structural refinement.
        - Adjust diffusion samples depending on system complexity.
        - Enable MSA for improved prediction quality when available.

        Interpretation
        --------------
        The resulting AtomStruct corresponds to the top-ranked IntelliFold
        prediction according to the ranking score metric.

        The output structures can be used for:

        - Structural visualization
        - Docking studies
        - Interaction analysis
        - Molecular dynamics simulations
        - Comparative structural modelling

        Warnings
        --------
        - Large sequences may require significant computational resources.
        - Increasing diffusion samples substantially increases runtime.
        - Ensure input sequences are biologically valid and correctly formatted.

        Final Perspective
        -----------------
        ProtIntelliFold provides an integrated interface for IntelliFold-based
        biomolecular modelling within Scipion. It automates sequence preparation,
        configuration building, model execution, confidence evaluation, and
        structured result retrieval, enabling reproducible structural prediction
        workflows.
    """
    _label = 'intellifold modelling'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addHidden('useGpu', params.BooleanParam, default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation. Choose one.")

        form.addHidden('gpuList', params.StringParam, default='0',
                       label="Choose GPU IDs",
                       help="Comma-separated GPU devices that can be used.")

        form.addSection(label='Input')
        form.addParam('inputOrigin', params.EnumParam, default=0,
                      label='Input origin: ', choices=['Sequence', 'SetOfSequences', 'AtomStruct', 'fasta file'],
                       help='Input origin to add to the set')
        ProtChai._addInputForm(form)

        form.addParam('inputList', params.TextParam, width=100, condition='inputOrigin!=3',
                      default='', label='List of inputs: ',
                      help='The list of input to use for the final output set.')

        form.addParam('file', params.FileParam, condition='inputOrigin == 3',
                      label='Sequence file: ',
                      help='Select the results fasta file.')

        group = form.addGroup('Parameters')
        group.addParam('recyclingSteps', params.IntParam, default=10, expertLevel=params.LEVEL_ADVANCED,
                        label='Recycling steps: ', help="Number of recycling steps for prediction.")
        group.addParam('samplingSteps', params.IntParam, default=200,
                        label='Sampling steps: ', help="Number of sampling steps for prediction.")
        group.addParam('diffusionSamples', params.IntParam, default=5,
                        label='Output models: ', help="Number of output models for prediction.")
        group.addParam('msa', params.BooleanParam, default=True,
                       label="Use MSA: ",
                       help='Choose whether to use multiple sequence alignment.')
        group.addParam('model', params.EnumParam, default=2,
                      label='Prediction model: ', choices=['v1', 'v2', 'v2-flash'],
                      help="Options are 'v1', 'v2', and 'v2-flash'. 'v2-flash' is the default and recommended model, "
                           "which is faster and more accurate than 'v1' and 'v2'. 'v1' is the original model used in "
                           "the IntelliFold paper, and 'v2' is an improved version of the model with better performance "
                           "but slower inference speed than 'v2-flash'. You can choose the model based on your needs and computational resources.")


        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        if self.inputOrigin.get() == 3:
            self._insertFunctionStep(self.createJsonFromFastaStep)
        else:
            self._insertFunctionStep(self.createInputFileStep)
        self._insertFunctionStep(self.createYamlFileStep)
        self._insertFunctionStep(self.runIntelliFoldStep)
        self._insertFunctionStep(self.createOutputStep)

    def createYamlFileStep(self):
        jsonPath = os.path.abspath(self._getPath("input.json"))
        yamlPath = os.path.abspath(self._getPath("input.yaml"))

        scriptPath = os.path.join(os.path.dirname(__file__), "..", "scripts", "buildYaml.py")

        Plugin.runCondaCommand(
            self,
            program="python",
            args=f"{scriptPath} {jsonPath} {yamlPath}",
            condaDic=INTELLIFOLD_DIC
        )

    def createJsonFromFastaStep(self):
        fastaPath = os.path.abspath(self.file.get())

        entities = []
        nextChain = iter(string.ascii_uppercase)

        currentHeader = None
        sequenceBuffer = ""

        with open(fastaPath) as f:
            for line in f:
                line = line.strip()

                if line.startswith(">"):

                    if currentHeader is not None:

                        entityType = self.guessEntityType(sequenceBuffer)

                        chainIds = self.extractChainsFromHeader(currentHeader)

                        if not chainIds:
                            chainIds = [next(nextChain)]

                        entityDict = {
                            "id": chainIds if len(chainIds) > 1 else chainIds[0],
                            "sequence": sequenceBuffer
                        }

                        entities.append({entityType: entityDict})

                    currentHeader = line[1:]
                    sequenceBuffer = ""

                else:
                    sequenceBuffer += line

        if currentHeader is not None:

            entityType = self.guessEntityType(sequenceBuffer)

            chainIds = self.extractChainsFromHeader(currentHeader)

            if not chainIds:
                chainIds = [next(nextChain)]

            entityDict = {
                "id": chainIds if len(chainIds) > 1 else chainIds[0],
                "sequence": sequenceBuffer
            }

            entities.append({entityType: entityDict})

        jsonPath = os.path.abspath(self._getPath("input.json"))

        with open(jsonPath, "w") as f:
            json.dump({"sequences": entities}, f, indent=2)


    def createInputFileStep(self):
        entities = []
        chainIdIiter = iter(string.ascii_uppercase)

        for inputLine in self.inputList.get().split('\n'):
            if not inputLine.strip():
                continue

            inpJson = json.loads(inputLine.split(')')[1].strip())
            seqDic = parseFasta(os.path.abspath(inpJson['seqFile']))
            _, sequence = next(iter(seqDic.items()))
            entity = inpJson.get('entity', 'protein')

            chainId = next(chainIdIiter)

            entity = IntelliFoldEntity(
                entity_type=entity,
                chain_id=chainId,
                sequence=sequence
            )
            entities.append(entity)

        merged = {}
        for e in entities:
            key = (
                e.entity_type,
                e.sequence,
                e.smiles,
                e.ccd,
                e.msa,
            )
            if key not in merged:
                merged[key] = e
            else:
                merged[key].ids.extend(e.ids)

        sequences = []
        for e in merged.values():
            body = {
                "id": e.ids if len(e.ids) > 1 else e.ids[0]
            }

            if e.entity_type in ("protein", "dna", "rna"):
                body["sequence"] = e.sequence
                if e.msa:
                    body["msa"] = "empty"

            sequences.append({e.entity_type: body})

        jsonPath = os.path.abspath(self._getPath("input.json"))

        with open(jsonPath, "w") as f:
            json.dump({"sequences": sequences}, f, indent=2)



    def runIntelliFoldStep(self):
        filePath = os.path.abspath(self._getPath("input.yaml"))
        args = [str(filePath)]

        if self.msa.get():
            args.append("--use_msa_server")

        args.extend([
            "--recycling_iters", str(self.recyclingSteps.get()),
            "--sampling_steps", str(self.samplingSteps.get()),
            "--num_diffusion_samples", str(self.diffusionSamples.get()),
            "--model", self.getEnumText('model'),
            "--num_workers", str(self.numberOfThreads.get()),
            "--out_dir", os.path.abspath(self._getPath('intellifold_output'))
        ])

        program = "intellifold predict"
        if self.useGpu.get():
            program = f"CUDA_VISIBLE_DEVICES={self.gpuList.get()} intellifold predict"

        Plugin.runCondaCommand(
            self,
            args=" ".join(args),
            condaDic=INTELLIFOLD_DIC,
            program=program,
            cwd=os.path.abspath(Plugin.getVar(INTELLIFOLD_DIC['home']))
        )

    def createOutputStep(self):
        resultsPath = os.path.join(self._getPath(), "intellifold_output/input/predictions/input")
        if not os.path.exists(resultsPath):
            raise Exception(f"Predictions folder does not exist: {resultsPath}")

        cifFiles = sorted([f for f in os.listdir(resultsPath) if f.lower().endswith('.cif')])
        if not cifFiles:
            raise Exception(f"No CIF files found in {resultsPath}")

        outPath = os.path.join(self._getExtraPath(), 'outputs')
        os.makedirs(outPath, exist_ok=True)

        outputSet = SetOfAtomStructs.create(self._getPath())
        self.meanConf = {}

        for cifName in cifFiles:
            cifPath = os.path.join(resultsPath, cifName)
            summaryFile = cifPath.replace('.cif', '_summary_confidences.json')
            if not os.path.exists(summaryFile):
                raise Exception(f"Missing summary confidence file for {cifName}")

            with open(summaryFile) as f:
                summary = json.load(f)

            self.meanConf[cifName] = summary.get('ranking_score', 0)

            dst = os.path.join(outPath, cifName)
            shutil.copy(cifPath, dst)
            atomStruct = AtomStruct(filename=dst)
            atomStruct.origin = String()
            atomStruct.setAttributeValue('origin', 'IntelliFold')
            outputSet.append(atomStruct)

        self.bestModel = max(self.meanConf, key=self.meanConf.get)
        bestStructPath = os.path.join(outPath, self.bestModel)
        bestStruct = AtomStruct(filename=bestStructPath)
        bestStruct.origin = String()
        bestStruct.setAttributeValue('origin', 'IntelliFold')

        resultsFile = os.path.join(self._getPath(), 'results.txt')
        with open(resultsFile, 'w') as f:
            f.write("Model\tRankingScore\tpLDDT\tPTM\n")
            for cifName in sorted(self.meanConf.keys()):
                summaryFile = os.path.join(resultsPath, cifName.replace('.cif', '_summary_confidences.json'))
                with open(summaryFile) as sf:
                    summary = json.load(sf)
                ranking_score = summary.get('ranking_score', 0)
                plddt = summary.get('plddt', 0)
                ptm = summary.get('ptm', 0)
                f.write(f"{cifName}\t{ranking_score:.3f}\t{plddt:.3f}\t{ptm:.3f}\n")

            f.write(f"BEST\t{self.bestModel}\t{self.meanConf[self.bestModel]:.3f}\n")

        self._defineOutputs(
            outputBestAtomStruct=bestStruct,
            outputSetOfAtomStructs=outputSet
        )

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        resultsFile = os.path.join(self._getPath(), 'results.txt')
        if not os.path.exists(resultsFile):
            return ["Results file does not exist."]

        summary = ["IntelliFold predictions summary:"]
        with open(resultsFile) as f:
            for line in f:
                summary.append(line.strip())

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
    def guessEntityType(self, sequence):
        """ Guess if a sequence is DNA, RNA or Protein """
        sequence = sequence.upper().strip()

        protein_only = re.compile(r'[DEFHIKLMPQRVWY]')

        if protein_only.search(sequence):
            return 'protein'

        if 'U' in sequence:
            return 'rna'

        return 'dna'

    def extractChainsFromHeader(self, header):
        """
        Extract chain IDs from RCSB FASTA headers.

        Examples
        --------
        >8ZB4_1|Chains A, B|...
            -> ['A', 'B']

        >1ABC_1|Chain C|...
            -> ['C']
        """
        m = re.search(r'\bChains?\s+([^|]+)', header)

        if not m:
            return []

        return [c.strip() for c in m.group(1).split(",") if c.strip()]
