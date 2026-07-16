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
import string
import re

import os
import pyworkflow.protocol.params as params
from biofold.objects import BoltzEntity
from pwem.protocols import EMProtocol
from pyworkflow.object import String
import shutil

from pwchem import Plugin
from biofold.constants import BOLTZ_DIC

from pwem.objects import  AtomStruct, SetOfAtomStructs
from pwem.objects import Sequence, SetOfSequences
from pwchem.protocols.Sequences.protocol_define_sequences import ProtDefineSetOfSequences
from pwchem.utils.utilsFasta import parseFasta



class ProtBoltz(EMProtocol):
    """
    Protocol to use Boltz-2 model.

    AI Generated:

    ProtBoltz - User Manual

    Overview
    --------
    The ProtBoltz protocol performs biomolecular structure prediction using the
    Boltz-2 diffusion-based modelling framework. It supports proteins, DNA, and RNA
    inputs and can optionally estimate binding affinity. The protocol generates
    3D structural models from sequence information and returns the best predicted
    structure as an AtomStruct object.

    This protocol integrates input preprocessing, configuration file generation,
    model execution, and structured output creation within the Scipion framework.

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
    - **cyclic**: Specify whether the molecule is cyclic.
    - **inputList**: Multiple inputs merged into a single modelling job.
    - **GPU usage**: Enable GPU acceleration and specify GPU IDs.

    Prediction Parameters
    ---------------------
    - **infPot**: Enable inference potentials to improve physical plausibility.
    - **recyclingSteps**: Number of recycling refinement steps.
    - **samplingSteps**: Number of diffusion sampling steps.
    - **diffusionSamples**: Number of diffusion samples for structure prediction.
    - **stepScale**: Step size scaling factor (controls diffusion temperature).
    - **affinityMWcorr**: Apply molecular weight correction to affinity prediction.
    - **diffusionSamplesAff**: Number of diffusion samples for affinity estimation.

    Workflow
    --------
    1. **Input Processing**:
       - Parses sequences from Sequence objects, AtomStruct objects, or FASTA files.
       - Automatically infers entity type if not specified.
       - Merges identical entities when applicable.
       - Generates a structured `input.json` file.

    2. **Configuration Generation**:
       - Converts `input.json` into `input.yaml` required by Boltz-2.

    3. **Model Execution**:
       - Executes `boltz predict` using CPU or GPU.
       - Applies selected diffusion, recycling, and affinity parameters.
       - Stores results in the protocol working directory.

    4. **Output Creation**:
       - Retrieves the best predicted CIF structure.
       - Wraps it into an AtomStruct object.

    Outputs
    -------
    - **AtomStruct**: Best predicted 3D structure in CIF format.
    - Prediction results stored in the protocol output directory.

    Practical Recommendations
    -------------------------
    - Use GPU acceleration for faster predictions when available.
    - Increase recycling steps for improved structural refinement.
    - Adjust sampling and diffusion parameters depending on system size.
    - Use inference potentials for better physical realism in complex systems.

    Interpretation
    --------------
    The resulting AtomStruct corresponds to the top-ranked Boltz-2 prediction.
    The output can be used for:

    - Structural visualization
    - Downstream docking
    - Interaction energy calculations
    - Molecular dynamics simulations
    - Comparative structural analysis

    Warnings
    --------
    - Large sequences may require substantial GPU memory.
    - Multiple entities increase computational cost.
    - Ensure input sequences are biologically valid and correctly formatted.

    Final Perspective
    -----------------
    ProtBoltz provides an integrated interface for advanced diffusion-based
    biomolecular modelling within Scipion. It streamlines sequence preparation,
    configuration building, execution, and structure retrieval, enabling
    reproducible structural prediction workflows.
    """
    _label = 'boltz-2 modelling'

    # -------------------------- DEFINE param functions ----------------------
    def _addInputForm(self, form):
        form.addParam('inputSequence', params.PointerParam,
                      pointerClass='Sequence', allowsNull=True,
                      label="Input sequence: ", condition='inputOrigin==0',
                      help='Select the sequence object to add to the set')

        form.addParam('inputSetOfSequences', params.PointerParam,
                      pointerClass='SetOfSequences', allowsNull=True,
                      label="Input set of sequences: ", condition='inputOrigin==1',
                      help='Select the sequence set to select the sequence from.')
        form.addParam('inputSequenceFromSet', params.StringParam, allowsNull=True,
                      label="Input sequence: ", condition='inputOrigin==1',
                      help='Select the sequence object to add to the set')

        form.addParam('inputAtomStruct', params.PointerParam,
                      pointerClass='AtomStruct', allowsNull=True,
                      label="Input structure: ", condition='inputOrigin==2',
                      help='Select the AtomStruct object whose sequence to add to the set')

        form.addParam('entityType', params.EnumParam, default=0, condition='inputOrigin != 3',
                      label='Input entity type: ', choices=['Protein', 'DNA', 'RNA'],
                      help='Input entity type to add to the set')

        form.addParam('inpChain', params.StringParam,
                      label='Input chain: ', condition='inputOrigin == 2',
                      help='Specify the protein chain to use as sequence.')

        form.addParam('inpPositions', params.StringParam, condition='inputOrigin != 3',
                      label='Input positions: ',
                      help='Specify the positions of the sequence to add in the output.')

        form.addParam('cyclic', params.BooleanParam, default=False,
                      label="Cyclic: ",
                      help='Choose whether the input is cyclic or not.')

        form.addParam('addInput', params.LabelParam, condition='inputOrigin != 3',
                      label='Add input: ',
                      help='Add sequence to the output set')

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
                       label='Input origin: ', choices=['Sequence','SetOfSequences','AtomStruct', 'fasta file'],
                       help='Input origin to add to the set')
        self._addInputForm(form)

        form.addParam('inputList', params.TextParam, width=100, condition='inputOrigin!=3',
                      default='', label='List of inputs: ',
                      help='The list of input to use for the final output set.')

        form.addParam('file', params.FileParam, condition='inputOrigin == 3',
                      label='Sequence file: ',
                      help='Select the results fasta file.')

        group = form.addGroup('Parameters')
        group.addParam('infPot', params.BooleanParam, default=False,
                        label="Inference potentials: ",
                        help='Choose whether to use inference potentials to improve physical plausibility of the predicted poses.')
        group.addParam('recyclingSteps', params.IntParam, default=3, expertLevel=params.LEVEL_ADVANCED,
                        label='Recycling steps: ', help="Number of recycling steps for prediction.")
        group.addParam('samplingSteps', params.IntParam, default=200,
                        label='Sampling steps: ', help="Number of sampling steps for prediction.")
        group.addParam('diffusionSamples', params.IntParam, default=5,
                        label='Output models: ', help="Number of output models for prediction.")
        group.addParam('stepScale', params.FloatParam, default=1.638,
                        label='Steps size: ', help="Number of step size. Its related to the temperature at which the diffusion process samples the distribution.")
        group.addParam('affinityMWcorr', params.BooleanParam, default=False,
                       label="Molecular weight correction: ",
                       help='Choose whether to add the molecular weight correction to the affinity prediction.')
        group.addParam('diffusionSamplesAff', params.IntParam, default=5, expertLevel=params.LEVEL_ADVANCED,
                       label='Diffusion samples for affinity: ', help="Number of diffusion samples for affinity.")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        if self.inputOrigin.get() == 3:
            self._insertFunctionStep(self.createJsonFromFastaStep)
        else:
            self._insertFunctionStep(self.createInputFileStep)
        self._insertFunctionStep(self.createYamlFileStep)
        self._insertFunctionStep(self.runBoltzStep)
        self._insertFunctionStep(self.createOutputStep)

    def createYamlFileStep(self):
        jsonPath = os.path.abspath(self._getPath("input.json"))
        yamlPath = os.path.abspath(self._getPath("input.yaml"))

        scriptPath = os.path.join(os.path.dirname(__file__), "..", "scripts", "buildYaml.py")

        Plugin.runCondaCommand(
            self,
            program="python",
            args=f"{scriptPath} {jsonPath} {yamlPath}",
            condaDic=BOLTZ_DIC
        )

    def createJsonFromFastaStep(self):
        fastaPath = os.path.abspath(self.file.get())
        seqDic = parseFasta(fastaPath)

        chainIdIiter = iter(string.ascii_uppercase)
        entities = []

        for seqName, sequence in seqDic.items():
            chainId = next(chainIdIiter)
            entity = self.guessEntityType(sequence)
            cyclic = self.cyclic.get()

            entityDict = {
                "id": chainId,
                "cyclic": cyclic,
                "sequence": sequence
            }

            entities.append({entity: entityDict})

        jsonPath = os.path.abspath(self._getPath("input.json"))
        with open(jsonPath, 'w') as f:
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
            cyclic = inpJson.get('cyclic', False)

            chainId = next(chainIdIiter)

            entity = BoltzEntity(
                entity_type=entity,
                chain_id=chainId,
                sequence=sequence,
                cyclic=cyclic
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
                e.cyclic
            )
            if key not in merged:
                merged[key] = e
            else:
                merged[key].ids.extend(e.ids)

        sequences = []
        for e in merged.values():
            body = {
                "id": e.ids if len(e.ids) > 1 else e.ids[0],
                "cyclic": e.cyclic
            }

            if e.entity_type in ("protein", "dna", "rna"):
                body["sequence"] = e.sequence
                if e.msa:
                    body["msa"] = e.msa

            sequences.append({e.entity_type: body})

        jsonPath = os.path.abspath(self._getPath("input.json"))

        with open(jsonPath, "w") as f:
            json.dump({"sequences": sequences}, f, indent=2)



    def runBoltzStep(self):
        filePath = os.path.abspath(self._getPath("input.yaml"))
        args = [str(filePath)]

        if self.infPot.get():
            args.append("--use_potentials")

        args.append("--use_msa_server --cache ./mol")
        args.append(f" --recycling_steps {self.recyclingSteps.get()}")
        args.append(f" --sampling_steps {self.samplingSteps.get()}")
        args.append(f" --diffusion_samples {self.diffusionSamples.get()}")
        args.append(f" --step_scale {self.stepScale.get()}")

        if self.affinityMWcorr.get():
            args.append(" --affinity_mw_correction")

        args.append(f" --diffusion_samples_affinity {self.diffusionSamplesAff.get()}")
        args.append(f" --out_dir {os.path.abspath(self._getPath())}")

        if self.useGpu.get():
            gpu_ids = self.gpuList.get()
            selected_gpu = gpu_ids.split(",")[0]
            os.environ["CUDA_VISIBLE_DEVICES"] = selected_gpu
            args.append("--accelerator gpu")
        else:
            args.append("--accelerator cpu")

        Plugin.runCondaCommand(
            self,
            args=" ".join(args),
            condaDic=BOLTZ_DIC,
            program="boltz predict",
            cwd=os.path.abspath(Plugin.getVar(BOLTZ_DIC['home']))
        )

    def createOutputStep(self):
        predictionsPath = os.path.join(
            os.path.abspath(self._getPath()),
            "boltz_results_input",
            "predictions"
        )

        inputFolders = [
            f for f in os.listdir(predictionsPath)
            if os.path.isdir(os.path.join(predictionsPath, f))
        ]
        if not inputFolders:
            raise Exception(f"No prediction folders found in {predictionsPath}")

        inputFolder = os.path.join(predictionsPath, inputFolders[0])

        cifFiles = sorted([f for f in os.listdir(inputFolder) if f.lower().endswith('.cif')])
        if not cifFiles:
            raise Exception(f"No CIF files found in {inputFolder}")

        # Create output directory
        outPath = os.path.join(self._getExtraPath(), 'outputs')
        os.makedirs(outPath, exist_ok=True)

        # Create set of structures
        outputSet = SetOfAtomStructs.create(self._getPath())

        scores = {}

        for cifName in cifFiles:
            cifPath = os.path.join(inputFolder, cifName)

            modelBase = os.path.splitext(cifName)[0]
            jsonName = f"confidence_{modelBase}.json"
            jsonFile = os.path.join(inputFolder, jsonName)
            score = 0

            if os.path.exists(jsonFile):
                try:
                    with open(jsonFile) as f:
                        data = json.load(f)
                    score = data.get("confidence_score", 0)
                except:
                    pass

            scores[cifName] = score

            # Copy to output folder
            dst = os.path.join(outPath, cifName)
            shutil.copy(cifPath, dst)

            atomStruct = AtomStruct(filename=dst)
            atomStruct.origin = String()
            atomStruct.setAttributeValue('origin', 'Boltz')

            outputSet.append(atomStruct)

        # Select best model
        bestModel = max(scores, key=scores.get)
        bestStructPath = os.path.join(outPath, bestModel)

        bestStruct = AtomStruct(filename=bestStructPath)
        bestStruct.origin = String()
        bestStruct.setAttributeValue('origin', 'Boltz')

        # Optional: write summary file
        resultsFile = os.path.join(self._getPath(), 'results.txt')
        with open(resultsFile, 'w') as f:
            f.write("Model\tConfidenceScore\n")
            for name in sorted(scores.keys()):
                f.write(f"{name}\t{scores[name]:.3f}\n")
            f.write(f"BEST\t{bestModel}\t{scores[bestModel]:.3f}\n")

        self._defineOutputs(
            outputBestAtomStruct=bestStruct,
            outputSetOfAtomStructs=outputSet
        )

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        resultsFile = os.path.join(self._getPath(), 'results.txt')

        if not os.path.exists(resultsFile):
            return ["Results file does not exist."]

        summary = ["Boltz predictions summary:"]

        try:
            with open(resultsFile) as f:
                for line in f:
                    summary.append(line.strip())
        except Exception as e:
            summary.append(f"Error reading results file: {e}")

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

