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

import os
import pyworkflow.protocol.params as params
from biofold.objects import BoltzEntity
from pwem.protocols import EMProtocol
from pwchem import Plugin
from biofold.constants import BOLTZ_DIC

from pwem.objects import  AtomStruct
from pwchem.protocols.Sequences.protocol_define_sequences import ProtDefineSetOfSequences
from pwchem.utils.utilsFasta import parseFasta



class ProtBoltz(EMProtocol):
    """
    Protocol to use Boltz-2 model.
    """
    _label = 'boltz-2 modelling'
    protSeq = ProtDefineSetOfSequences()

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

        form.addParam('entityType', params.EnumParam, default=0, condition='inputOrigin != 2',
                      label='Input entity type: ', choices=['Protein', 'DNA', 'RNA'],
                      help='Input entity type to add to the set')

        form.addParam('inpChain', params.StringParam,
                      label='Input chain: ', condition='inputOrigin == 1',
                      help='Specify the protein chain to use as sequence.')

        form.addParam('inpPositions', params.StringParam, condition='inputOrigin != 2',
                      label='Input positions: ',
                      help='Specify the positions of the sequence to add in the output.')

        form.addParam('cyclic', params.BooleanParam, default=False,
                      label="Cyclic: ",
                      help='Choose whether the input is cyclic or not.')

        form.addParam('addInput', params.LabelParam, condition='inputOrigin != 2',
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
                       label='Input origin: ', choices=['Sequence', 'AtomStruct', 'fasta file'],
                       help='Input origin to add to the set')
        self._addInputForm(form)

        form.addParam('inputList', params.TextParam, width=100, condition='inputOrigin in [0,1]',
                      default='', label='List of inputs: ',
                      help='The list of input to use for the final output set.')

        form.addParam('file', params.FileParam, condition='inputOrigin == 2',
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
        group.addParam('diffusionSamples', params.IntParam, default=1, expertLevel=params.LEVEL_ADVANCED,
                        label='Diffusion samples: ', help="Number of diffusion samples for prediction.")
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
        if self.inputOrigin.get() == 2:
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
                chainId=chainId,
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
        predictionsPath = os.path.join(os.path.abspath(self._getPath()), "boltz_results_input", "predictions")

        inputFolders = [f for f in os.listdir(predictionsPath) if os.path.isdir(os.path.join(predictionsPath, f))]
        if not inputFolders:
            raise Exception(f"No prediction folders found in {predictionsPath}")

        inputFolder = os.path.join(predictionsPath, inputFolders[0])

        cifFiles = sorted([f for f in os.listdir(inputFolder) if f.lower().endswith('.cif')])
        if not cifFiles:
            raise Exception(f"No CIF files found in {inputFolder}")

        cifPath = os.path.join(inputFolder, cifFiles[0])

        bestStruct = AtomStruct(filename=cifPath)

        self._defineOutputs(
            outputAtomStruct=bestStruct
        )

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
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

