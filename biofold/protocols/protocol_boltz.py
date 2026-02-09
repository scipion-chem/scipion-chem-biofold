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
import zipfile

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

        form.addParam('entityType', params.EnumParam, default=0,
                      label='Input entity type: ', choices=['Protein', 'DNA', 'RNA'],
                      help='Input entity type to add to the set')

        form.addParam('inpChain', params.StringParam,
                      label='Input chain: ', condition='inputOrigin == 1',
                      help='Specify the protein chain to use as sequence.')

        form.addParam('inpPositions', params.StringParam,
                      label='Input positions: ',
                      help='Specify the positions of the sequence to add in the output.')

        form.addParam('cyclic', params.BooleanParam, default=False,
                      label="Cyclic: ",
                      help='Choose whether the input is cyclic or not.')

        form.addParam('addInput', params.LabelParam,
                      label='Add input: ',
                      help='Add sequence to the output set')

    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='Input')
        form.addParam('inputOrigin', params.EnumParam, default=0,
                       label='Input origin: ', choices=['Sequence', 'AtomStruct'],
                       help='Input origin to add to the set')
        self._addInputForm(form)

        form.addParam('inputList', params.TextParam, width=100,
                       default='', label='List of inputs: ',
                       help='The list of input to use for the final output set.')

        form = form.addGroup('Parameters')
        form.addParam('infPot', params.BooleanParam, default=False,
                        label="Inference potentials: ",
                        help='Choose whether to use inference potentials to improve physical plausibility of the predicted poses.')
        form.addParam('recyclingSteps', params.IntParam, default=3, expertLevel=params.LEVEL_ADVANCED,
                        label='Recycling steps: ', help="Number of recycling steps for prediction.")
        form.addParam('samplingSteps', params.IntParam, default=200,
                        label='Sampling steps: ', help="Number of sampling steps for prediction.")
        form.addParam('diffusionSamples', params.IntParam, default=1, expertLevel=params.LEVEL_ADVANCED,
                        label='Diffusion samples: ', help="Number of diffusion samples for prediction.")
        form.addParam('stepScale', params.FloatParam, default=1.638,
                        label='Steps size: ', help="Number of step size. Its related to the temperature at which the diffusion process samples the distribution.")
        form.addParam('affinityMWcorr', params.BooleanParam, default=False,
                       label="Molecular weight correction: ",
                       help='Choose whether to add the molecular weight correction to the affinity prediction.')
        form.addParam('diffusionSamplesAff', params.IntParam, default=5, expertLevel=params.LEVEL_ADVANCED,
                       label='Diffusion samples for affinity: ', help="Number of diffusion samples for affinity.")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createInputFileStep)
        self._insertFunctionStep(self.runBoltzStep)
        self._insertFunctionStep(self.createOutputStep)

    def createInputFileStep(self):
        entities = []
        chain_id_iter = iter(string.ascii_uppercase)

        for inputLine in self.inputList.get().split('\n'):
            if not inputLine.strip():
                continue

            inpJson = json.loads(inputLine.split(')')[1].strip())
            seqDic = parseFasta(os.path.abspath(inpJson['seqFile']))
            _, sequence = next(iter(seqDic.items()))
            entity = inpJson.get('entity', 'protein')
            cyclic = inpJson.get('cyclic', False)

            chain_id = next(chain_id_iter)

            entity = BoltzEntity(
                entity_type=entity,
                chain_id=chain_id,
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
            #todo ligand when all normal options work
            elif e.entity_type == "ligand":
                if e.smiles:
                    body["smiles"] = e.smiles
                elif e.ccd:
                    body["ccd"] = e.ccd

            sequences.append({e.entity_type: body})

        json_path = os.path.abspath(self._getPath("input.json"))
        yaml_path = os.path.abspath(self._getPath("input.yaml"))

        with open(json_path, "w") as f:
            json.dump({"sequences": sequences}, f, indent=2)

        script_path = os.path.join(os.path.dirname(__file__), "..", "scripts", "buildYaml.py")

        Plugin.runCondaCommand(
            self,
            program="python",
            args=f"{script_path} {json_path} {yaml_path}",
            condaDic=BOLTZ_DIC
        )


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


