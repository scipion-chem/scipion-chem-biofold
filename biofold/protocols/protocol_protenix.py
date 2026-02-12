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
from biofold import PROTENIX_DIC
from biofold.objects import BoltzEntity
from pwem.protocols import EMProtocol
from pwchem import Plugin
from biofold.constants import BOLTZ_DIC

from pwem.objects import  AtomStruct
from pwchem.protocols.Sequences.protocol_define_sequences import ProtDefineSetOfSequences
from pwchem.utils.utilsFasta import parseFasta



class ProtProtenix(EMProtocol):
    """
    Protocol to use Protenix model.
    """
    _label = 'protenix modelling'
    protSeq = ProtDefineSetOfSequences()
    models = ['protenix_base_default_v1.0.0', 'protenix_base_20250630_v1.0.0', 'protenix_base_default_v0.5.0', 'protenix_mini_default_v0.5.0']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='Input')
        form.addParam('inputOrigin', params.EnumParam, default=0,
                       label='Input origin: ', choices=['AtomStruct', 'pdb file'],
                       help='Input origin of input.')
        form.addParam('inputAtomStruct', params.PointerParam,
                      pointerClass='AtomStruct',
                      label="Input structure: ", condition='inputOrigin==0',
                      help='Select the AtomStruct object whose sequence to add to the set')
        form.addParam('file', params.FileParam, condition='inputOrigin == 1',
                      label='PDB file: ',
                      help='Select the PDB file.')

        form.addParam('msa', params.BooleanParam, default=True,
                        label="Use MSA: ",
                        help='Choose whether to use MSA.')
        form.addParam('model', params.EnumParam, default=0,
                      label='Model: ', choices=self.models,
                      help='Model to use for prediction. \nprotenix_base_default_v1.0.0 : trained with data cutoff (2021-09-30) aligned with AF3. \nprotenix_base_20250630_v1.0.0 : trained with an updated data cutoff (2025-06-30) for better practical performance. \nprotenix_base_default_v0.5.0 : previous version of the model. \nprotenix_mini_default_v0.5.0 : fast and lightweight.')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createJSONStep)
        #self._insertFunctionStep(self.runProtenixStep)
        #self._insertFunctionStep(self.createOutputStep)

    def createJSONStep(self):
        input = self.inputAtomStruct.get().getFileName()
        output = self._getExtraPath()

        Plugin.runCondaCommand(
            self,
            program="protenix json",
            args=f" --input {os.path.abspath(input)} --out_dir {os.path.abspath(output)} --altloc first",
            condaDic=PROTENIX_DIC
        )

    def runProtenixStep(self):
        filePath = os.path.abspath(self._getPath("input.yaml"))
        args = [str(filePath)]


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

