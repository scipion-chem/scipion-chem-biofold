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
import zipfile

import os
import re
import shutil
import pyworkflow.protocol.params as params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwem.objects import AtomStruct, SetOfAtomStructs


class ProtAlphaFold3(EMProtocol):
    """
    Protocol to import AlphaFold3 predicted structures.
    """
    _label = 'alphafold3 import'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('folder', params.FileParam,
                      label='AlphaFold3 server results: ',
                      help='Select the results folder downloaded from the AlphaFold3 server.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertStep)
        self._insertFunctionStep(self.extractPlddtStep)
        self._insertFunctionStep(self.createOutputStep)

    def convertStep(self):
        """Copy CIF files into the protocol extra folder"""
        self.extraFiles = []

        filePath = self.folder.get()
        extraPath = self._getExtraPath()
        os.makedirs(extraPath, exist_ok=True)

        with zipfile.ZipFile(filePath, 'r') as zip_ref:
            zip_ref.extractall(extraPath)

        for name in sorted(os.listdir(extraPath)):
            if name.lower().endswith('.cif'):
                self.extraFiles.append(name)

        if not self.extraFiles:
            raise Exception("No CIF files found in the selected folder.")

    def extractPlddtStep(self):
        """Extract per-residue pLDDT and compute mean pLDDT per model"""
        extraPath = self._getExtraPath()
        self.meanPlddt = {}

        for cifName in self.extraFiles:
            cifPath = os.path.join(extraPath, cifName)
            modelName = os.path.splitext(cifName)[0]

            headers = []
            plddtValues = []
            seenResidues = set()

            with open(cifPath) as f:
                for line in f:
                    if line.startswith('_atom_site.'):
                        headers.append(line.strip())
                    elif line.startswith('ATOM'):
                        break

            colIndex = {h.split('.')[-1]: i for i, h in enumerate(headers)}

            if 'B_iso_or_equiv' not in colIndex:
                raise Exception(f"No pLDDT field in {cifName}")

            with open(cifPath) as f:
                for line in f:
                    if not line.startswith('ATOM'):
                        continue
                    cols = re.sub(r'\s+', ' ', line.strip()).split()
                    resnum = int(cols[colIndex['auth_seq_id']])
                    plddt = float(cols[colIndex['B_iso_or_equiv']])
                    if resnum not in seenResidues:
                        seenResidues.add(resnum)
                        plddtValues.append(plddt)

            self.meanPlddt[modelName] = sum(plddtValues) / len(plddtValues)

        self.bestModel = max(self.meanPlddt, key=self.meanPlddt.get)

        resultsFile = self._getPath('results.txt')
        with open(resultsFile, 'w') as f:
            for model, plddt in self.meanPlddt.items():
                f.write(f"{model}\t{plddt:.2f}\n")
            f.write(f"BEST\t{self.bestModel}\n")

        self._store()

    def createOutputStep(self):
        extraPath = self._getExtraPath()
        outPath = os.path.join(extraPath, 'outputs')
        os.makedirs(outPath, exist_ok=True)

        outputSet = SetOfAtomStructs.create(self._getPath())

        for cifName in self.extraFiles:
            src = os.path.join(extraPath, cifName)
            dst = os.path.join(outPath, cifName)
            shutil.copy(src, dst)

            atomStruct = AtomStruct(filename=dst)
            outputSet.append(atomStruct)

        bestSrc = os.path.join(extraPath, self.bestModel + '.cif')
        bestStruct = AtomStruct(filename=bestSrc)

        self._defineOutputs(
            outputBestAtomStruct=bestStruct,
            outputSetOfAtomStructs=outputSet
        )

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        resultsFile = self._getPath('results.txt')

        if os.path.exists(resultsFile):
            summary.append("Mean pLDDT per model:")
            with open(resultsFile) as f:
                for line in f:
                    if line.startswith("BEST"):
                        bestModel = line.split()[1]
                    else:
                        model, plddt = line.split()
                        summary.append(f"  {model}: {plddt}")

            summary.append(f"\nBest structure (highest mean pLDDT): {bestModel}")

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