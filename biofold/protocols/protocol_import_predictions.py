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
import tarfile
import zipfile

import os
import re
import shutil
import pyworkflow.protocol.params as params
from pyworkflow.object import String
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwem.objects import AtomStruct, SetOfAtomStructs


class ProtImportPredictions(EMProtocol):
    """
    Protocol to import predicted structures.
    AlphaFold3 server: https://alphafoldserver.com/
    Protenix server: https://protenix-server.com/
    Chai server: https://lab.chaidiscovery.com/
    Boltz server: https://app.tamarind.bio/boltz/
    """
    _label = 'Import structure predictions'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputOrigin', params.EnumParam, default=0,
                      label='Input origin: ', choices=['AlphaFold3', 'Protenix', 'Chai', 'Boltz'],
                      help='Input entity type to add to the set')

        form.addParam('folder', params.FileParam,
                      label='Results: ',
                      help='Select the results folder downloaded from the server.')

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

        lower = filePath.lower()

        if lower.endswith(".zip"):
            with zipfile.ZipFile(filePath, 'r') as zip_ref:
                zip_ref.extractall(extraPath)

        elif lower.endswith(".tar.gz") or lower.endswith(".tgz") or lower.endswith(".tar"):
            with tarfile.open(filePath, 'r:*') as tar_ref:
                tar_ref.extractall(extraPath)

        else:
            raise Exception("Unsupported file format. Please provide .zip or .tar.gz/.tgz archive.")


        if self.inputOrigin.get() == 1: #protenix
            for root, dirs, files in os.walk(extraPath):
                if os.path.basename(root).lower() == "predictions":
                    for name in files:
                        if name.lower().endswith(".cif"):
                            rel = os.path.relpath(os.path.join(root, name), extraPath)
                            self.extraFiles.append(rel)
        elif self.inputOrigin.get() == 0: #af3
            for root, dirs, files in os.walk(extraPath):
                if "templates" in root.lower().split(os.sep):
                    continue
                for name in files:
                    if name.lower().endswith(".cif"):
                        rel = os.path.relpath(os.path.join(root, name), extraPath)
                        self.extraFiles.append(rel)
        elif self.inputOrigin.get() == 3: #boltz
            for root, dirs, files in os.walk(extraPath):
                if os.path.basename(root).lower() == "result":
                    for name in files:
                        if name.lower().endswith(".pdb"):
                            rel = os.path.relpath(os.path.join(root, name), extraPath)
                            self.extraFiles.append(rel)

        if not self.extraFiles:
            raise Exception("No CIF/PDB files found in the selected folder.")

    def extractPlddtStep(self):
        """Extract per-residue pLDDT and compute mean pLDDT per model (supports CIF and PDB)."""
        extraPath = self._getExtraPath()
        self.meanPlddt = {}

        for fileName in self.extraFiles:
            filePath = os.path.join(extraPath, fileName)
            modelName = os.path.splitext(fileName)[0]
            plddtValues = []

            if fileName.lower().endswith(".cif"):
                headers = []
                seenResidues = set()
                with open(filePath) as f:
                    for line in f:
                        if line.startswith('_atom_site.'):
                            headers.append(line.strip())
                        elif line.startswith('ATOM'):
                            break
                colIndex = {h.split('.')[-1]: i for i, h in enumerate(headers)}
                if 'B_iso_or_equiv' not in colIndex:
                    raise Exception(f"No pLDDT field in {fileName}")
                with open(filePath) as f:
                    for line in f:
                        if not line.startswith('ATOM'):
                            continue
                        cols = re.sub(r'\s+', ' ', line.strip()).split()
                        resnum = int(cols[colIndex['auth_seq_id']])
                        if resnum in seenResidues:
                            continue
                        seenResidues.add(resnum)
                        plddt = float(cols[colIndex['B_iso_or_equiv']])
                        plddtValues.append(plddt)

            elif fileName.lower().endswith(".pdb"):
                seenResidues = set()
                with open(filePath) as f:
                    for line in f:
                        if not line.startswith("ATOM"):
                            continue
                        cols = re.sub(r'\s+', ' ', line[0:66].strip()).split()
                        try:
                            resnum = int(line[22:26].strip())
                            bfactor = float(line[60:66].strip())
                        except ValueError:
                            continue
                        if resnum in seenResidues:
                            continue
                        seenResidues.add(resnum)
                        plddtValues.append(bfactor)

            else:
                raise Exception(f"Unsupported file type for pLDDT extraction: {fileName}")

            if not plddtValues:
                raise Exception(f"No pLDDT values found in {fileName}")

            self.meanPlddt[modelName] = sum(plddtValues) / len(plddtValues)

        self.bestModel = max(self.meanPlddt, key=self.meanPlddt.get)

        # Write results to file
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

        if self.inputOrigin.get() == 0:
            origin = 'AlphaFold3'
        elif self.inputOrigin.get() == 1:
            origin = 'Protenix'
        elif self.inputOrigin.get() == 2:
            origin = 'Chai'
        else:
            origin = 'Boltz'

        for cifName in self.extraFiles:
            src = os.path.join(extraPath, cifName)

            base = os.path.basename(cifName)
            dst = os.path.join(outPath, base)

            shutil.copy(src, dst)

            atomStruct = AtomStruct(filename=dst)
            atomStruct.origin = String()
            atomStruct.setAttributeValue('origin', origin)
            outputSet.append(atomStruct)

        if origin != 'Boltz': bestSrc = os.path.join(extraPath, self.bestModel + '.cif')
        else: bestSrc = os.path.join(extraPath, self.bestModel + '.pdb')

        bestStruct = AtomStruct(filename=bestSrc)
        bestStruct.origin = String()
        bestStruct.setAttributeValue('origin', origin)

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