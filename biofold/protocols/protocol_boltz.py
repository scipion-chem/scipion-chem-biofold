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
from pwchem import Plugin
from biofold.constants import BOLTZ_DIC

from pwem.objects import  AtomStruct, SetOfAtomStructs


class ProtBoltz(EMProtocol):
    """
    Protocol to use Boltz-2 model.
    """
    _label = 'boltz-2 modelling'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        #todo build my own files to use as .yaml
        form.addSection(label='Input')

        form.addParam('oneFile', params.BooleanParam, default=True, label='One single file: ')
        form.addParam('file', params.FileParam, condition='oneFile',
                      label='Sequence file: ',
                      help='Select the results folder downloaded from the AlphaFold3 server.')
        form.addParam('filesPath', params.PathParam, condition='not oneFile',
                       label="Sequences files directory: ",
                       help="Directory with the files you want to import.\n\n"
                            "Files should be .yaml")

        form = form.addGroup('Parameters')
        form.addParam('infPot', params.BooleanParam, default=False,
                        label="Inference potentials: ",
                        help='Choose whether to use inference potentials to improve physical plausability of the predicted poses.')
        form.addParam('recyclingSteps', params.IntParam, default=3, expertLevel=params.LEVEL_ADVANCED,
                        label='Recycling steps: ', help="Number of recycling steps for prediction.")
        form.addParam('samplingSteps', params.IntParam, default=200,
                        label='Sampling steps: ', help="Number of sampling steps for prediction.")
        form.addParam('diffusionSamples', params.IntParam, default=1, expertLevel=params.LEVEL_ADVANCED,
                        label='Difussion samples: ', help="Number of diffusion samples for prediction.")
        form.addParam('stepScale', params.FloatParam, default=1.638,
                        label='Steps size: ', help="Number of step size. Its related to the temperature at which the diffusion process samples the distribution.")
        form.addParam('affinityMWcorr', params.BooleanParam, default=False,
                       label="Molecular weight correction: ",
                       help='Choose whether to add the molecular weight correction to the affinity prediction.')
        form.addParam('diffusionSamplesAff', params.IntParam, default=5, expertLevel=params.LEVEL_ADVANCED,
                       label='Difussion samples for affinity: ', help="Number of diffusion samples for affinity.")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.runBoltzStep)
        #self._insertFunctionStep(self.extractPlddtStep)
        #self._insertFunctionStep(self.createOutputStep)

    def runBoltzStep(self):
        if self.oneFile.get():
            filePath = os.path.abspath(self.file.get())
            args = [str(filePath)]
        else:
            filePath = os.path.abspath(self.filesPath.get())
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
                raise Exception(f"No score field in {cifName}")

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

    def createOutputStep(self):
        """Define Scipion outputs"""
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
            outputSetOfAtomStructs=outputSet,
            outputBestAtomStruct=bestStruct
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