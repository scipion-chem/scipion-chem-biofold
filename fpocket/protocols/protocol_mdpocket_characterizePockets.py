# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Lobna Ramadane Morchadi (lobna.ramadane@alumnos.upm.es)
# *          Blanca Pueche (blanca.pueche@cnb.csic.es)
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
This protocol is used to perform a pocket charcaterization on a protein structure using the FPocket software

"""

import os, shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfStructROIs, StructROI
from gromacs import Plugin
from pwchem.utils import *
from fpocket import Plugin
from fpocket.constants import *


class MDpocketCharacterize(EMProtocol):
    """
    Executes the mdpocket software to perform pocket characterization.
    """
    _label = 'MDPocket pocket characterization'
    _pocketTypes = ['Small molecule binding sites', 'Putative channels and small cavities', 'Water binding sites', 'Big external pockets']
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('useSystem', params.BooleanParam, deafult=True,
                      label='Use MD system as input: ',
                      help='Select input files, Yes = MD System, No = set of ROIs')
        form.addParam('inputSystem', params.PointerParam, condition='useSystem',
                       pointerClass='GromacsSystem', allowsNull=True,
                       label="Input gromacs system: ",
                       help='Select the MD system to search for pockets.')
        form.addParam('inputPDBs', params.PointerParam, condition='not useSystem',
                      pointerClass='SetOfStructROIs', allowsNull=True,
                      label="Input set of struct ROIs: ",
                      help='Select the structural ROIs to use as input.')

        #todo change and use only one ROI
        form.addParam('inputROIs', params.PointerParam,
                      pointerClass='SetOfStructROIs',
                      label="Input ROIs: ",
                      help='Select the ROI(s) to use as input for characterization.')
        form.addParam('pocket', params.StringParam,
                      label="Input specific ROI: ",
                      help='Select the specific ROI to use as input for characterization.')

        form.addParallelSection(threads=4)

    def _getMDpocketSystemArgs(self, pocket):
        pdbFile = self.moveFiles()

        trajFile = self.inputSystem.get().getTrajectoryFile()
        trajBasename = os.path.basename((trajFile))
        args = ['--trajectory_file', trajBasename]

        trajExt = os.path.splitext(trajFile)[1][1:]
        args.append('--trajectory_format')
        args.append(trajExt)

        args.append('-f')
        args.append(pdbFile)

        args.append('--selected_pocket')
        p = os.path.basename(pocket.getFileName())
        args.append(p)

        return args

    def _getMDpocketPDBsArgs(self, pocket2):
        routes = []
        for pocket in self.inputPDBs.get():
            if pocket.getClass() == 'FPocket':
                routes.append(os.path.abspath(str(pocket.getFileName())))
            else:
                routes.append(os.path.abspath(str(pocket._extraFile)))

        inputFile = self._getExtraPath("mdpocketInputFile.txt")
        with open(inputFile, 'w') as f:
            for pdbRoute in routes:
                f.write(pdbRoute + '\n')


        movedFile = self.moveFilePDB()
        pluginArgs = ["-L", movedFile]

        pluginArgs.append('--selected_pocket')
        p = os.path.basename(pocket2.getFileName())
        pluginArgs.append(p)

        return pluginArgs

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('mdPocketStep')
        self._insertFunctionStep('defineOutputStep')

    def mdPocketStep(self):
        specificPocket = self.getSpecifiedPocketFile()
        if self.useSystem.get():
            # use system
            #if isinstance(self.pockets.get(), StructROI):
                self.runMDPocket(specificPocket)
            #else:
            #    for pocket in self.pockets.get():
            #        self.runMDPocket(pocket)
        else:
            # use pdb files
            #if isinstance(self.pockets.get(), StructROI):
                self.runMDPocketPDB(specificPocket)
            #else:
            #    for pocket in self.pockets.get():
            #        self.runMDPocketPDB(pocket)

        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        self.cleanUp(mdpocketDir)


    def defineOutputStep(self):
        pocketsDir = self._getExtraPath()
        pocketFiles = os.listdir(pocketsDir)

        outPockets = SetOfStructROIs(filename=self._getPath('pockets.sqlite'))
        proteinFile = self._getExtraPath("mdpout_all_atom_pdensities.pdb")
        for pFile in pocketFiles:
            if '.pdb' in pFile:
                roi = StructROI(filename=os.path.join(pocketsDir, pFile), proteinFile=proteinFile, pClass='MDpocket')
                outPockets.append(roi)

        outPockets.buildPDBhetatmFile()
        self._defineOutputs(outputSet=outPockets)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = ["If different isovalues were selected, the new pdb files appear in the 'extra' folder."]
        return summary

    def _methods(self):
        methods = []
        return methods

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def getCoords(self):
        pdbFile = os.path.abspath(self._getExtraPath('mdpoutput-{}.pdb'.format(str(self.isoValue.get()))))
        coords = []
        for pdbLine in open(pdbFile):
            line = pdbLine.split()
            coord = line[5:8]
            coord = list(map(float,coord))
            coords.append(coord)
        return coords

    def createPocketFile(self, clust, i):
        outFile = self._getExtraPath('pocketFile_{}.pdb'.format(i+1))
        with open(outFile, 'w') as f:
            for j, coord in enumerate(clust):
                f.write(writePDBLine(['HETATM', str(j + 1), 'APOL', 'STP', 'C', '1', *coord, 1.0, 0.0, '', 'Ve']))

            f.write('\nEND')
        return outFile

    def moveFiles(self):
        # move files to path where mdpocket is, it is picky with where it is executed and they input files routes
        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        trajFile = self.inputSystem.get().getTrajectoryFile()
        trajectory = os.path.abspath((trajFile))
        pdbFile = self.inputSystem.get().getAttributeValue('pdbFile')
        #for p in self.pockets.get():
        specificPocket = self.getSpecifiedPocketFile()
        shutil.copy(str(specificPocket), os.path.join(mdpocketDir, os.path.basename(specificPocket)))
        shutil.copy(str(pdbFile), os.path.join(mdpocketDir, os.path.basename(pdbFile)))
        shutil.copy(str(trajectory), os.path.join(mdpocketDir, os.path.basename(trajectory)))

        return os.path.basename(pdbFile)

    def moveFilePDB(self):
        file = self._getExtraPath("mdpocketInputFile.txt")
        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        shutil.copy(str(file),  os.path.join(mdpocketDir, os.path.basename(file)))

        return os.path.basename(file)

    def cleanUp(self, mdpocketDir):
        outDir = self._getPath()
        pdbDir = self._getExtraPath()

        os.makedirs(outDir, exist_ok=True)
        os.makedirs(pdbDir, exist_ok=True)

        for f in os.listdir(mdpocketDir):
            src = os.path.join(mdpocketDir, f)
            if not os.path.isfile(src):
                continue

            if f.endswith((".txt")):
                shutil.move(src, os.path.join(outDir, f))
            elif f.endswith(".pdb"):
                shutil.move(src, os.path.join(pdbDir, f))

    def runMDPocket(self, pocket):
        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        Plugin.runMDpocket(
            self,
            './mdpocket',
            args=self._getMDpocketSystemArgs(pocket),
            cwd=mdpocketDir
        )
        pdbFile = os.path.basename(self.inputSystem.get().getAttributeValue('pdbFile'))
        trajFile = os.path.basename(self.inputSystem.get().getTrajectoryFile())
        for f in [pdbFile, trajFile]:
            path = os.path.join(mdpocketDir, f)
            if os.path.exists(path):
                os.remove(path)

    def runMDPocketPDB(self, pocket):
        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        Plugin.runMDpocket(
            self,
            './mdpocket',
            args=self._getMDpocketPDBsArgs(pocket),
            cwd=mdpocketDir
        )
        f = "mdpocketInputFile.txt"
        path = os.path.join(mdpocketDir, f)
        if os.path.exists(path):
            os.remove(path)

    def getSpecifiedPocketFile(self):
        myPocket = None
        for pocket in self.pockets.get():
            if pocket.__str__() == self.pocket.get():
                myPocket = pocket.clone()
                break
        if myPocket == None:
            print("Could not find the specified pocket.")
            return None
        else:
            pocketFile = myPocket.getFileName()
            return os.path.abspath(pocketFile)