# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
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
This protocol is used to perform a pocket search on a protein structure using the FPocket software

"""
import math
import os.path

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import *
from fpocket import Plugin
from fpocket.constants import *

class MDpocketAnalyze(EMProtocol):
    """
    Executes the mdpocket software to look for protein pockets.
    """
    _label = 'MDPocket pocket detection'
    _pocketTypes = ['Small molecule binding sites', 'Putative channels and small cavities', 'Water binding sites', 'Big external pockets']
    stepsExecutionMode = params.STEPS_PARALLEL
    _customDens = 'mdpout_dens_iso_custom.pdb (User-selected density isovalue)'
    _customFreq = 'mdpout_freq_iso_custom.pdb (User-selected frequency isovalue)'
    _inputFileTxt = 'mdpocketInputFile.txt'
    _inputSystemPDB = 'inputSystem.pdb'
    _inputSystemPDBOpenMM = ''
    _help='Choose which of the detected pockets to characterize.'
    _labelParam = 'Pocket to characterize:'
    _mdpocketprogram = './mdpocket'

    choices1Both = [
        'mdpout_dens_iso_8.pdb (Default density grid, isovalue 8.0)',
        'mdpout_freq_iso_0_5.pdb (Default frequency grid, isovalue 0.5)',
        _customDens,
        _customFreq
    ]
    choices2Both = [
        _customDens,
        _customFreq
    ]
    choices1Dens = [
        'mdpout_dens_iso_8.pdb (Default density grid, isovalue 8.0)',
        _customDens
    ]
    choices2Dens = [
        _customDens
    ]
    choices1Freq = [
        'mdpout_freq_iso_0_5.pdb (Default frequency grid, isovalue 0.5)',
        _customFreq
    ]
    choices2Freq = [
        _customFreq
    ]

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('useSystem', params.BooleanParam, deafult=True,
                      label='Use MD system as input: ',
                      help='Select input files, Yes = MD System, No = set of pdbs')
        form.addParam('inputSystem', params.PointerParam, condition='useSystem',
                       pointerClass='GromacsSystem, OpenMMSystem', allowsNull=True,
                       label="Input system: ",
                       help='Select the MD system to search for pockets.')
        form.addParam('inputPDBs', params.PointerParam, condition='not useSystem',
                      pointerClass='SetOfStructROIs,SetOfAtomStructs', allowsNull=True,
                      label="Input set of struct ROIs/set of pdbs: ",
                      help='Select the structural ROIs to use as input.')

        group = form.addGroup('Search parameters')
        group.addParam('transDruggable', params.BooleanParam, deafult=False,
                      label='Search transient druggable binding pockets: ',
                      help='Assess at what point the identified pocket is likely to bind drug like molecules.')
        group.addParam('choosePocket', params.BooleanParam, default=False,
                      label='Choose pocket type for advanced search: ',
                      help='Select type of pocket.')

        group = form.addGroup('Output generation')
        group.addParam('chooseOutput', params.EnumParam, deafult=2, choices=['Frequency','Density','Both'],
                       label='Output files: ',
                       help='Choose what outputs to keep.')
        group.addParam('densIsoValue', params.FloatParam, default=8.0, expertLevel=params.LEVEL_ADVANCED, condition='chooseOutput==1 or chooseOutput==2',
                       label='Selected density isovalue: ',
                       help='Creates a pdb file of all grid point positions corresponding to grid points having (isovalue) or more Voronoi vertices nearby per snapshot. A default file will also be created with default isoValue 8.')
        group.addParam('freqIsoValue', params.FloatParam, default=0.5, expertLevel=params.LEVEL_ADVANCED, condition='chooseOutput==0 or chooseOutput==2',
                       label='Selected frequency isovalue: ',
                       help='Creates a pdb file of all grid point positions corresponding to grid points that are above selected fraction of the trajectory overlapping with a pocket. A default file will also be created with default isoValue 0.5 (50%).')
        group.addParam('keepDefaultFiles', params.BooleanParam, default=True,
                       label='Keep result files with default parameters: ',
                       help='The default files are created always, choose whether to keep them or not.')
        group.addParam('pockType', params.EnumParam,
                   choices=self._pocketTypes, default=0,  condition='choosePocket',
                   label='Pocket type:',
                   help='Detect different type of pockets with a set of specific inner parameters.'
                   )

        group = form.addGroup('Pocket characterization')
        group.addParam('pockTypeDefBoth', params.EnumParam, condition='keepDefaultFiles and chooseOutput==2',
                       choices=self.choices1Both, default=0,
                       label=self._labelParam,
                       help=self._help
                       )
        group.addParam('pockTypeNotDefBoth', params.EnumParam, condition='not keepDefaultFiles and chooseOutput==2',
                       choices=self.choices2Both, default=0,
                       label=self._labelParam,
                       help=self._help
                       )
        group.addParam('pockTypeDefDens', params.EnumParam, condition='keepDefaultFiles and chooseOutput==1',
                       choices=self.choices1Dens, default=0,
                       label=self._labelParam,
                       help=self._help
                       )
        group.addParam('pockTypeNotDefDens', params.EnumParam, condition='not keepDefaultFiles and chooseOutput==1',
                       choices=self.choices2Dens, default=0,
                       label=self._labelParam,
                       help=self._help
                       )
        group.addParam('pockTypeDefFreq', params.EnumParam, condition='keepDefaultFiles and chooseOutput==0',
                       choices=self.choices1Freq, default=0,
                       label=self._labelParam,
                       help=self._help
                       )
        group.addParam('pockTypeNotDefFreq', params.EnumParam, condition='not keepDefaultFiles and chooseOutput==0',
                       choices=self.choices2Freq, default=0,
                       label=self._labelParam,
                       help=self._help
                       )
        group.addParam('distanceClustering', params.FloatParam, default=5.0,
                       label="Threshold for clustering: ",
                       help='Select the distance threshold to create individual pockets from mdpocket output.')

        form.addParallelSection(threads=4)

    def _getMDpocketDefSystemArgs(self):
        pdbFile = self.moveFiles()

        trajFile = self.inputSystem.get().getTrajectoryFile()
        trajBasename = os.path.basename((trajFile))
        args = ['--trajectory_file', trajBasename]

        trajExt = os.path.splitext(trajFile)[1][1:]
        args.append('--trajectory_format')
        args.append(trajExt)

        args.append('-f')
        args.append(pdbFile)

        if (self.transDruggable.get()):
            args.append('-S')

        if (self.choosePocket.get()):
            selPock = self.getEnumText('pockType')

            if selPock == 'Putative channels and small cavities':
                args.append(' -m 2.8 -M 5.5 -i 3')
            elif selPock == 'Water binding sites':
                args.append('-m 3.5 -M 5.5 -i 3')
            elif selPock == 'Big external pockets':
                args.append('-m 3.5 -M 10.0 -i 3')
        return args

    def _getMDpocketPDBsArgs(self):
        routes = []
        for pocket in self.inputPDBs.get():
            routes.append(os.path.abspath(str(pocket.getFileName())))

        inputFile = self._getExtraPath(self._inputFileTxt)
        with open(inputFile, 'w') as f:
            for pdbRoute in routes:
                f.write(pdbRoute + '\n')

        movedFile = self.moveFilePDB()
        pluginArgs = ["-L", movedFile]

        return pluginArgs


    def _getselIsovalueDensArgs(self):
        # python extractISOPdb.py path/my_dx_file.dx outputname.pdb isovalue
        volFile = os.path.abspath(self._getPath("mdpout_dens_grid.dx")) #To get the extra path between home path and protocol
        args = [volFile]

        outputName = 'mdpoutput_dens_iso_{}.pdb'.format(str(self.densIsoValue.get()).replace('.', '_'))
        args += [outputName]

        isoValue = self.densIsoValue.get()
        args += [isoValue]
        return args

    def _getselIsovalueFreqArgs(self):
        # python extractISOPdb.py path/my_dx_file.dx outputname.pdb isovalue
        volFile = os.path.abspath(self._getPath("mdpout_freq_grid.dx")) #To get the extra path between home path and protocol
        args = [volFile]

        outputName = 'mdpoutput_freq_iso_{}.pdb'.format(str(self.freqIsoValue.get()).replace('.', '_'))
        args += [outputName]

        isoValue = self.freqIsoValue.get()
        args += [isoValue]
        return args


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('mdPocketDefStep')
        self._insertFunctionStep('selIsovalue')
        self._insertFunctionStep('pocketDefOutputStep')
        self._insertFunctionStep('mdPocketCharactStep')
        self._insertFunctionStep('defineOutputStep')

    def mdPocketDefStep(self):
        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        if self.useSystem.get():
            #use system
            Plugin.runMDpocket(
                self,
                self._mdpocketprogram,
                args=self._getMDpocketDefSystemArgs(),
                cwd=mdpocketDir
            )
            if 'GromacsSystem' in type(self.inputSystem.get()).__name__:
                pdbFile = self._inputSystemPDB
            else:
                pdbFile = self._inputSystemPDBOpenMM
            trajFile = os.path.basename(self.inputSystem.get().getTrajectoryFile())
            for f in [pdbFile, trajFile]:
                path = os.path.join(mdpocketDir, f)
                if os.path.exists(path):
                    os.remove(path)
        else:
            #use pdb files
            Plugin.runMDpocket(
                self,
                self._mdpocketprogram,
                args=self._getMDpocketPDBsArgs(),
                cwd=mdpocketDir
            )
            f = self._inputFileTxt
            path = os.path.join(mdpocketDir, f)
            if os.path.exists(path):
                os.remove(path)

        self.cleanUp(mdpocketDir)

    def selIsovalue(self):
        scriptDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'scripts'))
        if ((self.chooseOutput.get() == 1 or self.chooseOutput.get()==2) and not math.isclose(self.densIsoValue.get(), 8.0, rel_tol=1e-9, abs_tol=1e-9)):
            Plugin.runScript(self, 'extractISOPdb.py', args=self._getselIsovalueDensArgs(), cwd=scriptDir)
        if ((self.chooseOutput.get() == 0 or self.chooseOutput.get()==2) and not math.isclose(self.freqIsoValue.get(), 0.5, rel_tol=1e-9, abs_tol=1e-9)):
            Plugin.runScript(self, 'extractISOPdb.py', args=self._getselIsovalueFreqArgs(), cwd=scriptDir)

        self.cleanUp(scriptDir)

    def pocketDefOutputStep(self):
        pocketsDir = self._getExtraPath('pocketDetection')
        if not self.keepDefaultFiles.get():
            defaultFiles = [f'{pocketsDir}/mdpout_dens_iso_8.pdb', f'{pocketsDir}/mdpout_freq_iso_0_5.pdb']
            for f in defaultFiles:
                if os.path.exists(f):
                    os.remove(f)

    def mdPocketCharactStep(self):
        if self.useSystem.get():
            self.runMDPocket()
        else:
            self.runMDPocketPDB()

        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        self.cleanUp2(mdpocketDir)

    def defineOutputStep(self):
        pocketsDir = self._getExtraPath('pocketCharacterization')
        pocketFiles = os.listdir(pocketsDir)

        self._cleanupDefaultFile(pocketsDir)
        proteinFile = self._getProteinFile()

        outputRois = SetOfStructROIs(filename=self._getPath('pockets.sqlite'))

        for pFile in pocketFiles:
            if 'mdpocket.pdb' not in pFile:
                continue
            self._processPocketFile(pFile, pocketsDir, proteinFile, outputRois)

        outputRois.buildPDBhetatmFile()
        self._defineOutputs(outputROIs=outputRois)



    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = ["If different isovalues were selected, the new pdb files appear in the 'extra/pocketDetection' folder."]
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
        trajFile = self.inputSystem.get().getTrajectoryFile()
        trajectory = os.path.abspath((trajFile))
        if 'GromacsSystem' in type(self.inputSystem.get()).__name__:
            pdbFile = self.getPath(self._inputSystemPDB)
            self.convertGroToPDB(self.inputSystem.get().getSystemFile(), pdbFile)
        else:
            systemPath = os.path.dirname(trajectory)
            pdbFile = None
            for f in os.listdir(systemPath):
                if f.endswith(".pdb"):
                    pdbFile = os.path.join(systemPath, f)
                    self._inputSystemPDBOpenMM = os.path.basename(pdbFile)
                    break
        # move files to path where mdpocket is, it is picky with where it is executed and they input files routes
        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        shutil.copy(str(pdbFile), os.path.join(self._getExtraPath(), os.path.basename(pdbFile)))
        copyPdbFile = self._getExtraPath(os.path.basename(pdbFile))
        shutil.copy(str(copyPdbFile), os.path.join(mdpocketDir, os.path.basename(copyPdbFile)))
        shutil.copy(str(trajectory), os.path.join(mdpocketDir, os.path.basename(trajectory)))

        return os.path.basename(pdbFile)

    def moveFilePDB(self):
        file = self._getExtraPath(self._inputFileTxt)
        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        shutil.copy(str(file),  os.path.join(mdpocketDir, os.path.basename(file)))

        return os.path.basename(file)

    def cleanUp(self, mdpocketDir):
        outDir = self._getPath()
        pdbDir = self._getExtraPath('pocketDetection')

        os.makedirs(outDir, exist_ok=True)
        os.makedirs(pdbDir, exist_ok=True)

        for f in os.listdir(mdpocketDir):
            src = os.path.join(mdpocketDir, f)
            if not os.path.isfile(src):
                continue

            if f.endswith((".dx", ".txt")):
                shutil.move(src, os.path.join(outDir, f))
            elif f.endswith(".pdb"):
                shutil.move(src, os.path.join(pdbDir, f))

    def cleanUp2(self, mdpocketDir):
        outDir = self._getPath()
        pdbDir = self._getExtraPath('pocketCharacterization')

        os.makedirs(outDir, exist_ok=True)
        os.makedirs(pdbDir, exist_ok=True)

        for f in os.listdir(mdpocketDir):
            src = os.path.join(mdpocketDir, f)
            if not os.path.isfile(src):
                continue

            if f.endswith((".dx", ".txt")):
                shutil.move(src, os.path.join(outDir, f))
            elif f.endswith(".pdb"):
                shutil.move(src, os.path.join(pdbDir, f))

    def runMDPocket(self):
        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        Plugin.runMDpocket(
            self,
            self._mdpocketprogram,
            args=self._getMDpocketCharactSystemArgs(),
            cwd=mdpocketDir
        )
        if 'GromacsSystem' in type(self.inputSystem.get()).__name__:
            pdbFile = self.getPath(self._inputSystemPDB)
        else:
            pdbFile = self._getExtraPath(self._inputSystemPDBOpenMM)

        trajFile = os.path.basename(self.inputSystem.get().getTrajectoryFile())
        for f in [os.path.basename(pdbFile), trajFile]:
            path = os.path.join(mdpocketDir, f)
            if os.path.exists(path):
                os.remove(path)

    def _getMDpocketCharactSystemArgs(self):
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
        specificPocket = self.moveChosenPocket()
        args.append(specificPocket)

        return args

    def moveChosenPocket(self):
        specificPocket = os.path.abspath(self.getSpecifiedPocketFile())
        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        shutil.copy(str(specificPocket), os.path.join(mdpocketDir, os.path.basename(specificPocket)))
        return os.path.basename(specificPocket)

    def getSpecifiedPocketFile(self):
        filesPath = self._getExtraPath('pocketDetection')

        selector, choices = self._getSelectorAndChoices()
        if selector is None or choices is None:
            return ''

        selected = selector.get()
        if selected < 0 or selected >= len(choices):
            return ''

        name = choices[selected].split()[0]
        fileName = self.checkCustomOrDef(name)
        pocketFile = os.path.join(filesPath, fileName)

        return os.path.abspath(pocketFile)

    def _getSelectorAndChoices(self):
        choose = self.chooseOutput.get()
        keep = self.keepDefaultFiles.get()

        if choose == 2:  # BOTH
            if keep:
                return self.pockTypeDefBoth, self.choices1Both
            else:
                return self.pockTypeNotDefBoth, self.choices2Both

        elif choose == 1:  # DENS
            if keep:
                return self.pockTypeDefDens, self.choices1Dens
            else:
                return self.pockTypeNotDefDens, self.choices2Dens

        elif choose == 0:  # FREQ
            if keep:
                return self.pockTypeDefFreq, self.choices1Freq
            else:
                return self.pockTypeNotDefFreq, self.choices2Freq

        return None, None

    def checkCustomOrDef(self, name):
        filename = ''
        if 'dens_iso_custom' in name:
            isoStr = str(self.densIsoValue.get()).replace('.', '_')
            filename = name.replace('custom', isoStr)
        elif 'freq_iso_custom' in name:
            isoStr = str(self.freqIsoValue.get()).replace('.', '_')
            filename = name.replace('custom', isoStr)
        else:
            filename = name
        return filename

    def runMDPocketPDB(self):
        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        Plugin.runMDpocket(
            self,
            self._mdpocketprogram,
            args=self._getMDpocketCharactPDBsArgs(),
            cwd=mdpocketDir
        )
        f = self._inputFileTxt
        path = os.path.join(mdpocketDir, f)
        if os.path.exists(path):
            os.remove(path)

    def _getMDpocketCharactPDBsArgs(self):
        routes = []
        for pocket in self.inputPDBs.get():
            routes.append((str(pocket.getFileName())))

        inputFile = self._getExtraPath(self._inputFileTxt)
        with open(inputFile, 'w') as f:
            for pdbRoute in routes:
                f.write(os.path.abspath(pdbRoute) + '\n')


        movedFile = self.moveFilePDB()
        specificPocket = self.moveChosenPocket()
        pluginArgs = ["-L",(movedFile)]

        pluginArgs.append('--selected_pocket')
        p = (specificPocket)
        pluginArgs.append(p)

        return pluginArgs

    def convertGroToPDB(self, input, output):
        scriptArgs = [os.path.abspath(input), os.path.abspath(output)]
        Plugin.runMyScript(self, "groToPdb.py", args=scriptArgs)


    def runClustering(self, pFile, dir):
        file = os.path.join(self._getExtraPath('pocketCharacterization'), pFile)
        scriptArgs = [os.path.abspath(file), self.distanceClustering.get(),
                       os.path.abspath(dir)]
        Plugin.runMyScript(self, "splitPockets.py", args=scriptArgs)

    def _cleanupDefaultFile(self, pocketsDir):
        if not self.chooseOutput.get():
            defaultFile = f'{pocketsDir}/mdpout_mdpocket_atoms.pdb'
            if os.path.exists(defaultFile):
                os.remove(defaultFile)

    def _getProteinFile(self):
        if self.useSystem.get():
            if 'GromacsSystem' in type(self.inputSystem.get()).__name__:
                return self.getPath(self._inputSystemPDB)
            else:
                trajFile = self.inputSystem.get().getTrajectoryFile()
                trajectory = os.path.abspath((trajFile))
                systemPath = os.path.dirname(trajectory)
                return os.path.join(systemPath, self._inputSystemPDBOpenMM)

        return self.inputPDBs.get().getFirstItem().getFileName()

    def _processPocketFile(self, pFile, pocketsDir, proteinFile, outputRois):
        clustersDir = self._getExtraPath('pocketCharacterization/pocketsFile')

        try:
            self.runClustering(pFile, clustersDir)
            self._addClusteredPockets(clustersDir, pFile, proteinFile, outputRois)

        except Exception as e:
            print(f"[WARNING] Clustering failed for {pFile}: {e}")
            self._addFallbackPocket(pocketsDir, pFile, proteinFile, outputRois)

    def _addClusteredPockets(self, clustersDir, pFile, proteinFile, outputRois):
        for c in os.listdir(clustersDir):
            if 'corrected' in c:
                continue

            outPocket = StructROI(
                filename=os.path.join(clustersDir, c),
                proteinFile=proteinFile,
                extraFile=os.path.abspath(pFile),
                pClass='MDPocket'
            )
            outPocket.setVolume(outPocket.getPocketVolume())
            outputRois.append(outPocket)

    def _addFallbackPocket(self, pocketsDir, pFile, proteinFile, outputRois):
        outPocket = StructROI(
            filename=os.path.join(pocketsDir, pFile),
            proteinFile=proteinFile,
            pClass='MDPocket'
        )
        outPocket.setVolume(outPocket.getPocketVolume())
        outputRois.append(outPocket)


