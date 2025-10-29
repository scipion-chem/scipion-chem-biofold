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
This protocol is used to perform a pocket search on a protein structure using the FPocket software

"""

import os, shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pyworkflow.object import String
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfStructROIs, PredictStructROIsOutput, StructROI
from gromacs import Plugin
from gromacs.objects import GromacsSystem
from pwchem.utils import *
from fpocket import Plugin
from fpocket.constants import *
from pwem.convert.atom_struct import cifToPdb

class MDpocketAnalyze(EMProtocol):
    """
    Executes the mdpocket software to look for protein pockets.
    """
    _label = 'Pocket detection'
    _pocketTypes = ['Small molecule binding sites', 'Putative channels and small cavities', 'Water binding sites', 'Big external pockets']
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSystem', params.BooleanParam, deafult=True,
                      label='Use MD system as input: ',
                      help='Select input files, Yes = MD System, No = set of ROIs')
        form.addParam('inputSystem', params.PointerParam, condition='inputSystem',
                       pointerClass='GromacsSystem', allowsNull=True,
                       label="Input gromacs system: ",
                       help='Select the MD system to search for pockets.')
        form.addParam('setOfROIs', params.PointerParam, condition= not 'characterize',
                      pointerClass='SetOfStructROIs', allowsNull=True,
                      label="Input ROIs: ",
                      help='Select the set of ROIs to use as input.')
        form.addParam('characterize', params.BooleanParam, deafult=False,
                      label='Perform pocket characterization: ',
                      help='Perform pocket characterization over a especific set of pockets.')



        group = form.addGroup('Search parameters')
        group.addParam('transDruggable', params.BooleanParam, deafult=False,
                      label='Search transient druggable binding pockets: ',
                      help='Assess at what point the identified pocket is likely to bind drug like molecules.')
        group.addParam('choosePocket', params.BooleanParam, default=False,
                      label='Choose pocket type for advanced search: ',
                      help='Select type of pocket.')

        group.addParam('densIsoValue', params.FloatParam, default=8.0, expertLevel=params.LEVEL_ADVANCED,
                       label='Selected density isovalue: ',
                       help='Creates a pdb file of all grid point positions corresponding to grid points having (isovalue) or more Voronoi vertices nearby per snapshot. A default file will also be created with default isoValue 8.')
        group.addParam('freqIsoValue', params.FloatParam, default=0.5, expertLevel=params.LEVEL_ADVANCED,
                       label='Selected frequency isovalue: ',
                       help='Creates a pdb file of all grid point positions corresponding to grid points that are above selected fraction of the trajectory overlapping with a pocket. A default file will also be created with default isoValue 0.5 (50%).')
        #todo do i need this?
        #group.addParam('maxIntraDistance', params.FloatParam, default='2.0',  condition='choosePocket' and not 'transDruggable',
        #               label='Maximum distance between pocket points (A): ',
        #               help='Maximum distance between two pocket atoms to considered them same pocket')
        group.addParam('pockType', params.EnumParam,
                   choices=self._pocketTypes, default=0,  condition='choosePocket',
                   label='Pocket type:',
                   help='Detect different type of pockets with a set of specific inner parameters.'
                   )

        form.addParallelSection(threads=4)

    def _getMDpocketArgs(self):
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

            if selPock == 'Small molecule binding sites':
                pass
            elif selPock == 'Putative channels and small cavities':
                args.append(' -m 2.8 -M 5.5 -i 3')
            elif selPock == 'Water binding sites':
                args.append('-m 3.5 -M 5.5 -i 3')
            elif selPock == 'Big external pockets':
                args.append('-m 3.5 -M 10.0 -i 3')
        return args

    def _getselIsovalueDensArgs(self):
    # python extractISOPdb.py path/my_dx_file.dx outputname.pdb isovalue

        volFile = os.path.abspath(self._getPath("mdpout_dens_grid.dx")) #To get the extra path between home path and protocol
        args = [volFile]

        outputName = 'mdpoutput_dens_iso_{}.pdb'.format(str(self.isoValue.get()))
        args += [outputName]

        isoValue = self.isoValue.get()
        args += [isoValue]
        return args

    def _getselIsovalueFreqArgs(self):
    # python extractISOPdb.py path/my_dx_file.dx outputname.pdb isovalue

        volFile = os.path.abspath(self._getPath("mdpout_freq_grid.dx")) #To get the extra path between home path and protocol
        args = [volFile]

        outputName = 'mdpoutput_freq_iso_{}.pdb'.format(str(self.isoValue.get()).replace('.', '_'))
        args += [outputName]

        isoValue = self.isoValue.get()
        args += [isoValue]
        return args


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('mdPocketStep')
        self._insertFunctionStep('selIsovalue')
        #self._insertFunctionStep('defineOutputStep')

    def mdPocketStep(self):
        mdpocketDir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))

        Plugin.runMDpocket(
            self,
            './mdpocket',
            args=self._getMDpocketArgs(),
            cwd=mdpocketDir
        )

        #move resulting files to specific paths
        outDir = self._getPath()
        for f in ["mdpout_dens_grid.dx",
                  "mdpout_freq_grid.dx",
                  "time.txt"]:
            src = os.path.join(mdpocketDir, f)
            if os.path.exists(src):
                shutil.move(src, os.path.join(outDir, f))
        pdbDir = self._getExtraPath()
        for f in ["mdpout_all_atom_pdensities.pdb",
                  "mdpout_dens_iso_8.pdb",
                  "mdpout_freq_iso_0_5.pdb"]:
            src = os.path.join(mdpocketDir, f)
            if os.path.exists(src):
                shutil.move(src, os.path.join(pdbDir, f))


        pdbFile = os.path.basename(self.inputSystem.get().getAttributeValue('pdbFile'))
        trajFile = os.path.basename(self.inputSystem.get().getTrajectoryFile())

        for f in [pdbFile, trajFile]:
            path = os.path.join(mdpocketDir, f)
            if os.path.exists(path):
                os.remove(path)



    def selIsovalue(self):
        if (self.densIsoValue.get() != 8.0):
            Plugin.runSelIsovalue(self, 'extractISOPdb.py', args=self._getselIsovalueDensArgs(), cwd=self._getExtraPath())
        if (self.freqIsoValue.get() != 0.5):
            Plugin.runSelIsovalue(self, 'extractISOPdb.py', args=self._getselIsovalueFreqArgs(), cwd=self._getExtraPath())

    def defineOutputStep(self):
        #todo change this, build StructROIs with pdb files and create SetOfStructROIs as output
        pocketsDir = self._getExtraPath()
        pocketFiles = os.listdir(pocketsDir)

        outPockets = SetOfStructROIs(filename=self._getExtraPath('pockets.sqlite'))
        for pFile in pocketFiles:
            if '.pdb' in pFile:
                #todo find nPoints, contactResidues, contactAtoms
                roi = StructROI(filename=os.path.join(pocketsDir, pFile), pClass='MDpocket')


        outPockets.buildPDBhetatmFile()
        self._defineOutputs()



    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = ["If different isovalue were selected, the new pdb files appear in the 'extra' folder."]
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
        pdbFile = self.inputSystem.get().getAttributeValue('pdbFile')
        # move files to path where mdpocket is, it is picky with where it is executed and they input files routes
        mdpocket_dir = os.path.abspath(os.path.join(Plugin.getVar(FPOCKET_DIC['home']), 'bin'))
        shutil.copy(str(pdbFile), os.path.join(mdpocket_dir, os.path.basename(pdbFile)))
        shutil.copy(str(trajectory), os.path.join(mdpocket_dir, os.path.basename(trajectory)))

        return os.path.basename(pdbFile)