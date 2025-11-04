# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb
from autodock.protocols import ProtChemADTPrepareReceptor
from gromacs.protocols import GromacsSystemPrep, GromacsMDSimulation
from ..protocols import MDpocketAnalyze
import os

workflow = '''{'simTime': 100.0, 'timeStep': 0.002, 'nStepsMin': 100, 'emStep': 0.002, 'emTol': 1000.0, 'timeNeigh': 10, 'saveTrj': False, 'trajInterval': 1.0, 'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1, 'restraints': 'Protein', 'restraintForce': 50.0, 'integrator': 'steep', 'ensemType': 'Energy min', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman'}
{'simTime': 0.1, 'timeStep': 0.002, 'nStepsMin': 100, 'emStep': 0.002, 'emTol': 1000.0, 'timeNeigh': 10, 'saveTrj': True, 'trajInterval': 0.05, 'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1, 'restraints': 'MainChain', 'restraintForce': 50.0, 'integrator': 'steep', 'ensemType': 'NVT', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman'}
{'simTime': 0.2, 'timeStep': 0.002, 'nStepsMin': 100, 'emStep': 0.002, 'emTol': 1000.0, 'timeNeigh': 10, 'saveTrj': True, 'trajInterval': 0.05, 'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1, 'restraints': 'None', 'restraintForce': 50.0, 'integrator': 'steep', 'ensemType': 'NPT', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman'}
'''
summary = '''1) Minimization (steep): 100 steps, 1000.0 objective force, restraint on Protein, 300.0 K
2) MD simulation: 0.1 ps, NVT ensemble, restraint on MainChain, 300.0 K
3) MD simulation: 0.2 ps, NPT ensemble, 300.0 K'''


class TestMDPocket(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)
        cls._runImportPDB()
        #cls._runPrepareReceptorADT()
        cls._runPrepareSystem()
        cls._runSimulation()

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/1ake_mut1.pdb'))
        cls.proj.launchProtocol(protImportPDB, wait=True)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runPrepareSystem(cls):
        protPrepare = cls.newProtocol(
            GromacsSystemPrep,
            inputStructure=cls.protImportPDB.outputPdb,
            boxType=1, sizeType=1, padDist=2.0,
            mainForceField=0, waterForceField=2,
            placeIons=1, cationType=7, anionType=1)

        cls.launchProtocol(protPrepare)
        cls.protPrepare = protPrepare

    @classmethod
    def _runPrepareReceptorADT(cls):
        cls.protPrepareReceptor = cls.newProtocol(
            ProtChemADTPrepareReceptor,
            inputAtomStruct=cls.protImportPDB.outputPdb,
            prepProg=1, nphs=False,
            HETATM=True, rchains=False)

        cls.protPrepareReceptor = cls.launchProtocol(cls.protPrepareReceptor)

    @classmethod
    def _runSimulation(self):
        protSim = self.newProtocol(
            GromacsMDSimulation,
            gromacsSystem=self.protPrepare.outputSystem, workFlowSteps=workflow, summarySteps=summary)
        protSim.setObjLabel('gromacs - gmx MD sim')

        outIndex = protSim.getCustomIndexFile()
        if os.path.exists(outIndex):
            protSim.parseIndexFile(outIndex)
        else:
            protSim.createIndexFile(self.protPrepare.outputSystem, inIndex=None, outIndex=protSim.getCustomIndexFile())

        self.launchProtocol(protSim)
        self.protSim = protSim

    def _runMDPocketFind(self):
        protMDPocket = self.newProtocol(
            MDpocketAnalyze,
            useSystem=True,
            inputSystem=self.protSim.outputSystem,
            characterize=False,
            transDruggable=False
        )

        self.launchProtocol(protMDPocket)
        pdbOut = getattr(protMDPocket, 'outputSet', None)
        self.assertIsNotNone(pdbOut)

    def testFpocket(self):
        self._runMDPocketFind()


class TestMDPocketAdvanced(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)
        cls._runImportPDB()
        #cls._runPrepareReceptorADT()
        cls._runPrepareSystem()
        cls._runSimulation()

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/1ake_mut1.pdb'))
        cls.proj.launchProtocol(protImportPDB, wait=True)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runPrepareSystem(cls):
        protPrepare = cls.newProtocol(
            GromacsSystemPrep,
            inputStructure=cls.protImportPDB.outputPdb,
            boxType=1, sizeType=1, padDist=2.0,
            mainForceField=0, waterForceField=2,
            placeIons=1, cationType=7, anionType=1)

        cls.launchProtocol(protPrepare)
        cls.protPrepare = protPrepare

    @classmethod
    def _runPrepareReceptorADT(cls):
        cls.protPrepareReceptor = cls.newProtocol(
            ProtChemADTPrepareReceptor,
            inputAtomStruct=cls.protImportPDB.outputPdb,
            prepProg=1, nphs=False,
            HETATM=True, rchains=False)

        cls.protPrepareReceptor = cls.launchProtocol(cls.protPrepareReceptor)

    @classmethod
    def _runSimulation(self):
        protSim = self.newProtocol(
            GromacsMDSimulation,
            gromacsSystem=self.protPrepare.outputSystem, workFlowSteps=workflow, summarySteps=summary)
        protSim.setObjLabel('gromacs - gmx MD sim')

        outIndex = protSim.getCustomIndexFile()
        if os.path.exists(outIndex):
            protSim.parseIndexFile(outIndex)
        else:
            protSim.createIndexFile(self.protPrepare.outputSystem, inIndex=None, outIndex=protSim.getCustomIndexFile())

        self.launchProtocol(protSim)
        self.protSim = protSim

    def _runMDPocketFindAdvanced(self):
        protMDPocketAdv = self.newProtocol(
            MDpocketAnalyze,
            useSystem=True,
            inputSystem=self.protSim.outputSystem,
            characterize=False,
            transDruggable=False,
            choosePocket=True,
            pockType=1,
            densIsoValue=7.0,
            freqIsoValue=0.4
        )

        self.launchProtocol(protMDPocketAdv)
        pdbOut = getattr(protMDPocketAdv, 'outputSet', None)
        self.assertIsNotNone(pdbOut)


    def testFpocket(self):
        self._runMDPocketFindAdvanced()




