# **************************************************************************
# *
# * Authors:   Blanca Pueche (blanca.pueche@cnb.csis.es)
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
import subprocess
import unittest

import os

from biofold.protocols import ProtChai, ProtBoltz, ProtIntelliFold, ProtProtenix, ProtBoltzCofolding
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb
from pwchem.protocols import ProtChemImportSmallMolecules


defSetASChain, defSetPDBChain = 'A', 'B'
defSetPDBFile = 'Tmp/5ni1_{}_FIRST-LAST.fa'.format(defSetPDBChain)

names = ['5ni1']
defSetChains = [None, defSetASChain, defSetPDBChain]
defSetFiles = [defSetPDBFile]

defSetSeqs = '''1) {"name": "%s", "chain": "%s", "index": "FIRST-LAST", "seqFile": "%s"}''' % \
                         (names[0], defSetPDBChain, defSetPDBFile)


class TestChai(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)

    def _runChai(self):
        protChai = self.newProtocol(
            ProtChai,
            inputOrigin=3,
            file=self.ds.getFile('Sequences/3lqd_B_mutated.fasta')
        )

        self.launchProtocol(protChai)
        best = getattr(protChai, 'outputBestAtomStruct', None)
        self.assertIsNotNone(best)
        all = getattr(protChai, 'outputSetOfAtomStructs', None)
        self.assertIsNotNone(all)

    def test(self):
        self._runChai()

class TestBoltz(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)

    def _runBoltz(self):
        protBoltz = self.newProtocol(
            ProtBoltz,
            inputOrigin=3,
            entityType=1,
            recyclingSteps=1,
            samplingSteps=50,
            file=self.ds.getFile('Sequences/3lqd_B_mutated.fasta')
        )

        self.launchProtocol(protBoltz)
        best = getattr(protBoltz, 'outputBestAtomStruct', None)
        self.assertIsNotNone(best)
        all = getattr(protBoltz, 'outputSetOfAtomStructs', None)
        self.assertIsNotNone(all)

    def test(self):
        self._runBoltz()

class TestIntelliFold(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)

    def _runIntelliFold(self):
        protIntelliFold = self.newProtocol(
            ProtIntelliFold,
            inputOrigin=3,
            entityType=1,
            recyclingSteps=1,
            samplingSteps=50,
            file=self.ds.getFile('Sequences/3lqd_B_mutated.fasta')
        )

        self.launchProtocol(protIntelliFold)
        best = getattr(protIntelliFold, 'outputBestAtomStruct', None)
        self.assertIsNotNone(best)
        all = getattr(protIntelliFold, 'outputSetOfAtomStructs', None)
        self.assertIsNotNone(all)

    def test(self):
        self._runIntelliFold()

class TestProtenix(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)

    def _runProtenix(self):
        protProtenix = self.newProtocol(
            ProtProtenix,
            inputOrigin=3,
            entityType=1,
            recyclingSteps=1,
            samplingSteps=50,
            file=self.ds.getFile('Sequences/3lqd_B_mutated.fasta'),
            model=1
        )

        self.launchProtocol(protProtenix)
        best = getattr(protProtenix, 'outputBestAtomStruct', None)
        self.assertIsNotNone(best)
        all = getattr(protProtenix, 'outputSetOfAtomStructs', None)
        self.assertIsNotNone(all)

    def test(self):
        self._runProtenix()


class TestBoltzCofolding(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet('smallMolecules')
        setupTestProject(cls)

    def test(self):
        protPdb = self.newProtocol(ProtImportPdb, inputPdbData=1,
                                   pdbFile=self.ds.getFile('PDBx_mmCIF/1ake_start.pdb'))
        self.launchProtocol(protPdb)
        protMols = self.newProtocol(ProtChemImportSmallMolecules,
                                    filesPath=self.dsLig.getFile('sdf'), filesPattern='200[01].sdf')
        self.launchProtocol(protMols)

        # Low settings (and a capped MSA) so it fits small GPUs
        protCofold = self.newProtocol(ProtBoltzCofolding,
                                      inputAtomStruct=protPdb.outputPdb,
                                      inputSmallMolecules=protMols.outputSmallMolecules,
                                      diffusionSamples=2, recyclingSteps=1, samplingSteps=50,
                                      maxParallelSamples=1, maxMsaSeqs=512, diffusionSamplesAff=2)
        self.launchProtocol(protCofold)

        out = getattr(protCofold, 'outputSmallMolecules', None)
        self.assertIsNotNone(out)
        self.assertTrue(out.isDocked())
        self.assertEqual(out.getSize(), 4)
        for mol in out:
            self.assertTrue(os.path.exists(mol.getPoseFile()))
            self.assertTrue(os.path.exists(mol.getProteinFile()))
            self.assertIsNotNone(mol._boltzAffinity.get())
            self.assertIsNotNone(mol._boltzIntProbability.get())
