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

from biofold.protocols import (ProtChai, ProtBoltz, ProtIntelliFold, ProtProtenix,
                               ProtBoltzCofolding, ProtIntelliFoldCofolding,
                               ProtProtenixCofolding, ProtChaiCofolding)
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


class TestCofoldingBase(BaseTest):
    """Shared setup and checks for the cofolding protocols.

    Every backend consumes the same receptor and ligands and must produce a docked
    SetOfSmallMolecules.
    """
    #: set by each subclass
    protClass, protKwargs, extraAttrs = None, {}, ()

    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet('smallMolecules')
        setupTestProject(cls)

    def _runInputs(self):
        protPdb = self.newProtocol(ProtImportPdb, inputPdbData=1,
                                   pdbFile=self.ds.getFile('PDBx_mmCIF/1ake_start.pdb'))
        self.launchProtocol(protPdb)
        protMols = self.newProtocol(ProtChemImportSmallMolecules,
                                    filesPath=self.dsLig.getFile('sdf'), filesPattern='200[01].sdf')
        self.launchProtocol(protMols)
        return protPdb, protMols

    def _checkCofolding(self):
        """Run the backend on 2 ligands and check the docked output it must produce"""
        protPdb, protMols = self._runInputs()
        prot = self.newProtocol(self.protClass,
                                inputAtomStruct=protPdb.outputPdb,
                                inputSmallMolecules=protMols.outputSmallMolecules,
                                **self.protKwargs)
        self.launchProtocol(prot)

        out = getattr(prot, 'outputSmallMolecules', None)
        self.assertIsNotNone(out, 'no outputSmallMolecules was produced')
        self.assertTrue(out.isDocked())
        # 2 ligands x 2 poses each
        self.assertEqual(out.getSize(), 4)
        self.assertTrue(os.path.exists(out.getProteinFile()))

        for mol in out:
            self.assertTrue(os.path.exists(mol.getPoseFile()))
            # each pose carries its own predicted receptor, for the viewer and for MD prep
            self.assertTrue(os.path.exists(mol.getProteinFile()))
            for attr in ('_cofoldLigPLDDT', '_cofoldLigIpTM', '_cofoldConfidence',
                         '_cofoldRecRMSD') + self.extraAttrs:
                self.assertIsNotNone(getattr(mol, attr, None), f'{attr} missing')
            # the ligand pLDDT must be the ligand's own, not the complex-wide one
            self.assertTrue(0 <= mol._cofoldLigPLDDT.get() <= 1)


class TestBoltzCofolding(TestCofoldingBase):
    protClass = ProtBoltzCofolding
    # low settings, and a capped MSA, so it fits small GPUs
    protKwargs = dict(diffusionSamples=2, recyclingSteps=1, samplingSteps=50,
                      maxParallelSamples=1, maxMsaSeqs=512, diffusionSamplesAff=2,
                      samplingStepsAff=50, infPot=False)
    # Boltz-2 is the only backend that predicts affinity
    extraAttrs = ('_boltzAffinity', '_boltzIntProbability')

    def test(self):
        self._checkCofolding()


class TestIntelliFoldCofolding(TestCofoldingBase):
    protClass = ProtIntelliFoldCofolding
    protKwargs = dict(diffusionSamples=2, recyclingSteps=1, samplingSteps=50)

    def test(self):
        self._checkCofolding()


class TestProtenixCofolding(TestCofoldingBase):
    protClass = ProtProtenixCofolding
    protKwargs = dict(diffusionSamples=2, recyclingSteps=1, samplingSteps=50)

    def test(self):
        self._checkCofolding()


class TestChaiCofolding(TestCofoldingBase):
    protClass = ProtChaiCofolding
    protKwargs = dict(diffusionSamples=2, trunkRecycles=1, samplingSteps=50)

    def test(self):
        self._checkCofolding()
