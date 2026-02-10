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
import zipfile

import os
import re
import shutil
import pyworkflow.protocol.params as params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwchem import Plugin
from biofold.constants import CHAI_DIC

from pwem.objects import  AtomStruct, SetOfAtomStructs


class ProtChai(EMProtocol):
    """
    Protocol to use Chai-1 model.
    """
    _label = 'chai-1 modelling'

    # -------------------------- DEFINE param functions ----------------------
    def _addInputForm(self, form):
        form.addParam('inputSequence', params.PointerParam,
                      pointerClass='Sequence', allowsNull=True,
                      label="Input sequence: ", condition='inputOrigin==0',
                      help='Select the sequence object to add to the set')

        form.addParam('inputAtomStruct', params.PointerParam,
                      pointerClass='AtomStruct', allowsNull=True,
                      label="Input structure: ", condition='inputOrigin==1',
                      help='Select the AtomStruct object whose sequence to add to the set')

        form.addParam('entityType', params.EnumParam, default=0, condition='inputOrigin in [0,1]',
                      label='Input entity type: ', choices=['Protein', 'DNA', 'RNA'],
                      help='Input entity type to add to the set')

        form.addParam('inpChain', params.StringParam,
                      label='Input chain: ', condition='inputOrigin == 1',
                      help='Specify the protein chain to use as sequence.')

        form.addParam('inpPositions', params.StringParam,
                      label='Input positions: ', condition='inputOrigin in [0,1]',
                      help='Specify the positions of the sequence to add in the output.')

        form.addParam('addInput', params.LabelParam,
                      label='Add input: ', condition='inputOrigin in [0,1]',
                      help='Add sequence to the output set')

    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='Input')
        form.addParam('inputOrigin', params.EnumParam, default=0,
                      label='Input origin: ', choices=['Sequence', 'AtomStruct', 'fasta file'],
                      help='Input origin to add to the set')
        self._addInputForm(form)

        form.addParam('inputList', params.TextParam, width=100, condition='inputOrigin in [0,1]',
                      default='', label='List of inputs: ',
                      help='The list of input to use for the final output set.')

        form.addParam('file', params.FileParam, condition='inputOrigin == 2',
                      label='Sequence file: ',
                      help='Select the fasta file. The format of the fasta file can be seen in https://github.com/chaidiscovery/chai-lab/tree/main.')


        form = form.addGroup('Parameters')
        form.addParam('msa', params.BooleanParam, default=True,
                      label="Run with MSAs: ",
                      help='Choose whether to run with MSAs for improved performance.')
        form.addParam('trunkRecycles', params.IntParam, default=3, expertLevel=params.LEVEL_ADVANCED,
                        label='Recycling steps: ', help="Number of recycling steps for prediction.")
        form.addParam('timeSteps', params.IntParam, default=200,
                        label='Sampling steps: ', help="Number of sampling steps for prediction.")
        form.addParam('trunkSamples', params.IntParam, default=1, expertLevel=params.LEVEL_ADVANCED,
                        label='Trunk samples: ', help="Number of trunk samples for prediction.")
        form.addParam('diffNsamples', params.IntParam, default=5, expertLevel=params.LEVEL_ADVANCED,
                       label='Difussion samples for affinity: ', help="Number of diffusion samples for affinity.")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        if (self.inputOrigin.get() != 2):
            self._insertFunctionStep(self.createInputFileStep)
        self._insertFunctionStep(self.runChaiStep)
        self._insertFunctionStep(self.extractScoreStep)
        self._insertFunctionStep(self.createOutputStep)

    def createInputFileStep(self):
        path = os.path.abspath(self._getPath())
        if not os.path.exists(path):
            os.makedirs(path)

        fastaPath = os.path.join(path, "input.fasta")
        with open(fastaPath, 'w') as f:
            counter = 1
            for line in self.inputList.get().splitlines():
                line = line.strip()
                if not line:
                    continue

                try:
                    if ')' in line:
                        json_part = line.split(')', 1)[1].strip()
                    else:
                        json_part = line
                    inpDict = json.loads(json_part)
                except json.JSONDecodeError:
                    continue

                base_name = inpDict.get("name", "unknown")
                seqFile = inpDict.get("seqFile")
                entity = inpDict.get("entity", "protein").lower()

                sequence = ""
                with open(seqFile) as sf:
                    for l in sf:
                        if l.startswith(">"):
                            continue
                        sequence += l.strip()

                unique_name = f"{base_name}_{counter}"
                counter += 1
                f.write(f">{entity}|name={unique_name}\n")
                f.write(f"{sequence}\n")

    def runChaiStep(self):
        if (self.inputOrigin.get() == 2):
            filePath = os.path.abspath(self.file.get())
        else:
            filePath = os.path.abspath(self._getPath('input.fasta'))
        args = [str(filePath)]

        args.append(os.path.join(os.path.abspath(self._getPath()), "chai_results"))

        if self.msa.get():
            args.append("--use-msa-server")

        args.append(f" --num-trunk-recycles {self.trunkRecycles.get()}")
        args.append(f" --num-diffn-timesteps {self.timeSteps.get()}")
        args.append(f" --num-trunk-samples {self.trunkSamples.get()}")
        args.append(f" --num-diffn-samples {self.diffNsamples.get()}")

        Plugin.runCondaCommand(
            self,
            args=" ".join(args),
            condaDic=CHAI_DIC,
            program="chai-lab fold",
            cwd=os.path.abspath(Plugin.getVar(CHAI_DIC['home']))
        )

    def extractScoreStep(self):
        """Extract per-residue score and compute mean score per model"""
        resultsPath = os.path.join(os.path.abspath(self._getPath()), "chai_results")
        self.meanScore = {}

        extraFiles = self.getExtraFiles()

        for cifName in extraFiles:
            cifPath = os.path.join(resultsPath, cifName)
            modelName = os.path.splitext(cifName)[0]

            headers = []
            scoreValues = []
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
                    score = float(cols[colIndex['B_iso_or_equiv']])
                    if resnum not in seenResidues:
                        seenResidues.add(resnum)
                        scoreValues.append(score)

            self.meanScore[modelName] = sum(scoreValues) / len(scoreValues)

        self.bestModel = max(self.meanScore, key=self.meanScore.get)

    def createOutputStep(self):
        resultsPath = os.path.join((self._getPath()), "chai_results")
        extraFiles = self.getExtraFiles()
        outputSet = SetOfAtomStructs.create(self._getPath())
        for cifName in extraFiles:
            atomStruct = AtomStruct(filename=os.path.join(resultsPath, cifName))
            outputSet.append(atomStruct)

        bestSrc = os.path.join(resultsPath, self.bestModel + '.cif')
        bestStruct = AtomStruct(filename=bestSrc)

        self._defineOutputs(
            outputBestAtomStruct=bestStruct,
            outputSetOfAtomStructs=outputSet
        )

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []

        resultsPath = os.path.join(os.path.abspath(self._getPath()), "chai_results")

        try:
            extraFiles = self.getExtraFiles()
        except Exception as e:
            summary.append(f"No CIF files found in {resultsPath}: {e}")
            return summary

        if not extraFiles:
            summary.append(f"No CIF files found in {resultsPath}.")
            return summary

        scores = {}
        logFile = os.path.join(self._getPath('logs'), "run.stdout")
        if os.path.exists(logFile):
            with open(logFile) as f:
                for line in f:
                    m = re.search(r"Score=([\d.]+), writing output to (.+\.cif)", line)
                    if m:
                        score = float(m.group(1))
                        cifPath = m.group(2)
                        modelName = os.path.splitext(os.path.basename(cifPath))[0]
                        scores[modelName] = score

        summary.append("Scores per model:")
        for modelName in sorted(scores.keys()):
            summary.append(f"  {modelName}: {scores[modelName]}")

        if scores:
            bestModel = max(scores, key=scores.get)
            summary.append(f"\nBest structure (highest score): {bestModel}.cif")

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
    def getExtraFiles(self):
        extraFiles = []
        resultsPath = os.path.join(os.path.abspath(self._getPath()), "chai_results")

        for name in sorted(os.listdir(resultsPath)):
            if name.lower().endswith('.cif'):
                extraFiles.append(name)

        if not extraFiles:
            raise Exception("No CIF files found in the selected folder.")

        return extraFiles