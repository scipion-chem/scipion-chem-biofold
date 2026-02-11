# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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
import os

from biofold.protocols import ProtBoltz, ProtChai
from pwem.objects import AtomStruct, Sequence
from pwem.wizards import SelectResidueWizard
from pyworkflow.object import Pointer

from pwchem.wizards.wizard_select_chain import SelectChainWizardQT, SelectResidueWizardQT

SelectChainWizardQT().addTarget(protocol=ProtBoltz,
                                targets=['inpChain'],
                                inputs=[{'inputOrigin': ['inputSequence',
                                                         'inputAtomStruct']}],
                                outputs=['inpChain'])

SelectResidueWizardQT().addTarget(protocol=ProtBoltz,
                                  targets=['inpPositions'],
                                  inputs=[{'inputOrigin': ['inputSequence', 'inputAtomStruct']},
                                          'inpChain'],
                                  outputs=['inpPositions'])

SelectChainWizardQT().addTarget(protocol=ProtChai,
                                targets=['inpChain'],
                                inputs=[{'inputOrigin': ['inputSequence',
                                                         'inputAtomStruct']}],
                                outputs=['inpChain'])

SelectResidueWizardQT().addTarget(protocol=ProtChai,
                                  targets=['inpPositions'],
                                  inputs=[{'inputOrigin': ['inputSequence', 'inputAtomStruct']},
                                          'inpChain'],
                                  outputs=['inpPositions'])

class AddSequenceWizardBoltz(SelectResidueWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol
        inputParams, outputParam = self.getInputOutput(form)

        # StructName
        inputObj = getattr(protocol, inputParams[0]).get()
        pdbFile, AS, addPointer = '', False, True
        if issubclass(type(inputObj), str):
            outStr = [inputObj]
            AS, addPointer = True, False
        elif issubclass(type(inputObj), AtomStruct):
            pdbFile = inputObj.getFileName()
            outStr = [os.path.splitext(os.path.basename(pdbFile))[0]]
            AS = True
        elif issubclass(type(inputObj), Sequence):
            outStr = [inputObj.getId().replace("|", "_")]

        if AS:
            # Chain
            chainJson = getattr(protocol, inputParams[1]).get()
            chainId = json.loads(chainJson)['chain']
        else:
            chainJson = ''

        # Positions
        posJson = getattr(protocol, inputParams[2]).get()
        if posJson:
            posIdxs = json.loads(posJson)['index']
            seq = json.loads(posJson)['residues']
            outStr += [posIdxs]
        else:
            outStr += ['FIRST-LAST']
            finalResiduesList = self.getResidues(form, inputObj, chainJson)
            idxs = [json.loads(finalResiduesList[0].get())['index'], json.loads(finalResiduesList[-1].get())['index']]
            seq = self.getSequence(finalResiduesList, idxs)

        chainStr, chainFileId = '', ''
        if AS:
          chainStr = ', "chain": "{}"'.format(chainId)
          chainFileId = '_{}'.format(chainId)

        prevStr = getattr(protocol, outputParam[0]).get()
        lenPrev = len(prevStr.strip().split('\n')) + 1
        if prevStr.strip() == '':
          lenPrev -= 1
        elif not prevStr.endswith('\n'):
          prevStr += '\n'

        cyclic = getattr(protocol, 'cyclic').get()
        val = getattr(protocol, 'entityType').get()
        if val==0: entity = 'protein'
        elif val==1: entity = 'dna'
        else: entity = 'rna'

        seqFile = protocol.getProject().getTmpPath('{}{}_{}.fa'.format(outStr[0], chainFileId, outStr[1]))
        with open(seqFile, 'w') as f:
            f.write('>{}\n{}\n'.format(outStr[0], seq))

        jsonStr = '%s) {"name": "%s"%s, "index": "%s", "seqFile": "%s", "entity": "%s", "cyclic": "%s"}\n' % \
                  (lenPrev, outStr[0], chainStr, outStr[1], seqFile, entity, cyclic)
        form.setVar(outputParam[0], prevStr + jsonStr)


AddSequenceWizardBoltz().addTarget(protocol=ProtBoltz,
                              targets=['addInput'],
                              inputs=[{'inputOrigin': ['inputSequence', 'inputAtomStruct', 'inputPDB']},
                                      'inpChain', 'inpPositions'],
                              outputs=['inputList', 'inputPointers'])

class AddSequenceWizardChai(SelectResidueWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol
        inputParams, outputParam = self.getInputOutput(form)

        # StructName
        inputObj = getattr(protocol, inputParams[0]).get()
        pdbFile, AS, addPointer = '', False, True
        if issubclass(type(inputObj), str):
            outStr = [inputObj]
            AS, addPointer = True, False
        elif issubclass(type(inputObj), AtomStruct):
            pdbFile = inputObj.getFileName()
            outStr = [os.path.splitext(os.path.basename(pdbFile))[0]]
            AS = True
        elif issubclass(type(inputObj), Sequence):
            outStr = [inputObj.getId().replace("|", "_")]

        if AS:
            # Chain
            chainJson = getattr(protocol, inputParams[1]).get()
            chainId = json.loads(chainJson)['chain']
        else:
            chainJson = ''

        # Positions
        posJson = getattr(protocol, inputParams[2]).get()
        if posJson:
            posIdxs = json.loads(posJson)['index']
            seq = json.loads(posJson)['residues']
            outStr += [posIdxs]
        else:
            outStr += ['FIRST-LAST']
            finalResiduesList = self.getResidues(form, inputObj, chainJson)
            idxs = [json.loads(finalResiduesList[0].get())['index'], json.loads(finalResiduesList[-1].get())['index']]
            seq = self.getSequence(finalResiduesList, idxs)

        chainStr, chainFileId = '', ''
        if AS:
          chainStr = ', "chain": "{}"'.format(chainId)
          chainFileId = '_{}'.format(chainId)

        prevStr = getattr(protocol, outputParam[0]).get()
        lenPrev = len(prevStr.strip().split('\n')) + 1
        if prevStr.strip() == '':
          lenPrev -= 1
        elif not prevStr.endswith('\n'):
          prevStr += '\n'

        val = getattr(protocol, 'entityType').get()
        if val==0: entity = 'protein'
        elif val==1: entity = 'dna'
        else: entity = 'rna'

        seqFile = protocol.getProject().getTmpPath('{}{}_{}.fa'.format(outStr[0], chainFileId, outStr[1]))
        with open(seqFile, 'w') as f:
            f.write('>{}\n{}\n'.format(outStr[0], seq))

        jsonStr = '%s) {"name": "%s"%s, "index": "%s", "seqFile": "%s", "entity": "%s"}\n' % \
                  (lenPrev, outStr[0], chainStr, outStr[1], seqFile, entity)
        form.setVar(outputParam[0], prevStr + jsonStr)


AddSequenceWizardChai().addTarget(protocol=ProtChai,
                              targets=['addInput'],
                              inputs=[{'inputOrigin': ['inputSequence', 'inputAtomStruct']},
                                      'inpChain', 'inpPositions'],
                              outputs=['inputList', 'inputPointers'])