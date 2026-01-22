import os

from pwem.objects import SetOfAtomStructs, AtomStruct
from pwem.viewers import Chimera
from pyworkflow.viewer import Viewer

from biofold.protocols import ProtAlphaFold3


class ProtAlphaFold3Viewer(Viewer):
    """ Viewer for ChimeraProtAlphaFold3 protocol output. """
    _label = 'viewer discrepancies'
    _targets = [AtomStruct]

    def visualize(self, obj, **args):
        # Create Chimera command file
        fnCmd = self.protocol._getExtraPath("chimera_output.cxc")
        with open(fnCmd, 'w') as f:
            # Process protocol outputs
            for output in self.protocol._outputs:
                # If the file is an atomic structure (.cif or .pdb), open it in Chimera
                fileName = os.path.abspath(eval(f'self.protocol.{output}.getFileName()'))
                if fileName.endswith(".cif") or fileName.endswith(".pdb"):
                    f.write(f"open {fileName}\n")
                    f.write("color bfactor palette alphafold\n")

        # Run Chimera with the generated command file
        Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")