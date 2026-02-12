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

from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin
from .constants import *

_references = ['']


class Plugin(pwchemPlugin):
    @classmethod
    def defineBinaries(cls, env):
        cls.addBoltzPackage(env)
        cls.addChaiPackage(env)
        cls.addProtenixPackage(env)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(BOLTZ_DIC['home'], cls.getEnvName(BOLTZ_DIC))
        cls._defineEmVar(CHAI_DIC['home'], cls.getEnvName(CHAI_DIC))
        cls._defineEmVar(PROTENIX_DIC['home'], cls.getEnvName(PROTENIX_DIC))

    @classmethod
    def addBoltzPackage(cls, env, default=True):
        installer = InstallHelper(
            BOLTZ_DIC['name'],
            packageHome=cls.getVar(BOLTZ_DIC['home']),
            packageVersion=BOLTZ_DIC['version']
        )

        installer.getCondaEnvCommand(
            BOLTZ_DIC['name'],
            binaryVersion=BOLTZ_DIC['version'],
            pythonVersion='3.11'
        ).addCommand(
            f"{cls.getEnvActivationCommand(BOLTZ_DIC)} && "
            "git clone --branch v2.2.1 --depth 1 https://github.com/jwohlwend/boltz.git && "
            "pip install --editable ./boltz[cuda]",
            f"{BOLTZ_DIC['name']}_installed"
        )

        installer.addPackage(
            env,
            dependencies=['conda', 'pip', 'git'],
            default=default
        )

    @classmethod
    def addProtenixPackage(cls, env, default=True):
        installer = InstallHelper(
            PROTENIX_DIC['name'],
            packageHome=cls.getVar(PROTENIX_DIC['home']),
            packageVersion=PROTENIX_DIC['version']
        )

        installer.getCondaEnvCommand(
            PROTENIX_DIC['name'],
            binaryVersion=PROTENIX_DIC['version'],
            pythonVersion='3.11'
        ).addCommand(
            f"{cls.getEnvActivationCommand(PROTENIX_DIC)} && "
            "conda install -y -c nvidia cuda-toolkit=12.1 && "
            "echo 'export CUDA_HOME=$CONDA_PREFIX' >> $CONDA_PREFIX/etc/conda/activate.d/protenix_cuda.sh && "
            
            "mkdir -p $CONDA_PREFIX/etc/conda/activate.d && "
            "echo 'export PATH=$CUDA_HOME/bin:$PATH' >> $CONDA_PREFIX/etc/conda/activate.d/protenix_cuda.sh && "
            "echo 'export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH' >> $CONDA_PREFIX/etc/conda/activate.d/protenix_cuda.sh && "
            "echo 'export CPATH=$CUDA_HOME/include:$CPATH' >> $CONDA_PREFIX/etc/conda/activate.d/protenix_cuda.sh && "

            "git clone --branch v1.0.4 --depth 1 https://github.com/bytedance/Protenix.git && "
            "cd Protenix && "
            "pip install -e . ",
            f"{PROTENIX_DIC['name']}_installed"
        )

        installer.addPackage(
            env,
            dependencies=['conda', 'pip', 'git'],
            default=default
        )

    @classmethod
    def addChaiPackage(cls, env, default=True):
        installer = InstallHelper(
            CHAI_DIC['name'],
            packageHome=cls.getVar(CHAI_DIC['home']),
            packageVersion=CHAI_DIC['version']
        )

        installer.getCondaEnvCommand(
            CHAI_DIC['name'],
            binaryVersion=CHAI_DIC['version'],
            pythonVersion='3.11'
        ).addCommand(
            f"{cls.getEnvActivationCommand(CHAI_DIC)} && "
            "pip install chai_lab==0.6.1",
            f"{CHAI_DIC['name']}_installed"
        )

        installer.addPackage(
            env,
            dependencies=['conda', 'pip', 'git'],
            default=default
        )




