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
import os

_references = ['']


class Plugin(pwchemPlugin):
    @classmethod
    def defineBinaries(cls, env):
        cls.addBoltzPackage(env)
        cls.addChaiPackage(env)
        cls.addIntelliFoldPackage(env)
        cls.addProtenixPackage(env)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(BOLTZ_DIC['home'], cls.getEnvName(BOLTZ_DIC))
        cls._defineEmVar(CHAI_DIC['home'], cls.getEnvName(CHAI_DIC))
        cls._defineEmVar(INTELLIFOLD_DIC['home'], cls.getEnvName(INTELLIFOLD_DIC))
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
            # boltz only pins torch>=2.2, so installing boltz[cuda] first pulls the
            # latest wheel (currently cu130), which needs a very recent driver
            # (CUDA 13) and fails CUDA init on common 12.x drivers. Install a
            # CUDA 12.4 torch first -- it runs on any driver supporting CUDA >=12.4
            # -- so boltz then sees torch>=2.2 satisfied and keeps this build.
            "pip install torch --index-url https://download.pytorch.org/whl/cu124 && "
            "pip install --editable ./boltz[cuda]",
            f"{BOLTZ_DIC['name']}_installed"
        )

        installer.addPackage(
            env,
            dependencies=['conda', 'pip', 'git'],
            default=default
        )

    @classmethod
    def addIntelliFoldPackage(cls, env, default=True):
        installer = InstallHelper(
            INTELLIFOLD_DIC['name'],
            packageHome=cls.getVar(INTELLIFOLD_DIC['home']),
            packageVersion=INTELLIFOLD_DIC['version']
        )

        installer.addCommand(
            "git clone --branch v2.0.2 --depth 1 https://github.com/IntelliGen-AI/IntelliFold.git",
            f"{INTELLIFOLD_DIC['name']}_cloned"
        ).addCommand(
            f"cd IntelliFold && conda env create -f environment.yaml -n {INTELLIFOLD_DIC['name']}-{INTELLIFOLD_DIC['version']}",
            f"{INTELLIFOLD_DIC['name']}_env_created"
        ).addCommand(
            f"conda install -n {INTELLIFOLD_DIC['name']}-{INTELLIFOLD_DIC['version']} -c nvidia cuda-toolkit=12.4 -y",
            f"{INTELLIFOLD_DIC['name']}_cuda_installed"
        ).addCommand(
            f"cd IntelliFold && conda run -n {INTELLIFOLD_DIC['name']}-{INTELLIFOLD_DIC['version']} pip install -e . && touch ../intellifold_installed",
            f"{INTELLIFOLD_DIC['name']}_installed"
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

    @classmethod
    def addProtenixPackage(cls, env, default=True):
        installer = InstallHelper(
            PROTENIX_DIC['name'],
            packageHome=cls.getVar(PROTENIX_DIC['home']),
            packageVersion=PROTENIX_DIC['version']
        )

        pkgHome = cls.getVar(PROTENIX_DIC['home'])
        commonDir = os.path.join(pkgHome, "common")

        installer.getCondaEnvCommand(
            PROTENIX_DIC['name'],
            binaryVersion=PROTENIX_DIC['version'],
            pythonVersion='3.11'
        ).addCommand(
            f"{cls.getEnvActivationCommand(PROTENIX_DIC)} && "
            "pip install protenix==2.0.0 && "
            "conda install -y -c nvidia cuda-nvcc cuda-toolkit && "
            f"mkdir -p {commonDir} && "
            f"wget -P {commonDir} https://files.wwpdb.org/pub/pdb/data/monomers/components.cif && "
            f"python -c \"import pickle; from rdkit import Chem; from protenix.data.core.ccd import biotite_load_ccd_cif; print('Processing CCD... this takes a minute'); " 
            f"from protenix.data.core.ccd import get_all_ccd_code;\" && "
            "export PATH=$CONDA_PREFIX/bin:$PATH && "
            "export CUDA_HOME=$CONDA_PREFIX && "
            "export CPATH=$CONDA_PREFIX/include:$CPATH && "
            "export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH && "
            f"touch {pkgHome}/{PROTENIX_DIC['name']}_installed",
            f"{PROTENIX_DIC['name']}_installed"
        )


        installer.addPackage(
            env,
            dependencies=['conda', 'pip', 'git'],
            default=default
        )




