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
        cls.addIntelliFoldPackage(env)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(BOLTZ_DIC['home'], cls.getEnvName(BOLTZ_DIC))
        cls._defineEmVar(CHAI_DIC['home'], cls.getEnvName(CHAI_DIC))
        cls._defineEmVar(INTELLIFOLD_DIC['home'], cls.getEnvName(INTELLIFOLD_DIC))

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
    def addIntelliFoldPackage(cls, env, default=True):
        installer = InstallHelper(
            INTELLIFOLD_DIC['name'],
            packageHome=cls.getVar(INTELLIFOLD_DIC['home']),
            packageVersion=INTELLIFOLD_DIC['version']
        )

        installer.addCommand(
            "git clone --branch v2.0.0 --depth 1 https://github.com/IntelliGen-AI/IntelliFold.git",
            f"{INTELLIFOLD_DIC['name']}_cloned"
        ).addCommand(
            f"cd IntelliFold && conda env create -f environment.yaml -n {INTELLIFOLD_DIC['name']}-{INTELLIFOLD_DIC['version']}",
            f"{INTELLIFOLD_DIC['name']}_env_created"
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




