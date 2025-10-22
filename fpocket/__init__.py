# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

from os.path import join, exists

import pwem
from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin
from .constants import *

_version_ = '0.1'
_logo = "fpocket_logo.png"
_references = ['']


class Plugin(pwchemPlugin):
    _homeVar = FPOCKET_DIC['home']
    _pathVars = [FPOCKET_DIC['home']]
    _supportedVersions = [FPOCKET_DIC['version']]

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(FPOCKET_DIC['home'], FPOCKET_DIC['name'] + '-' + FPOCKET_DIC['version'])

    @classmethod
    def defineBinaries(cls, env, default=True):
        installer = InstallHelper(FPOCKET_DIC['name'], packageHome=cls.getVar(FPOCKET_DIC['home']),
                                  packageVersion=FPOCKET_DIC['version'])

        fpocketPath = join(pwem.Config.EM_ROOT, cls.getEnvName(FPOCKET_DIC))
        installer.addCommand(f'conda create -y -c conda-forge fpocket -p {fpocketPath}',
                             f'{FPOCKET_DIC["name"]}_installed').\
            addPackage(env, dependencies=['conda'], default=default)

    @classmethod
    def runFpocket(cls, protocol, program, args, cwd=None):
        """ Run Fpocket command from a given protocol. """
        protocol.runJob(join(cls.getVar(FPOCKET_DIC['home']), 'bin/{}'.format(program)), args, cwd=cwd)

    # ---------------------------------- Utils functions  -----------------------

