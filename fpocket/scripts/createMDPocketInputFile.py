# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


import fnmatch
import os,re,sys

if len(sys.argv)>2 :
	d=sys.argv[1]
	if not os.path.exists(d):
		sys.exit(" ERROR : the directory you provided as data directory does not exist. Breaking up.")
else :
	sys.exit(" ERROR : please provide the name of the data directory (where you have all your snapshots) and the name of the output file like : python createMDPocketInputFile.py data_dir myOutput.txt\nThe snapshots must have the pdb extension.")

prefix=sys.argv[1]
outputFile=sys.argv[2]

#just a helper function to do some stupid concatenation
def getFname(s):
    return prefix+os.sep+s


snapshots=fnmatch.filter(os.listdir(prefix),"*.pdb") #get the fpocket output folders

snapshots=[os.path.abspath(getFname(sn)) for sn in snapshots]

RE_DIGIT = re.compile(r'(\d+)')     #set a pattern to find digits
ALPHANUM_KEY = lambda s: [int(g) if g.isdigit() else g for g in RE_DIGIT.split(s)]      #create on the fly function (lambda) to return numbers in filename strings
snapshots.sort(key=ALPHANUM_KEY)      #sort by these numbers in filenames

fout=open(outputFile,"w")
[fout.write(sn+"\n") for sn in snapshots]
fout.close()

