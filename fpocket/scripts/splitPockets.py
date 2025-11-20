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
import sys
import os
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from scipy.cluster.hierarchy import fcluster, linkage
import numpy as np

pdb_file = sys.argv[1]
distance_threshold = float(sys.argv[2])
output_folder = sys.argv[3]

os.makedirs(output_folder, exist_ok=True)


# ------------------------------------------------------
#  PREPROCESS PDB: add MODEL headers, fix numbering/occupancy/B-factor
# ------------------------------------------------------
def preprocessPdbAddModels(inputPdb, outputPdb):
    """
    Rewrite PDB to:
    - Add MODEL X headers if missing (before each block ending with ENDMDL)
    - Renumber atoms sequentially within each model
    - Fix occupancy/B-factor
    """
    lines = []
    atomCounter = 1
    modelCounter = 1
    inside_model = False

    with open(inputPdb) as f:
        for line in f:
            if line.startswith("ENDMDL"):
                lines.append(line)
                inside_model = False
                continue

            if line.startswith(("ATOM", "HETATM")):
                if not inside_model:
                    # Start a new model
                    lines.append(f"MODEL        {modelCounter}\n")
                    modelCounter += 1
                    atomCounter = 1
                    inside_model = True

                # Ensure line is long enough
                line = line.rstrip().ljust(66)

                # Renumber atom
                line = f"{line[:6]}{atomCounter:5d}{line[11:]}"
                atomCounter += 1

                # Fix occupancy
                if line[54:60].strip() == "":
                    line = line[:54] + "  0.00 " + line[60:]
                # Fix B-factor
                if line[60:66].strip() == "":
                    line = line[:60] + "  0.00" + line[66:]

                lines.append(line + "\n")

    with open(outputPdb, "w") as f_out:
        f_out.writelines(lines)

    return outputPdb


# ------------------------------------------------------
#  SAFE LOAD OF PDB (handles corrupted PDBs gracefully)
# ------------------------------------------------------
corrected_pdb = os.path.join(output_folder, "corrected.pdb")
pdb_file = preprocessPdbAddModels(pdb_file, corrected_pdb)

try:
    u = mda.Universe(pdb_file)
except Exception as e:
    print("\n[ERROR] Could not read PDB. Probably malformed or corrupted.\n")
    print("MDAnalysis message:")
    print(e)

    marker = os.path.join(output_folder, "INVALID_PDB.txt")
    with open(marker, "w") as f:
        f.write("ERROR: PDB could not be read.\n")

    sys.exit(0)


# ------------------------------------------------------
#  EXTRACT COORDINATES
# ------------------------------------------------------
atoms = u.atoms
coords = atoms.positions

# ------------------------------------------------------
#  CLUSTERING
# ------------------------------------------------------
dist_matrix = distance_array(coords, coords)
condensed = dist_matrix[np.triu_indices(len(coords), k=1)]
Z = linkage(condensed, method="single")
clusters = fcluster(Z, t=distance_threshold, criterion="distance")

num_clusters = clusters.max()
print(f"Found {num_clusters} clusters.")

# ------------------------------------------------------
#  SAVE EACH POCKET
# ------------------------------------------------------
for i in range(1, num_clusters + 1):
    cluster_atoms = atoms[clusters == i]
    if len(cluster_atoms) == 0:
        continue

    cluster_atoms.occupancies[:] = 0.0
    cluster_atoms.tempfactors[:] = 0.0

    out_path = os.path.join(output_folder, f"pocket_{i}.pdb")
    with mda.Writer(out_path) as writer:
        writer.write(cluster_atoms)

    print(f"Saved: {out_path}  ({len(cluster_atoms)} atoms)")

    # ------------------------------------------------------
    #  CLEANUP: REMOVE 'SYST' AND FIX OCC/TEMPFACTOR
    # ------------------------------------------------------
    cleaned_lines = []
    with open(out_path) as f_in:
        for line in f_in:
            if line.startswith(("ATOM", "HETATM")):

                # Remove 'SYST' after column 54
                prefix = line[:54]
                rest = line[54:].replace("SYST", "    ")

                # Clean occupancy + tempFactor block
                numblock = rest[:12]
                numblock = "".join(
                    ch if ch in "0123456789.+- " else "0"
                    for ch in numblock
                ).ljust(12)

                rest = numblock + rest[12:]
                line = prefix + rest

            cleaned_lines.append(line)

    with open(out_path, "w") as f_out:
        f_out.writelines(cleaned_lines)