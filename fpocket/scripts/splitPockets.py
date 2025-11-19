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

u = mda.Universe(pdb_file)
atoms = u.atoms
coords = atoms.positions

dist_matrix = distance_array(coords, coords)
condensed = dist_matrix[np.triu_indices(len(coords), k=1)]
Z = linkage(condensed, method="single")
clusters = fcluster(Z, t=distance_threshold, criterion="distance")

num_clusters = clusters.max()
print(f"Found {num_clusters} clusters.")

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
