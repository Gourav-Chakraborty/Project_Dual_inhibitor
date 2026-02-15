import os
os.environ["OPENMM_PLUGIN_DIR"] = ""
import MDAnalysis as mda
import numpy as np
from tqdm import tqdm

# ======== USER INPUTS ========
topology_file = "analysis.prmtop"
trajectory_file = "analysis.nc"
output_npy = "ort_inter_protein_distance.npy"
output_pairs = "ort_residue_pairs.npy"
n_sample_frames = 2000
resid_prot1 = (1, 229)     # 14-3-3ζ
resid_prot2 = (230, 326)   # α-synuclein
# ==============================

# Load trajectory
u = mda.Universe(topology_file, trajectory_file)
n_total_frames = len(u.trajectory)

# Select CA atoms from both proteins
protein1_ca = u.select_atoms(f"protein and name CA and resid {resid_prot1[0]}-{resid_prot1[1]}")
protein2_ca = u.select_atoms(f"protein and name CA and resid {resid_prot2[0]}-{resid_prot2[1]}")

n_res1 = len(protein1_ca)
n_res2 = len(protein2_ca)

# Sanity checks
assert n_res1 == (resid_prot1[1] - resid_prot1[0] + 1), f"Expected {resid_prot1[1] - resid_prot1[0] + 1} residues for protein 1, got {n_res1}"
assert n_res2 == (resid_prot2[1] - resid_prot2[0] + 1), f"Expected {resid_prot2[1] - resid_prot2[0] + 1} residues for protein 2, got {n_res2}"

# Compute evenly spaced frame indices
frame_indices = np.linspace(0, n_total_frames - 1, n_sample_frames, dtype=int)
assert len(frame_indices) == n_sample_frames, "Sampling frame count mismatch."

# Pre-allocate distance matrix
inv_distances = np.zeros((n_sample_frames, n_res1 * n_res2), dtype=np.float32)

# Create residue pair list
res_pairs = [(i, j) for i in range(resid_prot1[0], resid_prot1[1] + 1)
                   for j in range(resid_prot2[0], resid_prot2[1] + 1)]
np.save(output_pairs, res_pairs)
print(f"Saved residue-pair index map to: {output_pairs} (length: {len(res_pairs)})")

# Compute inverse distances
print(f"\nProcessing {n_sample_frames} evenly spaced frames from {n_total_frames} total frames...")
for frame_idx, ts in enumerate(tqdm(frame_indices)):
    u.trajectory[ts]  # go to the selected frame

    dists = np.linalg.norm(
        protein1_ca.positions[:, None, :] - protein2_ca.positions[None, :, :],
        axis=-1
    )

    with np.errstate(divide='ignore'):
        inv_dists = 1.0 / dists
        inv_dists[np.isinf(inv_dists)] = 0.0

    inv_distances[frame_idx, :] = inv_dists.flatten()

# Save final distance matrix
np.save(output_npy, inv_distances)
print(f"\nSaved inverse distances to: {output_npy}")
print(f"Shape: {inv_distances.shape} (frames × features)")


