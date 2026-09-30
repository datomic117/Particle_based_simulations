import numpy as np

filename = "b6_start.pdb"
L = 42.4612104

positions = []

with open(filename, "r") as f:
    for line in f:
        if line.startswith("HETATM"):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            positions.append([x, y, z])

positions = np.array(positions)

print("Number of sites:", len(positions))

sigma = {
    0: 3.75,   # CH3
    1: 3.95    # CH2
}

def site_type(i):
    p = i % 5
    return 0 if p in (0, 4) else 1

def excluded_same_molecule(i, j):
    # Different molecules always interact.
    if i // 5 != j // 5:
        return False

    # Same pentane molecule:
    # 1-2, 1-3 and 1-4 are excluded.
    separation = abs((i % 5) - (j % 5))

    if separation in (1, 2, 3):
        return True

    # 1-5, separation = 4, DOES interact.
    return False

minimum_r = float("inf")
minimum_pair = None

N = len(positions)

for i in range(N):
    for j in range(i + 1, N):

        if excluded_same_molecule(i, j):
            continue

        dr = positions[i] - positions[j]

        # Minimum-image convention
        dr -= L * np.rint(dr / L)

        r = np.linalg.norm(dr)

        if r < minimum_r:
            minimum_r = r
            minimum_pair = (i, j)

i, j = minimum_pair

ti = site_type(i)
tj = site_type(j)

sigma_ij = 0.5 * (sigma[ti] + sigma[tj])

names = {0: "CH3", 1: "CH2"}

print()
print("Closest interacting pair:")
print("particles:", i, j)
print("types:", names[ti], "-", names[tj])
print("distance =", minimum_r, "Angstrom")
print("sigma_ij =", sigma_ij, "Angstrom")
print("r / sigma_ij =", minimum_r / sigma_ij)