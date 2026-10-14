# MD_student_v7.0.0, extended to an n-pentane simulator

Molecular dynamics of 400 united-atom n-pentane molecules (2000 sites, CH3-CH2-CH2-CH2-CH3) in a
periodic cubic box at 626 kg/m^3 and 293 K, for Assignment 2 of 6EMA02 Particle-based Simulations
(group 5). The starting point is the course code (Lennard-Jones fluid); everything below describes
what the code does *now*. The report is `../A2.ipynb` (this folder is `code/` in the hand-in zip).

## Build and run

```
gcc -O3 code/*.c -o md -lm                      # production; run from  Assignment_3/
gcc -g -Wall -Wextra code/*.c -o md_debug -lm   # debugging; run from Assignment_3/
./md > run.csv
```


## Units



## Source files



## How the data in `../data/` were made
All of B1: Random & Dissipative force = 0, chain_length = 1, aAA = aAB = 25, box = 10^3
B1_0005.csv: dt = 0.005 & num_steps = 1600
B1_001.csv: dt = 0.01 & num_steps = 800
B1_002.csv: dt = 0.02 & num_steps = 400
B1_004.csv: dt = 0.04 & num_steps = 200

All of B2: chain_length = 1, Conservative force = 1, aAA = aAB = 25, box-10^3
B2_ideal_gas.csv: dt = 0.04 & num_steps = 20000
B2_rdf_ideal_gas.csv: dt = 0.04, num_steps = 20000 and analysis_on = 1

All of B3: dt = 0.04, num_steps = 20000, all forces on, aAA = aAB = 25, box = 10^3
B3_monomer_fluid.csv: see above
B3_velocity_off_...: conservative forces off, velocity_analysis off
B3_velocity_on_...: conservative forces on, velocity_analysis on





D1: dt = 0.04, num_steps = 20000, all forces on, aAA = aAB = 25
Chain_size_analysis ON
D1_chain_size.csv: the average of the squared end-to-end vector of every chain
The ideal chain reference is: 3* (N-1)*(kT/kb) as each harmonic bond has r^2 = 3kT/kb when the rest length = 0
D1_chain_size.pdb: preserves the final state for reference in OVITO

D2: box: 20x10x10, 6000 beads, num_steps = 60000, aAB = 28, aAA = 25, all forces on, demixed start.
Build separately with chain_length = 3, 5, and 10. The PDB snapshots/trajectories and restart files use chain-length-specific names.
`composition_profile_on` writes `data/D2_N<N>_composition_blocks.csv` every 5000 steps, averaging positions sampled every 100 steps. Each block has 80 x-slabs (width 0.25) and contains x, rho_A, rho_B, and phi_A, so equilibration and block-to-block variation can be inspected.
After all three runs, execute `python plot_d2_profiles.py` from `Assignment_3/`. It averages the post-20,000-step blocks and writes the combined figure to `data/D2_composition_profiles.png`.


## What was added to the course code



## Who did what
