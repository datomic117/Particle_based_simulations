# Dissipative particle dynamics (DPD) code, Assignment 3 of 6EMA02 Particle-based Simulations (group A3_6)

A dissipative particle dynamics code for beads and bead-spring chains in a periodic box, written in C.
It was converted from our Assignment 2 molecular dynamics code (itself built on the course's
MD_student v7.0.0 Lennard-Jones code) and is used to reproduce results of Groot and Warren,
*J. Chem. Phys.* **107**, 4423 (1997). Equation and figure numbers below refer to that paper.
The report is `../A3.ipynb` (this folder is `code/` in the hand-in zip).

## Build and run

Everything is plain C99 with the standard library; the only dependency is `gcc` and `libm`.

```
gcc -O3 code/*.c -o dpd -lm                      # production; run from Assignment_3/
gcc -g -Wall -Wextra code/*.c -o dpd_debug -lm   # debugging; run from Assignment_3/
./dpd > run.csv                                  # one line per step: step,time,Epot,Ekin,Etot
```

There are no input files. As delivered, `setparameters.c` is set for the D2 run ($N = 3$, $a_{AB} = 28$, output prefix `data/D2_N<N>`, trajectory `data/D2_N<N>.pdb`, restart file `restart_files/restart_D2_N<N>.dat`). **All settings live in `setparameters.c`**; a different run means editing
the values there and recompiling. The long runs of Parts C and D3 are done by the script
`run_c.sh`, which compiles one binary per run from a patched copy of `setparameters.c` (in
`build_c/<run>/`, deleted afterwards; `code/setparameters.c` is never modified) and runs them in
parallel:

```
export PATH=/c/msys64/ucrt64/bin:$PATH            # Windows/MSYS2: make gcc visible to bash
bash code/run_c.sh                                # default runs C1, C2, C3 (4x) and C4, 60 000 steps
RUNS="D3_N5_da3.5 D3_N10_da1" bash code/run_c.sh  # chain runs: D3_N<chain length>_da<Delta a>
STEPS=20000 BLOCK=2000 RUNS="C2_a37" bash code/run_c.sh   # shorter test run
```

Run names: `C<task>_a<a_AB>` (monomers) and `D3_N<N>_da<Delta a>` (chains, `a_AB = 25 + Delta a`).
The restart files are written to `restart_files/`. On a 20-core machine the seven monomer runs take
about 10 minutes together and the nine chain runs about 23 minutes together (chain initialisation
is slow, see "Known limitations").

## Units

Reduced DPD units: the cut-off $r_c = 1$ is the unit of length, the bead mass $m = 1$ the unit of
mass and $k_BT = 1$ the unit of energy. Time is in $r_c\sqrt{m/k_BT}$. All files in `data/` are in
these units. Density $\rho = 3$ (num_part = $3 L_xL_yL_z$), $a_{AA} = a_{BB} = 25$,
$\gamma = 4.5$, $\sigma = \sqrt{2\gamma k_BT}$, $\Delta t = 0.04$, spring $C = 2$ with rest length 0.

## Source files

| file | content |
|---|---|
| `main.c` | set-up, restart handling and the velocity-Verlet loop (unchanged from A2: the force routine is called after the first half kick and the position update, which is Eq. (9) with $\lambda = 1/2$) |
| `setparameters.c/.h` | **all run settings** (box, $a_{ij}$, $\gamma$, $\Delta t$, steps, force switches, chain length, start type, analysis switches, output names) |
| `forces.c/.h` | `calculate_forces`; `calculate_forces_nb`: DPD conservative, dissipative and random pair force; `calculate_forces_bond`: harmonic spring. Accumulates the energy and the (conservative + spring) virial, also per axis, for the pressure tensor |
| `dynamics.c/.h` | half-step velocity update (also gives the kinetic part of the pressure tensor per axis), position update, periodic boundary conditions. The Berendsen thermostat of A2 is no longer called |
| `initialise.c/.h` | `initialise_types` (A or B per chain), `initialise_bond_connectivity` (bead $k$ to $k+1$ within a chain), `initialise_positions` (random-walk chains with steps of the mean bond length, demixed start option), `initialise_velocities` |
| `nbrlist.c/.h` | cell list and Verlet neighbour list (unchanged from A2) |
| `analysis.c/.h` | on-the-fly analysis: $g(r)$, velocity histogram, chain size, composition profile along $x$ (two writers, see below), pressure tensor |
| `fileoutput.c/.h`, `memory.c/.h`, `random.c/.h`, `structs.h`, `constants.h` | PDB and restart files, allocation, uniform random numbers in $[-\sqrt3, \sqrt3]$ (`gauss`, mean 0 and variance 1), data structures |
| `force_test.c/.h` | finite-difference force and virial test (switch `force_test`; valid for the conservative force and the springs with the dissipative and random force off) |
| `run_c.sh` | script for the Part C and D3 production runs, see above |

## Switches in `setparameters.c`

| setting | meaning |
|---|---|
| `a[0]`, `a[1]` | $a_{AA} = a_{BB}$ and $a_{AB}$ |
| `conservative_force_on`, `dissipative_random_force_on` | the two switches of Part B (B1: random and dissipative off; B2: conservative off) |
| `chain_length` | beads per chain $N$; 1 means no bonds. `num_part` must be a multiple of it |
| `demix_initialisation` | 1: all A chains in $x < L_x/2$ and all B chains in $x > L_x/2$ (every demixing run) |
| `load_restart`, `reseed_velocities` | start from a restart file (B1 starts from the equilibrated monomer fluid) |
| `analysis_on`, `num_dt_phi`, `analysis_start_step` | $g(r)$ (bin 0.05, up to $r = 3$), sampling interval, steps skipped |
| `velocity_analysis_on` | velocity histogram and variance per type and component (B3) |
| `chain_size_analysis_on` | mean-square end-to-end distance (D1) |
| `profile_on`, `num_dt_profile`, `num_dt_pressure`, `num_dt_block` | composition profile along $x$ (bin 0.25) averaged per block, and the pressure tensor time series (C2 to C4, D3) |
| `composition_profile_on` | the second, older composition-profile writer used for D2 |
| `num_dt_output`, `num_dt_pdb`, `num_dt_restart` | output intervals in steps (one step is 0.04 time units) |

## Output of the analysis (`filename_analysis` + suffix)

| file | content |
|---|---|
| `<name>_rdf.csv` | $g(r)$ per pair (AA, AB, BB) with pair counts and expected counts; normalised by the ideal-gas expectation for the shell volume |
| `<name>_rdf_blocks.csv` | the same, per block of `num_dt_block` steps (convergence of $g(r)$, C1) |
| `<name>_profile.csv` | block, step, $x$, $\rho_A(x)$, $\rho_B(x)$ in bins of 0.25 per block of `num_dt_block` steps |
| `<name>_pressure.csv` | step, time, $P_{xx}$, $P_{yy}$, $P_{zz}$ of the whole box: $P_{\alpha\alpha} = [\sum m v_\alpha^2 + \sum r_{ij,\alpha}F_{ij,\alpha}]/V$ with the conservative and spring forces only |
| `<name>_velocity_hist.csv`, `<name>_velocity_stats.csv` | velocity histogram with the Maxwell-Boltzmann curve, and measured $k_BT$ per type and component |
| `<name>_chain_size.csv` | mean-square end-to-end distance of the chains |
| `<name>_composition_blocks.csv` | blocks of the older composition-profile writer: block, start step, end step, samples, x, rho_A, rho_B, phi_A |

The surface tension of C4 is $\sigma_s = \frac{L_x}{2}\big[P_{xx} - \tfrac12(P_{yy}+P_{zz})\big]$ (Eq. (31) integrated
over the whole periodic box, two interfaces), computed in the notebook from `_pressure.csv`.

## How the data in `../data/` were made

| data | settings |
|---|---|
| `B1_004.csv`, `B1_002.csv`, `B1_001.csv`, `B1_0005.csv` | dissipative and random force off, $N = 1$, $a_{AA} = a_{AB} = 25$, box $10^3$ (3000 beads), start from `restart_files/restart_monomer_fluid.dat`; $\Delta t$ = 0.04 / 0.02 / 0.01 / 0.005 with 200 / 400 / 800 / 1600 steps (all to $t = 8$) |
| `B2_ideal_gas.csv`, `B2_rdf_ideal_gas.csv` | conservative force off, $N = 1$, $a = 25$, box $10^3$, $\Delta t = 0.04$, 20 000 steps; `analysis_on = 1` for the $g(r)$ |
| `B3_monomer_fluid.csv` (energy log), `B3_velocity_off_*`, `B3_velocity_on_*` | $\Delta t = 0.04$, 20 000 steps, all forces on (conservative force off for `off`), $a = 25$, box $10^3$; velocity analysis on |
| `B2_analysis_rdf.csv` | energy log of the B2 run (not a $g(r)$ file despite its name; the notebook uses `B2_rdf_ideal_gas.csv`) |
| `C1_a25_*` | $N = 1$, box $10^3$, $a = 25$, fresh start, 20 000 steps; $g(r)$ sampled every 20 steps from step 4000, blocks of 4000 steps |
| `C2_a37_*`, `C3_a35_*`, `C3_a39_*`, `C3_a43_*`, `C3_a47_*`, `C4_a25_*` | $N = 1$, box $20 \times 10 \times 10$ (6000 beads), $a_{AB}$ = 37, 35, 39, 43, 47, 25 (control), demixed start, 60 000 steps, profile blocks of 5000 steps, pressure tensor every 10 steps; `_traj.pdb` frames every 10 000 (C2) or at the end |
| `D1_chain_size.csv`, `D1_chain_snapshot.pdb` | $N = 5$, $a = 25$ for all pairs, box $10^3$, 20 000 steps, all forces on, chain-size analysis on. The ideal reference is $3(N-1)k_BT/C$ (each bond has $\langle r^2\rangle = 3k_BT/C$) |
| `D2_N3/N5/N10_composition.csv`, `D2_composition_profiles.png` | box $20\times10\times10$, 6000 beads, $a_{AB} = 28$, $a_{AA} = 25$, demixed start, 60 000 steps, `chain_length` = 3, 5, 10, built separately. Averages of the trajectory frames after the first third, made by the D2 cell of the notebook |
| `D3_N<N>_da<Da>_*` | $N$ = 3 ($\Delta a$ = 3, 5, 7), 5 (2, 3.5, 5), 10 (1, 1.5, 2); box $20\times10\times10$, 6000 beads, demixed start, 60 000 steps, `run_c.sh` |

## What was added to the course code

* DPD pair forces (conservative, dissipative, random) in `calculate_forces_nb` with $\omega^D = (\omega^R)^2 = (1-r)^2$, $\sigma = \sqrt{2\gamma k_BT}$, one random number per pair and step, divided by $\sqrt{\Delta t}$, for both beads of the pair; potential energy and virial from the conservative force only.
* Type-dependent $a_{ij}$ and bead-spring chains of $N$ beads (`initialise_*`, springs with rest length 0, non-bonded interaction of all pairs including bonded ones: `factor_12_nb = factor_13_nb = factor_14_nb = 1`); angle and torsion forces and the Berendsen thermostat are removed.
* Initialisation of random-walk chains, optionally with all A chains in the left and all B chains in the right half of the box.
* Analysis: $g(r)$ per pair, velocity statistics, chain size, composition profile in blocks, diagonal pressure tensor.
* `run_c.sh` for the production runs.

## Known limitations

* Two composition-profile writers exist in `analysis.c` (`composition_profile_on`, used for D2, and `profile_on`, used for C and D3). They give the same quantity; the second writes block files that the notebook uses for the error estimates.
* The initialisation of chains places beads one by one with a rejection test, which is slow for chain length 10 (about a minute).
* The D2 cell of the notebook reads the trajectory files `data/D2_N<N>.pdb`; if these are not present it uses the saved `D2_N<N>_composition.csv`.
* `force_test` does not test the dissipative and random forces, which derive from no potential; these are checked with the energy, $g(r)$ and velocity tests of Part B.

## Who did what

* Job Scheepens (jobscheepens): DPD forces, chains and initialisation, $g(r)$, Part B, D1 and D2.
* Tomas Dadikozyan (Code-Tomas-Dadikozyan): Part C (C1 to C4), D3, report text and housekeeping sections.
* Claude Code (Anthropic) assisted with the Part C and D3 analysis code and text, as declared in the Sources section of the notebook.
