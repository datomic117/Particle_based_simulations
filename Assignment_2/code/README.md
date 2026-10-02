# MD_student_v7.0.0, extended to an n-pentane simulator

Molecular dynamics of 400 united-atom n-pentane molecules (2000 sites, CH3-CH2-CH2-CH2-CH3) in a
periodic cubic box at 626 kg/m^3 and 293 K, for Assignment 2 of 6EMA02 Particle-based Simulations
(group 5). The starting point is the course code (Lennard-Jones fluid); everything below describes
what the code does *now*. The report is `../A2.ipynb` (this folder is `code/` in the hand-in zip).

## Build and run

```
gcc -O3 *.c -o md -lm                      # production
gcc -g -Wall -Wextra *.c -o md_debug -lm   # debugging
./md > run.csv                             # run from this folder
```

- **Windows:** use MSYS2 (`gcc` is in `C:\msys64\ucrt64\bin`) or the VS Code tasks in `.vscode/`.
  PowerShell redirection (`>`) writes UTF-16 files; the notebook's `load_csv` reads both encodings.
- **There are no input files.** All run settings are in `setparameters.c`; change a value and recompile.
- **Output.** The program prints a CSV to standard output every `num_dt_output` steps:
  `step,time,Epot,Ekin,Etot,T,T_CH3_K,T_CH2_K` (internal units, see below). `.pdb` trajectories and
  binary restart files go to the current folder / `restart_files/`.
- The `setparameters.c` in this folder holds the settings of the last B9 run (NVT, tau = 1 ps,
  20 ps, `analysis_on = 0`). The C1/C2 runs are made by `run_c.sh`, which patches a copy of it.

## Units

Energy unit k_B x 1 K, length 1 Angstrom, mass 1 amu; time follows: 1 internal time unit = 1.0967 ps,
so the 1 fs time step is `dt = 0.0009118`. Temperature in kelvin equals the energy per k_B.
Derivation in the notebook (A1).

## Source files

| File | What it does | Notebook task |
|---|---|---|
| `main.c` | control flow: setup, velocity-Verlet loop, CSV output, finite-difference test mode, torsion net force/torque check | B1, B5, B7 |
| `setparameters.c` | **all** run settings: force-field parameters, box, time step, thermostat, run length, output, restart, analysis | all |
| `constants.h`, `structs.h` | constants, data types (`Parameters`, `Vectors`, ...) and the `v3_*` vector helpers (the `Analysis` state is declared in `analysis.h`) | |
| `memory.c` | allocation and freeing of all arrays | |
| `initialise.c` | site types (CH3/CH2), topology (bonds, angles, dihedrals, exclusions), molecule packing, Maxwell-Boltzmann velocities | B2, B6 |
| `nbrlist.c` | cell-linked list and Verlet neighbour list (1-2, 1-3, 1-4 pairs are excluded) | |
| `forces.c` | LJ with Lorentz-Berthelot mixing and shifted cut-off; bond, angle and Ryckaert-Bellemans dihedral forces; virials | B2, B4 |
| `dynamics.c` | velocity-Verlet halves (type-dependent masses), periodic boundaries, Berendsen thermostat | B2, B8 |
| `force_test.c` | finite-difference test of forces and virial (`force_test > 0` in `setparameters.c`) | B3, B5 |
| `analysis.c`, `analysis.h` | on-the-fly analysis: dihedral histograms and residence times (C1), mean-square displacement of molecular centres of mass (C2) | C1, C2 |
| `fileoutput.c` | `.pdb`/`.xyz` trajectories, restart files | |
| `random.c` | uniform and (sum of 12 uniforms) Gaussian random numbers | |
| `run_c.sh` | produces the C1/C2 data (see below) | C1, C2 |
| `tools/` | small helper scripts, see below | B6, B7 |
| `restart_files/` | binary states (positions, velocities, forces) saved by the runs; `restart_b9_tau001.dat` is the start of the C runs | |
| `html/`, `Doxyfile`, `mainpage.dox`, `DoxygenLayout.xml` | Doxygen documentation of the course code (open `html/index.html`); generated before `analysis.c` was added, rebuild with `doxygen Doxyfile` | |
| `b6_render.png` | rendering of the packed start configuration (OVITO), used in B6 | B6 |

## How the data in `../data/` were made

All runs: `gcc -O3 *.c -o md -lm`, edit `setparameters.c`, `./md > ../data/NAME.csv`.
The exact settings of the B1-B9 runs were not saved (see the notebook's *Build & run log*); the table
gives the type of each run.

| Data | Run |
|---|---|
| `b3_*` | 2000 sites, NVE at 626 kg/m^3, 50 ps with dt = 1 fs and 0.5 fs; `b3_force_test.txt` is `force_test` output |
| `b5_*` | one molecule (`num_part = 5`): thermalisation, finite-difference tests, 100 ps NVE |
| `b7_nve.csv`, `b9_*` | 400 molecules from the packed start (B6), 20 ps: NVE, NVT with tau = 0.01 ps and 1 ps |
| `c_equil_k*`, `c_prod_k*` | C1/C2, written by `run_c.sh` (below) |

### C1/C2 production (`run_c.sh`)

```
NREP=4 EQUIL_STEPS=100000 PROD_STEPS=300000 bash run_c.sh      # run from this folder, ~100 min
```

For each of 4 independent replicas it builds two binaries from a patched copy of `setparameters.c`
(which is restored afterwards) and runs them in parallel:

1. **Equilibration**, 100 ps NVT (Berendsen, tau = 1 ps), from `restart_files/restart_b9_tau001.dat`
   with new velocities (`reseed_velocities`). Writes `../data/c_equil_k.csv` and the dihedral
   histograms in 10 ps blocks (`c_equil_k_hist2d.csv`, ...). Restart output: `restart_files/restart_c_equil_k.dat`.
2. **Production**, 300 ps NVE from that restart. Writes `../data/c_prod_k.csv` and the analysis files:

| File (`c_prod_k_...`) | Content |
|---|---|
| `hist2d.csv` | joint histogram of (phi1, phi2), 72 x 72 bins of 5 degrees, per 100 ps block |
| `dwell.csv` | samples in and transitions out of the trans and gauche states (mean residence times) |
| `phi_trace.csv` | phi1 and phi2 of molecule 0 against time |
| `msd.csv` | MSD of the molecular centres of mass against lag (0-100 ps), total and per block |

`restart_c_*.dat` are regenerated by the script and are not kept.

The notebook discards the first 100 ps of each production run (the CH3 ends are still hotter than the CH2 interior there) and uses the remaining 200 ps per replica; the files themselves contain all 300 ps.

## Analysis (`analysis.c`)

- **Unwrapped positions.** `boundary_conditions` wraps positions into the box every step, so the
  MSD uses a second array `r_unwrapped` that accumulates the per-step displacement `dr`.
- **Dihedrals.** IUPAC angle phi from the minimum-image bond vectors (`dihedral_angle`):
  trans = +-180 degrees, g+ = +60. States: trans for |phi| > 120, g+ for 0 < phi < 120, g- otherwise.
  Sampled every `num_dt_phi` steps (0.1 ps).
- **MSD.** Mass-weighted molecular centres of mass from `r_unwrapped`; sampled every `num_dt_msd` steps
  (10 fs); a new time origin every `num_dt_msd_origin` steps (1 ps), the last 100 origins are kept,
  so lags go up to 100 ps. Averaged over the 400 molecules and over all origins.
- **Parameters** (all in `setparameters.c`): `analysis_on`, `msd_on`, `num_dt_phi`, `num_dt_msd`,
  `num_dt_msd_origin`, `num_dt_block` (steps per error-estimation block), `filename_analysis` (path prefix).
  Files are rewritten at every block boundary and at the end.
- **Replicas.** `reseed_velocities > 0` redraws the velocities after loading a restart file, which
  gives independent replicas of the same equilibrated state.

## Helper scripts (`tools/`)

Run from inside `tools/`; they are not needed to build or run the simulation.

- `check_b6.py`: closest interacting pair of the B6 start configuration, from `../b6_start.pdb`
  (written by the B6 run; not kept in the zip).
- `check_b7.py`: drift and transient statistics of `../../data/b7_nve.csv`.
- `render_b6.py`: renders `b6_start.pdb` with OVITO's Python module into `../b6_render.png`.

## What was added to the course code

The course code is MD_student_v7.0.0 (Dr. Ir. E.A.J.F. Peters, TU/e). Our work: typed masses and LJ
(B2), bonded forces (B4), `check_torsion_force_and_torque` and the per-term test wrappers in `main.c`
(B5), molecule packing and typed Maxwell-Boltzmann velocities (B6), the Berendsen thermostat (B8),
`analysis.c` and the hooks for it in `main.c` (C1, C2), and a per-call table of the LJ pair parameters
in `calculate_forces_nb`. The notebook's *Sources* section states which parts were written with AI help.

## Who did what

Olivier: A1, A2, B1-B5, C2. Thomas: B6-B9, C1.
