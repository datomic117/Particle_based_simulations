# DEM silo-discharge code, Assignment 4 of 6EMA02 Particle-based Simulations (group A4_8)

A discrete element code (linear spring-dashpot normal force, Cundall-Strack tangential force) in C, from the course's `DEM_student_v7.1.0` starting code. The report is `../A4.ipynb` (this folder is `code/` in the hand-in zip). The starter `README.txt` stays as delivered.

## Build and run

Plain C99 with the standard library; the only dependencies are `gcc` and `libm`. Run from `Assignment_4/`, one directory per run:

```
gcc -O3 code/*.c -o dem -lm                       # production, default settings of setparameters.c
gcc -g -Wall -Wextra code/*.c -o dem_debug -lm    # debugging
mkdir -p runs/A_1_v0.3_f0.05 && gcc -O3 -DRUN=1 -DNUM_TEST=1 -DV_IMPACT=0.3 -DDT_FRACTION=0.05 code/*.c -o runs/A_1_v0.3_f0.05/dem.exe -lm \
  && (cd runs/A_1_v0.3_f0.05 && ./dem.exe > run.dat)    # one Part A run
```

The settings that change from run to run (`RUN`, `NUM_TEST`, `V_IMPACT`, `DT_FRACTION`, `TILT_RATE`, `D_ORIFICE`, `T_DISCHARGE`) are macros at the top of `setparameters.c` and can be given with `-D`. The program prints a table (step, time, Z, Epot, Ekin, Etot, N_silo) to the screen and writes `trajectories.xyz` and `restart.dat` into the directory it runs in. Windows/MSYS2: `export PATH=/c/msys64/ucrt64/bin:$PATH` first.

## Source files

| file | content |
|---|---|
| `main.c` | set-up and the velocity-Verlet loop (unchanged from the starting code) |
| `setparameters.c/.h` | run selection, contact constants derived from $k_n$, restitution and contact time, time step (starting code, with the Part B and C stubs filled in by later parts) |
| `initialise.c/.h` | particle properties, positions and velocities; `initialise_collision` (Part A) |
| `forces.c/.h`, `dynamics.c/.h`, `nbrlist.c/.h`, `walls.c/.h` | contact forces, integrator, neighbour list and contact list, walls (starting code; later parts add to `walls.c`) |
| `fileoutput.c/.h`, `memory.c/.h`, `random.c/.h`, `structs.h`, `constants.h` | trajectory and restart files, allocation, random numbers, data structures (starting code) |
| `run_a.sh` | the 24 restitution runs of Part A (table) |
| `run_a_scan.sh` | 248 runs of Part A at 31 impact speeds (time-step explanation) |
| `analyse_a.py` | Part A analysis: $e$ from the energy table, writes `../data/A_*.csv` and `../data/A_A?_v0.3.dat` |

## What was added to the course code

* **Part A:** `initialise_collision` in `initialise.c`: equal radii $R_\mathrm{min}$ and the matching mass, zero spin; one particle at height $1.05R$ moving down (A1), or two particles on a line along $x$ in the box middle, centres $2\cdot1.05R$ apart, moving towards each other (A2), each with `v_impact`.

## How the data in `../data/` were made

| data | settings |
|---|---|
| `A_restitution.csv` | `bash code/run_a.sh` then `python code/analyse_a.py`: `RUN=1`, `NUM_TEST` 1 and 2, `V_IMPACT` 0.1, 0.3, 1.0, `DT_FRACTION` 0.2, 0.1, 0.05, 0.025 (24 runs, 2 ms each). $E_\mathrm{kin,before}$ = table value at the last step before the first step with $E_\mathrm{pot}>0$; $E_\mathrm{kin,after}$ = at the first step after the last step with $E_\mathrm{pot}>0$; $e=\sqrt{E_\mathrm{kin,after}/E_\mathrm{kin,before}}$. `phase` is the fraction of a time step between the first touch and the first grid time with overlap |
| `A_A1_v0.3.dat`, `A_A2_v0.3.dat` | the printed table of the 0.3 m/s, `DT_FRACTION` 0.05 runs (copied from `runs/A_?_v0.3_f0.05/run.dat`) |
| `A_phase_scan.csv` | `bash code/run_a_scan.sh` then `python code/analyse_a.py`: both tests, the four `DT_FRACTION`, `V_IMPACT` 0.10 to 1.00 in steps of 0.03 (248 runs), same definition of $e$ |

## Known limitations

* The Part A analysis uses the energy table at the grid times; $E_\mathrm{kin,before}$ needs at least one printed step before the first contact (checked by an assertion in `analyse_a.py`).
* `analyse_a.py` hard-codes $R_\mathrm{min}=0.975$ mm (from `setparameters.c`) to compute the contact phase; it asserts that the phase predicts the first overlap step of every run.
* Part A tests only the normal contact (head-on, no spin, no gravity); the tangential force and friction are not tested by it.
