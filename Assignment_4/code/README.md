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
| `extract_frame.py` | cuts one frame (by time or index) out of a `trajectories.xyz` into a small text file (x y z radius speed); used by Part B, reused by C and D |
| `render_ovito.py` | OVITO snapshot of one frame of a `trajectories.xyz`: slab slice through the axis, coloured by speed with a legend, Tachyon renderer, optional crop; used by Part C, reused by D (needs `pip install ovito`) |
| `tests/wall_check.c` | standalone check of `cylinder_wall` and `silo_bottom_wall` against a brute-force closest point; links the real `walls.c`, has its own `main`, so it lives in a subdirectory and `gcc code/*.c` does not see it |
| `run_c_noedge.sh` | Part C2: copies the sources to `runs/C_noedge/src`, switches the edge of the hole off in the copy (`return false`) and runs the silo; `code/` itself is not changed |

## What was added to the course code

* **Part A:** `initialise_collision` in `initialise.c`: equal radii $R_\mathrm{min}$ and the matching mass, zero spin; one particle at height $1.05R$ moving down (A1), or two particles on a line along $x$ in the box middle, centres $2\cdot1.05R$ apart, moving towards each other (A2), each with `v_impact`.
* **Part B:** the `RUN_SLOPE` branch of `update_parameters` in `setparameters.c`: the slope angle is zero until `t_tilt`, then `tilt_rate * (time - t_tilt)`; it is stored in `slope_angle` (rad) and gravity is set to $g(\sin\theta, 0, -\cos\theta)$. The walls do not move. `extract_frame.py` (see above).

* **Part C:** `cylinder_wall` (C1) and the edge of the hole in `silo_bottom_wall` (C2) in `walls.c`: the closest wall point is found in closed form, the rim circle for the hole, the cylinder surface at the height of the particle (derivation in the notebook, check in `tests/wall_check.c`). `update_parameters` (RUN_SILO branch) sets `orifice_open = 1` from `t_open` on. `initialise_positions_silo` (C3) in `initialise.c`: loose square lattice of spacing 2.2 R_max inside R_cyl - 1.1 R_max, layer by layer from z_silo + 1.1 R_max, shifts of at most 0.04 R_max per coordinate. `count_in_silo` (Task D2, one loop over the particles with r.z > z_silo) was also written in Part C, because the table of the silo run is the data of Part D and must contain `N_silo`. `render_ovito.py` (see above).

## How the data in `../data/` were made

| data | settings |
|---|---|
| `A_restitution.csv` | `bash code/run_a.sh` then `python code/analyse_a.py`: `RUN=1`, `NUM_TEST` 1 and 2, `V_IMPACT` 0.1, 0.3, 1.0, `DT_FRACTION` 0.2, 0.1, 0.05, 0.025 (24 runs, 2 ms each). $E_\mathrm{kin,before}$ = table value at the last step before the first step with $E_\mathrm{pot}>0$; $E_\mathrm{kin,after}$ = at the first step after the last step with $E_\mathrm{pot}>0$; $e=\sqrt{E_\mathrm{kin,after}/E_\mathrm{kin,before}}$. `phase` is the fraction of a time step between the first touch and the first grid time with overlap |
| `A_A1_v0.3.dat`, `A_A2_v0.3.dat` | the printed table of the 0.3 m/s, `DT_FRACTION` 0.05 runs (copied from `runs/A_?_v0.3_f0.05/run.dat`) |
| `A_phase_scan.csv` | `bash code/run_a_scan.sh` then `python code/analyse_a.py`: both tests, the four `DT_FRACTION`, `V_IMPACT` 0.10 to 1.00 in steps of 0.03 (248 runs), same definition of $e$ |
| `B_rate20.dat`, `B_rate10.dat` | the printed tables of `-DRUN=2` (20 °/s) and `-DRUN=2 -DTILT_RATE=10`, copied from `runs/B_rate20/run.dat` and `runs/B_rate10/run.dat` (one directory per run, as in Part A). Columns: step, time, $Z$, $E_\mathrm{pot}$, $E_\mathrm{kin}$, $E_\mathrm{tot}$, $N_\mathrm{silo}$ |
| `B_frame_20deg.dat` | `python code/extract_frame.py runs/B_rate20/trajectories.xyz data/B_frame_20deg.dat --time 1.3`: the frame nearest to $t = 1.3$ s (frame 59, $t = 1.2917$ s, slope 19.83°), columns x y z radius speed |
| `D7.dat` | the printed table of the silo run `-DRUN=3` (defaults: 7500 particles, $D = 7\,d_p$, 2.3 s, the orifice opens at 0.3 s), copied from `runs/D7/run.dat`. The base of Part C and Part D |
| `C_frame_before_open.dat` | `python code/extract_frame.py runs/D7/trajectories.xyz data/C_frame_before_open.dat --frame 13`: the last frame before $t_\mathrm{open}$ (frame 13, $t = 0.2846$ s), columns x y z radius speed |
| `C_snapshot_settled.png` | `python code/render_ovito.py runs/D7/trajectories.xyz data/C_snapshot_settled.png --time 0.286 --zrange 51 113 --xhalf 20 --vmax 0.001 --size 800 1000`: slab of 4 mm through the axis at frame 13, speed scale 0 to 1 mm/s |
| `C_wall_check.txt` | `gcc -O2 -Wall -Wextra -Icode code/tests/wall_check.c code/walls.c -o runs/_check/wall_check.exe -lm && runs/_check/wall_check.exe > data/C_wall_check.txt` (30 s): about 8e5 random particles against a brute-force closest point |
| `C_noedge_D7.dat` | `bash code/run_c_noedge.sh`, then `cp runs/C_noedge/run.dat data/C_noedge_D7.dat`: the table of the silo run with the edge of the hole switched off (used in the C2 answer only) |

## Known limitations

* The Part A analysis uses the energy table at the grid times; $E_\mathrm{kin,before}$ needs at least one printed step before the first contact (checked by an assertion in `analyse_a.py`).
* `analyse_a.py` hard-codes $R_\mathrm{min}=0.975$ mm (from `setparameters.c`) to compute the contact phase; it asserts that the phase predicts the first overlap step of every run.
* Part A tests only the normal contact (head-on, no spin, no gravity); the tangential force and friction are not tested by it.
* Part B: the onset angle $\theta_c$ is defined from the table ($E_\mathrm{kin} > 10^{-6}$ J, first time after `t_tilt`), so it depends on that threshold by about $\pm0.6^\circ$; each tilt rate was run once (same initial bed, no statistics over beds). `E_kin` of the table contains translation and rotation (`dynamics.c`) and `E_pot` is the elastic energy of the contact springs only, with no gravity term, so the table's `E_tot` is not the energy of the bed in the tilted gravity field. The trajectory file has no spins.
* Part C: the settled bed is one fill (one seed, one run, one frame); the error of $\phi$ comes from varying the core bounds and from four azimuthal quarters, not from independent beds. `d_p` is taken as 2.0 mm (task text; the mean diameter of the radii is 1.9988 mm). `Z` in the table counts particle-particle contacts only and over all particles, so wall and surface particles pull it down (the particles of the core have 4.8). The fit of the oscillation of the profile is descriptive ($\chi^2$/dof about 2.6, it underestimates the height of the second peak). The no-edge run is one run with the same seed.
