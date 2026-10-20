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
| `run_c_seeds.sh` | Part C follow-up: the silo run for the seeds 1 to 4 (`SEEDS="1 2 3 4"`), each with the edge on and off, 8 runs in parallel; every run in `runs/C_seed<s>_<edge or noedge>/` with a scratch copy of the sources whose `constants.h` has the new `SEED` (and, for `noedge`, `walls.c` with the edge off); `code/` is not changed |
| `analyse_c.py` | Part C follow-up reductions: `settling` (core packing and bed height per frame before $t_\mathrm{open}$), `crossing` (radius at which each grain crosses the plane of the plate), `seeds` (the rate $Q$ of the ten edge-on/off runs, with the definition of $Q$ of Part D2) |
| `run_d.sh` | Part D: the four extra silo runs D6, D8, D9 (Beverloo) and D3 (small orifice), in parallel, one directory each |
| `collect_times.sh` | Part D3: writes the wall-clock times and step counts of the five silo runs to `../data/D_walltimes.csv` (cost estimate) |
| `analyse_d.py` | Part D: reduces a silo `trajectories.xyz` to small files: exit time of every grain (`../data/D_exit_D7.csv`), mean speed map in (r, z) (`../data/D_speedmap_D7.csv`), or with `--stats` the per-frame speeds of the grains above the plate (`../data/D_speedstats_D3.csv`) |

## What was added to the course code

* **Part A:** `initialise_collision` in `initialise.c`: equal radii $R_\mathrm{min}$ and the matching mass, zero spin; one particle at height $1.05R$ moving down (A1), or two particles on a line along $x$ in the box middle, centres $2\cdot1.05R$ apart, moving towards each other (A2), each with `v_impact`.
* **Part B:** the `RUN_SLOPE` branch of `update_parameters` in `setparameters.c`: the slope angle is zero until `t_tilt`, then `tilt_rate * (time - t_tilt)`; it is stored in `slope_angle` (rad) and gravity is set to $g(\sin\theta, 0, -\cos\theta)$. The walls do not move. `extract_frame.py` (see above).

* **Part C:** `cylinder_wall` (C1) and the edge of the hole in `silo_bottom_wall` (C2) in `walls.c`: the closest wall point is found in closed form, the rim circle for the hole, the cylinder surface at the height of the particle (derivation in the notebook, check in `tests/wall_check.c`). `update_parameters` (RUN_SILO branch) sets `orifice_open = 1` from `t_open` on. `initialise_positions_silo` (C3) in `initialise.c`: loose square lattice of spacing 2.2 R_max inside R_cyl - 1.1 R_max, layer by layer from z_silo + 1.1 R_max, shifts of at most 0.04 R_max per coordinate. `count_in_silo` (Task D2, one loop over the particles with r.z > z_silo) was also written in Part C, because the table of the silo run is the data of Part D and must contain `N_silo`. `render_ovito.py` (see above).

* **Part C follow-up:** no change to the C code. Added `run_c_seeds.sh` and `analyse_c.py` (see the table above); the margin scans of $\phi$ and the comparison of the five seeds are in the notebook cell, so that the numbers in the text can be read from the cell output.
* **Part D:** no change to the C code (`count_in_silo` was written in Part C). Added only scripts: `run_d.sh`, `collect_times.sh`, `analyse_d.py` (see the table above). The fits (rate, Beverloo, cost) are in the notebook cell, not in a script, so that the numbers in the text can be read from the cell output.

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
| `C_seed1_edge.dat` ... `C_seed4_noedge.dat` | `bash code/run_c_seeds.sh`, then `cp runs/C_seed<s>_edge/run.dat data/C_seed<s>_edge.dat` (and `noedge`): the tables of the silo run (`-DRUN=3`, defaults) for the seeds 1 to 4, edge on and off. Seed 12345 is `D7.dat` (edge on) and `C_noedge_D7.dat` (off) |
| `C_edge_seeds.csv` | `python code/analyse_c.py seeds data/C_edge_seeds.csv`: per seed and edge setting the rate $Q$ (straight-line fit of $N(t)$ from 10 % to 60 % gone, the definition of D2) with its error (scatter of 4 pieces of the window), the time of 10 % gone, the fraction gone and $N$ at the end, $Z$ and $E_\mathrm{kin}$ at $t_\mathrm{open}$ |
| `C_settling_frames.csv` | `python code/analyse_c.py settling runs/D7/trajectories.xyz data/D7.dat data/C_settling_frames.csv`: for the frames from 0.15 s to $t_\mathrm{open}$ the core $\phi$ (core of C4, taken anew in each frame), $\rho_b$, $Z$ of the core particles, the highest centre and the mean height |
| `C_crossing_edge.csv`, `C_crossing_noedge.csv` | `python code/analyse_c.py crossing runs/D7/trajectories.xyz data/D7.dat data/C_crossing_edge.csv` (and `runs/C_noedge/trajectories.xyz` for `noedge`): per grain that leaves, the time of the first frame below the plate and the radius at which it crossed the plane of the plate, interpolated between the frames (0.022 s apart: an estimate good to about a radius) |
| `D6.dat`, `D8.dat`, `D9.dat`, `D3.dat` | `bash code/run_d.sh`: `-DRUN=3` with `-DD_ORIFICE=6 -DT_DISCHARGE=1.5`, `=8 ... =1.0`, `=9 ... =0.8`, `-DD_ORIFICE=3 -DT_DISCHARGE=0.8`, each in `runs/<name>/`; then `cp runs/<name>/run.dat data/<name>.dat`. Same columns as `D7.dat` (`N_silo` = $N(t)$) |
| `D_walltimes.csv` | `bash code/collect_times.sh`: run, $D/d_p$, wall-clock (s), simulated time (s), steps of D3, D6, D7, D8, D9 |
| `D_D7_early.png`, `D_D7_middle.png`, `D_D7_late.png` | `python code/render_ovito.py runs/D7/trajectories.xyz data/D_D7_early.png --time 0.5824 --vmax 0.1 --zrange 25 110 --xhalf 22 --size 640 1000` (and `--time 1.3815`, `--time 2.1762`): frames at 90 %, 60 %, 30 % of the particles left, 4 mm slab through the axis, speed scale 0 to 0.1 m/s |
| `D_exit_D7.csv` | `python code/analyse_d.py runs/D7/trajectories.xyz D7 --z-silo 0.0530997`: per grain (line = particle index) the position in the frame nearest $t_\mathrm{open}$ and the time of the first frame with $z < z_\mathrm{silo}$ (nan if it never leaves) |
| `D_speedmap_D7.csv` | same command: mean speed and mean $v_z$ of the grains in bins of 1 mm in $(r, z - z_\mathrm{silo})$, averaged over the frames of 0.58 to 1.92 s (the fit window of D2) |
| `D_speedstats_D3.csv` | `python code/analyse_d.py runs/D3/trajectories.xyz D3 --z-silo 0.0530997 --stats`: per frame, the median, 99th percentile and maximum speed of the grains above the plate and the position of the fastest one |
| `D_D3_blocked.png` | `python code/render_ovito.py runs/D3/trajectories.xyz data/D_D3_blocked.png --time 0.70 --vmax 0.001 --zrange 45 77 --xhalf 16 --size 900 900 --slab 4`: the blocked orifice of the $D = 3\,d_p$ run at frame 32 ($t = 0.7006$ s), speed scale 0 to 1 mm/s |

## Known limitations

* The Part A analysis uses the energy table at the grid times; $E_\mathrm{kin,before}$ needs at least one printed step before the first contact (checked by an assertion in `analyse_a.py`).
* `analyse_a.py` hard-codes $R_\mathrm{min}=0.975$ mm (from `setparameters.c`) to compute the contact phase; it asserts that the phase predicts the first overlap step of every run.
* Part A tests only the normal contact (head-on, no spin, no gravity); the tangential force and friction are not tested by it.
* Part B: the onset angle $\theta_c$ is defined from the table ($E_\mathrm{kin} > 10^{-6}$ J, first time after `t_tilt`), so it depends on that threshold by about $\pm0.6^\circ$; each tilt rate was run once (same initial bed, no statistics over beds). `E_kin` of the table contains translation and rotation (`dynamics.c`) and `E_pot` is the elastic energy of the contact springs only, with no gravity term, so the table's `E_tot` is not the energy of the bed in the tilted gravity field. The trajectory file has no spins.
* Part C: the settled bed of $\phi$, $Z$ and the profile is one fill (seed 12345, one run, one frame); the statistical error of $\phi$ comes from four azimuthal quarters of that bed, the systematic one from the scan of the wall margin, not from independent beds. The other four seeds are used only for the discharge rate with and without the edge. `d_p` is taken as 2.0 mm (task text; the mean diameter of the radii is 1.9988 mm). `Z` in the table counts particle-particle contacts only and over all particles, so wall and surface particles pull it down (the particles of the core have 4.8). The fit of the oscillation of the profile is descriptive ($\chi^2$/dof about 2.6, it underestimates the height of the second peak). The edge-on and edge-off runs are compared for five seeds (one run each); the radius at which the grains cross the plate plane (`C_crossing_*.csv`) comes from frames 0.022 s apart and is used for seed 12345 only.
* Part D: each orifice size was run once (one fill, one seed). The error of $Q$ is the scatter of 4 equal pieces of the fit window (3 degrees of freedom), so the $\chi^2$ of the Beverloo fit with 2 degrees of freedom is only a rough guide. The windows of D6 and D8 end at 34 % and 57 % gone (the runs end before 60 %). $\rho_b$ is the density of the settled core (C4), not of the flowing material. The blocked orifice of D4 is one realisation: no blocking probability can be given. The wall-clock times were measured with the four runs in parallel, so the cost per particle-step varies by a factor 1.5 between runs.
