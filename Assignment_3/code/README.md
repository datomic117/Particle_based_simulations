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


## Units



## Source files



## How the data in `../data/` were made



### C1/C2 production (`run_c.sh`)





## Analysis (`analysis.c`)



## Helper scripts (`tools/`)

Run from inside `tools/`; they are not needed to build or run the simulation.

- `check_b6.py`: closest interacting pair of the B6 start configuration, from `../b6_start.pdb`
  (written by the B6 run; not kept in the zip).
- `check_b7.py`: drift and transient statistics of `../../data/b7_nve.csv`.
- `render_b6.py`: renders `b6_start.pdb` with OVITO's Python module into `../b6_render.png`.

## What was added to the course code



## Who did what


