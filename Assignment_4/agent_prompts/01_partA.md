Read agent_prompts/00_COMMON.md first, then do Part A of A4 and nothing else.

# Part A: Verification, restitution

## Goal
Implement `initialise_collision` (RUN=1), run the 24 restitution runs, fill the Part A answer cells (11, 12), then the bookkeeping.

## Steps
1. Read notebook cells 3-6 and 10 (physics, model table, how to work, run table, Part A) from `A4.ipynb`; read `code/initialise.c` (`initialise_collision` and its `\todo`), `code/setparameters.c` (RUN 1 block; how spring/dashpot constants are derived from restitution `e` and contact time `t_c`; wall vs particle-particle; the meaning of `DT_FRACTION`, `NUM_TEST`, `V_IMPACT`), and the contact force in `code/forces.c`. You need this to answer the oral questions from the code, not from memory.
2. Implement `initialise_collision` exactly as its `\todo` says (equal radii R_min, mass = density*4/3*pi*R^3, zero angular velocity; NUM_TEST=1: one particle at height 1.05 R in the box middle moving down with V_IMPACT; NUM_TEST=2: two particles on a line along x in the box middle, centres 2*1.05 R apart, each moving towards the other with V_IMPACT). Check that anything the contact model needs (e.g. tangential displacements / collision list) is consistently initialised by the surrounding code.
3. Sanity test one run each (A1, A2) and look at the table (Step, Time, Z, Epot, Ekin, Etot): the particle must reach the contact, Epot rises and returns to 0, and Ekin after / Ekin before gives e ~ 0.85 (A1) and ~ 0.92 (A2) at the default dt. Anything else: debug before continuing.
4. Write `code/run_a.sh`: the 24 runs (NUM_TEST in {1,2} x V_IMPACT in {0.1,0.3,1.0} x DT_FRACTION in {0.2,0.1,0.05,0.025}), one `runs/A_<test>_v<v>_f<frac>/` each, as in the run table. Run it. Write `code/analyse_a.py` that extracts, per run, Ekin before the contact and Ekin after the contact has ended (define both precisely, e.g. last step before Epot > 0 and first step after Epot returns to 0 for good) and writes `data/A_restitution.csv` (num_test, v, dt_fraction, steps_per_contact = 1/dt_fraction, Ekin_before, Ekin_after, e = sqrt(after/before), e_tabulated, relative error). Copy the table of the 0.3 m/s, dt = t_c/20 run of A1 and A2 to `data/A_A1_v0.3.dat` and `data/A_A2_v0.3.dat` (small files) for the figure.
5. Time-step explanation: do not guess. Look at the dashpot force in the first step(s) of the contact (add a temporary debug print or a small separate test, remove it afterwards) and explain from what you see why the error in e is not a smooth function of dt (hypothesis to test: the contact starts at an arbitrary moment between grid points, so the first overlap and the dashpot force with the full approach velocity depend on the phase of the collision relative to the time grid). Verify that this explanation fits the numbers before writing it.
6. Notebook. Cell 12 (code): figure with two panels (A1, A2) of E_kin and E_pot against time for 0.3 m/s, the contact interval marked (e.g. `axvspan` where E_pot > 0); the table of measured e (rows: impact speed x steps per contact; columns A1, A2) next to the tabulated 0.85 / 0.92; print the number of steps per contact needed for all three speeds to be within 1 %, and whether t_c/20 suffices. Cell 11 (text): results A1/A2 (measured e vs tabulated, spread over speed), the time-step statement, the dashpot-at-onset explanation, the two sentences on why wall and particle-particle contacts need different spring/dashpot constants (derive from the code: effective mass m vs m/2, different e; same t_c), then **Oral prep** with the three Q/A pairs of the Part A task cell (where the constants come from and why dt is tied to t_c; why recovering e closes the verification loop; why wall and pair constants differ).

## Deliverables checklist (from the task cell, all required)
A1 and A2 measured e against tabulated | figure E_kin and E_pot vs time at 0.3 m/s with the contact marked | table of e at 0.1, 0.3, 1.0 m/s and dt/t_c = 0.2, 0.1, 0.05, 0.025 next to the tabulated value | statement: steps per contact for 1 % at all three speeds, is t_c/20 enough, why the error is not smooth | two sentences wall vs pair constants | the oral Q/A.

## Part-specific bookkeeping (do last, see COMMON)
- Create `code/README.md` (see COMMON c). `Assignment_4/.gitignore` already exists.
- Build & run log: `run_a.sh` row (24 runs, total wall-clock) and the machine line.
- Code map row "A test collision": `initialise_collision`, `run_a.sh`, `analyse_a.py`, the data files.
- Changelog block "Part A".
