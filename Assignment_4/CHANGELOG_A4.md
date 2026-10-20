# A4 changelog (append-only; one block per part; keep it to ~15 lines)
Legend: C = code, R = run, D = data/figure, N = notebook cell, ! = open issue / surprise

Template for each block:

```
## Part X: <title>  (<date>, <who: Claude Code session>)
- C: <file>:<function> — <what, 1 line>
- R: <name> | <flags> | <wall-clock> | -> <data file>
- D: <data/ files created>
- N: <cell indices edited>
- Result: <key numbers, 1 line>
- !: <surprises / things the team must double-check; "none">
- Checks: verify_tasks OK | compile -Wall -Wextra: <n> warnings | notebook executes: yes/no
```

---
(blocks below, newest last)

## Part A: Verification, restitution  (2026-10-20, Claude Code session)
- C: code/initialise.c:initialise_collision — equal radii R_min, mass, zero spin; A1 particle at 1.05 R above the wall moving down, A2 two particles 2.1 R apart along x moving together (v_impact each). Nothing else in code/ changed.
- R: run_a.sh (24 runs, 43 s incl. compilation, ~0.04 s per run) | -DRUN=1 -DNUM_TEST={1,2} -DV_IMPACT={0.1,0.3,1.0} -DDT_FRACTION={0.2,0.1,0.05,0.025} | -> A_restitution.csv; run_a_scan.sh (248 runs, 16 parallel, 50 s) -> A_phase_scan.csv
- D: data/A_restitution.csv, A_phase_scan.csv, A_A1_v0.3.dat, A_A2_v0.3.dat; code/analyse_a.py, run_a.sh, run_a_scan.sh, README.md (new)
- N: cells 7 (row A), 9, 11, 12 (executed, 2 figures + tables)
- Result: t_c/20: e = 0.848/0.849/0.848 (A1, tab. 0.85) and 0.919/0.920/0.919 (A2, tab. 0.92) at 0.1/0.3/1.0 m/s; within 1 % at all three speeds from 10 steps/contact (A1 worst 0.90 %); t_c/20 worst 0.24 %.
- !: Error of e is not smooth in dt: contact start and end are force jumps (dashpot eta*v) sampled at an arbitrary phase of the time grid; error = smooth function of phase + jump of +1.3..+1.5 % (A1, t_c/20) when the contact has one overlap step fewer. Scan over 31 speeds: at t_c/20 single collisions reach +1.47 % (A1, 19 % of speeds >1 %) and +0.77 % (A2); mean over phase ~0. 40 steps: all < 0.8 %. Not asked for by the task, added because "is t_c/20 enough" depends on it. The dashpot-at-onset numbers in the text come from a temporary debug print in a scratch copy (removed; code/ has none).
- Checks: verify_tasks OK | compile -Wall -Wextra: 0 warnings | notebook executes: yes | sources: none (damped-oscillator derivation read from the comments in setparameters.c)
