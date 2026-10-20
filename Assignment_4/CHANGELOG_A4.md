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

## Part B: Validation, granular slope stability  (2026-10-20, Claude Code session)
- C: code/setparameters.c:update_parameters (RUN_SLOPE branch) — theta = 0 until t_tilt, then tilt_rate*(time - t_tilt); stores slope_angle (rad), sets g = GRAVITY*(sin, 0, -cos). RUN_SILO branch untouched. New code/extract_frame.py (frame by time or index from a trajectories.xyz; reused by C, D).
- R: B_rate20 | -DRUN=2 | 29 s (1.5 s simulated) | -> data/B_rate20.dat ; B_rate10 | -DRUN=2 -DTILT_RATE=10 | 52 s (2.7 s) | -> data/B_rate10.dat ; extract_frame.py --time 1.3 | 0.5 s | -> data/B_frame_20deg.dat (frame 59, t = 1.2917 s, 19.83 deg)
- D: data/B_rate20.dat, B_rate10.dat, B_frame_20deg.dat; code/README.md rows added
- N: cells 7 (row B), 9, 14, 15 (executed: one figure, E_kin(angle) + speed by height)
- Result: theta_c (E_kin > 1e-6 J) = 11.7 deg (20 deg/s), 12.1 deg (10 deg/s); thresholds 5e-7..2e-6 J give -0.5/+0.6 deg; rate effect +0.4 deg, same size as the threshold effect. Prediction arctan(0.45) = 24.2 deg: bed fails at half of it. Frame 19.8 deg: bottom layer (z < 2 mm, 361 particles) 0.137 m/s, above 0.240 m/s, ratio 0.57 (rolling 0.5, block sliding 1).
- !: (1) The notebook says "rolling" is inferred from the speed ratio, because the trajectory file has no spins. Extra check in a scratch copy of the code (t_end cut to 19.83 deg, restart.dat read; not in code/ or data/, not reproducible from the hand-in): bottom-layer mean v_x = 0.135 m/s, mean omega_y*R = 0.106 m/s, so the floor contact point slips at 0.03 m/s on average (22 % of v_x, spread 0.05 m/s); upper bed has mean omega_y*R ~ 0. So the bottom layer rolls with some slip, it does not slide: consistent with the notebook text, but "pure rolling" would be too strong. (2) "Packing, not mu_w, sets theta_c" is an interpretation: mu_w was not varied (one extra run pair with muf_w changed would test it; not asked). (3) E_pot in the table is the contact-spring energy only (no gravity term), E_kin includes rotation; the notebook uses only E_kin. (4) First attempt at the runs failed (no /usr/bin/time) and produced empty files; re-run, tables bit-identical on a second run. (5) `python -m jupyter nbconvert` fails silently in this env; `python -m nbconvert` works (cell 15's outputs were empty until then).
- Checks: verify_tasks OK | compile -Wall -Wextra: 0 warnings | notebook executes: yes (python -m nbconvert) | sources: none
