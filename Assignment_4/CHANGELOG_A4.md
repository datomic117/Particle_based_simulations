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

## Part C: Building the silo  (2026-10-20, Claude Code session)
- C: code/walls.c:cylinder_wall (C1), silo_bottom_wall edge branch (C2: closest point = rim point in the particle's azimuth), **count_in_silo (Task D2, done here on purpose)**: the table of the D7 run is Part D's data and needs N_silo; Part D only explains it. code/setparameters.c:update_parameters (RUN_SILO: orifice_open = 1 for time >= t_open). code/initialise.c:initialise_positions_silo (C3). New: code/tests/wall_check.c, code/render_ovito.py, code/run_c_noedge.sh.
- R: C_smoke | -DRUN=3 -DT_DISCHARGE=0.0 | 27 s | not kept ; D7 | -DRUN=3 | 3 min 54 s (other jobs running) | -> data/D7.dat ; C_noedge | -DRUN=3, edge branch off in a scratch copy | 4 min 45 s | -> data/C_noedge_D7.dat ; wall_check | 30 s | -> data/C_wall_check.txt ; render_ovito 6 s (ovito 3.17.0 pip-installed into pbs-env; requirements.txt not touched)
- D: data/D7.dat, C_frame_before_open.dat (frame 13, t = 0.2846 s), C_snapshot_settled.png, C_wall_check.txt, C_noedge_D7.dat; code/README.md rows added
- N: cells 7 (rows C1, C2, C3, D2), 9, 17, 18 (executed: sketch, settling figure, OVITO image, profile figure)
- Result: phi_core = 0.597 +- 0.011 (std over 48 core-bound variants), rho_b = 1493 +- 27 kg/m^3; Z = 4.51 (table; core particles 4.82); E_kin at t_open = 3.9e-9 J (5.6e-7 of max, rms speed <= 0.32 mm/s); profile: wall layer 4.13 at 0.5 d_p, oscillation period 0.86 +- 0.01 d_p, decay length 1.3 d_p, gone in the noise at ~3.5 d_p; wall_check: ~8e5 random tests, 0 failures, r_iw within 1e-7 R of brute force.
- !: (1) Without the edge (scratch run, one seed) the discharge is only 17 % faster (3292 vs 2804 /s, window 10-60 % gone), not the 47-58 % of the naive "hole wider by d_p" Beverloo estimate; the notebook says so and does not explain it. (2) Oscillation period is 0.86 d_p, the task text says "~1 d_p"; reported as measured. (3) phi of the specified core (0.597) is the top of the range of core-bound variants (0.565-0.598; wall margin 1.5 d_p gives 0.570), so the error is one-sided in practice. (4) Fill is at rest on the scale of the flow, not frozen: E_kin still falling and Z still creeping (+0.012 in the last 50 ms) at t_open; t_open not changed. (5) The frame nearest to 0.286 s is frame 13 at 0.2846 s. (6) First wall_check reference (golden section from 3600 samples) "failed" 4 checks at 2e-6 R for centres a few 1e-4 R from the axis: the azimuth is not resolvable there in double precision (shown numerically); reference changed to compare the vector only for r_c >= 0.1 R_orifice and, everywhere, that the code's point is on the rim and not farther than any sampled point. (7) OVITO's legend title is not drawn in this layout; a separate text label gives the unit.
- Checks: verify_tasks OK | compile -Wall -Wextra: 0 warnings | notebook executes: yes (python -m nbconvert, 8 s) | sources: none (closest-point argument and the layer spacing sqrt(2/3) d_p of close-packed layers are textbook; the latter is only mentioned, not tested)
