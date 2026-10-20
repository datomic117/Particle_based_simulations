Read agent_prompts/00_COMMON.md first, then do Part B of A4 and nothing else. Part A is done (see CHANGELOG_A4.md and code/README.md; append to them, do not rewrite).

# Part B: Validation, granular slope stability

## Goal
Implement the RUN_SLOPE branch of `update_parameters`, run both tilt rates, extract the onset angle and the bottom-layer analysis, fill cells 14, 15, then the bookkeeping.

## Steps
1. Read notebook cell 13 (Part B) and `code/setparameters.c` (RUN 2 block, `tilt_rate` is stored in rad/s, the ramp stops at 24 degrees, `t_tilt` = 0.3 s, print every 100 steps, trajectory frame every 0.022 s). Implement the branch as the `\todo` says (theta = 0 before t_tilt, then tilt_rate*(time - t_tilt); store `slope_angle` in radians; set `g`). Leave the RUN_SILO branch untouched (Part C).
2. **Prediction first.** Before looking at any result, write the prediction into the notebook (cell 14): if the bed slides as one block over the floor, it starts to move at tan(theta) = mu of the relevant contact (the particle-wall friction coefficient 0.45 from the model table, since the bottom layer slides on the floor) -> theta ~ 24.2 degrees. Confirm which coefficient the model table gives and say why the others (0.50 particle-particle) do not enter.
3. Run `-DRUN=2` (20 deg/s) and `-DRUN=2 -DTILT_RATE=10` (names `runs/B_rate20`, `runs/B_rate10`). Save the tables as `data/B_rate20.dat`, `data/B_rate10.dat` (small).
4. Onset: define it explicitly (first time after t_tilt that E_kin exceeds 1e-6 J, about 1 cm/s mean speed for 2000 grains of 1.046e-5 kg; convert to angle = rate*(t - t_tilt)); compute for both rates; state if/how theta_c depends on the rate and by how much, with a reasonable uncertainty (e.g. sensitivity to the threshold 5e-7 ... 2e-6 J).
5. Frame at ~20 deg: write `code/extract_frame.py` (generic: `trajectories.xyz`, a time or frame index, output file; reads the xyz format given in cell 5 of the notebook; also reused by Parts C and D) and save the frame nearest to a 20 degree slope (for the 20 deg/s run t = 0.3 + 1.0 s) as `data/B_frame_20deg.dat` (x y z radius speed). Compute the mean speed of the bottom layer (centres < 2 mm above the floor; check where z = 0 is) and of all particles above it, and the ratio. State the expected ratio if the bottom layer rolls without slipping under a bed that slides over it (the layer's centres move at half the speed of the bed above) and what the measured ratio says about how the bed fails; use the data to explain why the friction coefficient then does not set theta_c. Only claim what the data support.
6. Notebook. Cell 15 (code): figure of E_kin (log) against slope angle for both rates on one plot, theta_c marked (and the 1e-6 J line); the bottom-layer vs upper-layer speeds as a small bar/figure; print theta_c for both rates, the ratio and the prediction. Cell 14 (text): prediction (written before the results), theta_c number (both rates) and rate dependence, comparison with the prediction, the bottom-layer analysis, **Oral prep** with the Part B Q/A pairs: why rotating gravity is equivalent to tilting the floor (rotated frame; periodic sides so no walls to move); the angle found, explained from the way the particles move; what static friction requires of the contact model (a tangential displacement stored per contact and the Coulomb cap) and what would happen to this test without it (friction would only act on sliding speed, so the bed would creep at any angle).

## Deliverables checklist
figure E_kin(log) vs angle, both rates, theta_c marked, onset definition stated | number theta_c and whether it depends on the tilt rate | comparison with the prediction | bottom layer vs upper layers mean speed at ~20 deg, expected rolling ratio, what the measured ratio says, why mu is not what sets theta_c | the oral Q/A.

## Part-specific bookkeeping (do last, see COMMON)
- Build & run log rows for the two runs (+ `extract_frame.py`).
- Code map row "B rotating gravity": `update_parameters` (setparameters.c), data files.
- `code/README.md` "Source files" / "What was added": `update_parameters` RUN_SLOPE branch, `extract_frame.py`.
- Changelog block "Part B".
