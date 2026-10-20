Read agent_prompts/00_COMMON.md first, then do this Part C follow-up and nothing else. Parts A-D and the Part C block of CHANGELOG_A4.md exist; append a new block "Part C follow-up" to the changelog and to `code/README.md` (do not rewrite their earlier content). Do NOT touch Part D cells (20, 21) or Part A/B cells.

# Part C follow-up: resolve the Part C agent's own concerns

## Why
The Part C block of the changelog lists open points (`!` 1 to 7). The team wants them closed or stated precisely before packaging. Four jobs, in this order. Budget: one batch of parallel runs (about 10 min wall-clock), the rest is analysis and text. No new physics, no tuning.

## Job 1: no-edge discharge over several seeds (concern 1)
Today "without the edge the silo discharges 17 % faster" rests on ONE seed and the notebook admits it does not know why it is only a third of the naive 47-58 %.
1. Seeds: `SEED` is a plain `#define SEED 12345u` in `code/constants.h`. Do NOT edit `code/`. Per run make a scratch source copy in `runs/<name>/src/` and `sed` the seed value in its `constants.h` (check that the fill really differs: the first rows of the first trajectory frame must differ between seeds). Use 4 new seeds (e.g. 1, 2, 3, 4), each with the edge ON (the shipped code) and edge OFF (the `return false` trick of `code/run_c_noedge.sh`), all other settings default (`-DRUN=3`). Seed 12345 is already done (`data/D7.dat`, `data/C_noedge_D7.dat`; do not rerun it). That is 8 runs of about 4-5 min: start them all in parallel in the background right away via a new `code/run_c_seeds.sh` (model: `run_c_noedge.sh` and `run_d.sh`), record each wall-clock.
2. Meanwhile do Jobs 2 to 4. When the runs finish, copy each table to `data/C_seed<s>_edge.dat` / `data/C_seed<s>_noedge.dat`.
3. Analysis with exactly the D2 definition of Q (straight-line fit of N(t), window from 10 % gone to 60 % gone, same code path as the Part D cell so the numbers are comparable; if the run ends first, say so). Write `data/C_edge_seeds.csv` (seed, edge, Q in particles/s, time at 10 % gone, N at the end, Z and E_kin at t_open). Report in the notebook: a small table per seed, then the paired ratio Q_noedge/Q_edge per seed, its mean and spread over the 5 seeds, and the spread of Q_edge alone over seeds (this is also the seed-to-seed scatter of the D7 rate: say if it is compatible with the Part D error of Q).
4. Replace every "17 %" in cells 17 (text) and its Oral prep by the multi-seed result with its spread. Keep the sentence that the naive estimate gives 47-58 %; say what the data now show.
5. Optional, only if cheap (one small script, no new runs): test the hypothesis that the gap is dynamical, not geometric: from `runs/D7` and `runs/C_noedge` trajectories, the radial position r_c at which grains first cross below z_silo (frame resolution 0.022 s is coarse; say so) and what fraction cross in the outer radius R_orifice - R .. R_orifice. If it clearly supports or contradicts the idea, say what you see; if not conclusive, keep "not investigated" and do not claim a mechanism.

## Job 2: settling at t_open (concern 4)
The fill is "at rest on the scale of the flow, not frozen" (E_kin still falling, Z creeping +0.012 in the last 50 ms). Close it with numbers, without changing t_open:
1. From `runs/D7/trajectories.xyz` take the frames from t ~ 0.15 s to t_open (every 0.022 s) and compute with a script `code/analyse_c.py` (new, model: `analyse_d.py`) phi of the specified core, Z of the core particles, and the top height of the bed, per frame; write `data/C_settling_frames.csv`.
2. State in the notebook: how much phi, rho_b and the bed height change over the last 4 to 5 frames (compare with the stated error 0.011), and that the discharge starts from this state (E_kin 5.6e-7 of its maximum). If the change is larger than the stated error of phi, the error must grow: say so and update the number everywhere it is quoted in cells 17-18 (and list the consequence for Part D in the changelog, see Job 3).

## Job 3: honest error on phi (concern 3)
phi = 0.597 +- 0.011 uses the standard deviation of 48 core-bound variants (0.565-0.598); the specified core sits at the top of that range, so the error is one-sided, and the lower values come from variants whose wall margin (1.5 d_p) lies inside the layered zone.
1. Recompute from the frame: phi vs wall margin (1.5 to 4 d_p in steps), phi vs bottom margin, phi vs top margin, each varied alone, plot or print in the existing C4 cell. Identify which variants are legitimate (margin at least the measured layering range ~3.5 d_p is the clean choice; the task's 2 d_p is the specified one) and state the number as: phi (specified core) with an explicit asymmetric or systematic error and a separate statistical one, whichever the data support. No new definition replaces the task's core: the task's number stays the headline.
2. Propagate to rho_b = phi * 2500. Check how Part D uses rho_b (cell 21 and `CHANGELOG_A4.md` Part D: 1493 +- 27 kg/m^3 gives C = 0.518 +- 0.015). If your rho_b or its error changes, do NOT edit Part D; write the exact new values and the effect on C in the changelog under `!` for the user. If nothing changes, say "Part D unaffected" in the changelog.

## Job 4: wording and oral polish (concern 2, 5, 6, 7 and readability)
Re-read cells 17 and 18 as the examinee. Cell 17 is long; the notebook must be easy to read aloud.
1. Cut it by about a quarter without dropping any deliverable of the Part C task cell (C1, C2 incl. formula/why closest/how checked/without the edge, C3, C4 numbers/figure/discussion): shorter wall_check paragraph (keep the headline: ~8e5 tests, 0 failures, 1e-7 R; details stay in `data/C_wall_check.txt`), shorter profile paragraph, no repeated numbers between Key numbers and the sections. Key numbers first, explanations after, Oral prep last.
2. Fix over-claims and loose ends: the oscillation period 0.86 d_p vs the task's "about 1 d_p" (state it as measured, one sentence on why it is below 1 d_p, marked as untested if untested); "we did not investigate/test" lines become either a tested statement (Job 1 optional step) or a short, clearly marked hypothesis; the damped-cosine fit (chi2/dof 2.6, misses the second peak) is described as descriptive only; OVITO legend note belongs in a caption, not in the argument.
3. Oral prep: keep the five existing Q/A pairs (they map to the task cell's four oral bullets) and add short pairs for the likely follow-up questions: why only X % faster without the edge (Job 1), is the fill really settled (Job 2), how sure is phi / why that core (Job 3), why is the period 0.86 d_p and over how many d_p does the layering extend.
4. Check every number in the text against the cell output after your last edit.

## Code changes
None expected. If you find a real bug in `walls.c`, `initialise.c` or `setparameters.c`, DO NOT fix it silently: stop, describe it, and say which data files (`D7.dat` and everything from Parts C/D) it would invalidate. The user decides.

## Deliverables checklist
Job 1: 5-seed table (edge on/off), paired ratio with spread, text updated everywhere "17 %" appeared | Job 2: `data/C_settling_frames.csv` + statement with numbers | Job 3: phi vs margin variants, explicit error statement, rho_b consequence for Part D recorded | Job 4: cell 17 shortened, over-claims fixed, Oral prep extended, numbers re-checked | all Part C task-cell deliverables still present.

## Bookkeeping (do last, see COMMON)
- Build & run log (cell 9): add the 8 seed runs (exact command from `run_c_seeds.sh`, wall-clock, output files) and `analyse_c.py`.
- Code map (cell 7): add `run_c_seeds.sh`, `analyse_c.py`, the new data files to the C rows.
- `code/README.md`: new rows in the source table and in "How the data were made"; update the Part C "Known limitations" (the no-edge run is no longer one seed).
- Changelog block "Part C follow-up" (about 12 lines): what changed in numbers (17 % -> ?, phi error, settling), `!` for anything Part D must be told, checks line.
- Self-check from COMMON (verify_tasks, compile, `python -m nbconvert --to notebook --execute --inplace A4.ipynb` NOT `python -m jupyter nbconvert`, which fails silently here).
