# A4 common rules. Read fully first; your part prompt refers to it.

## Project
6EMA02 Particle-based Simulations, Assignment 4 (DEM silo discharge). Group A4_8 (Tomas Dadikozyan, Shree Patil). Work in `Assignment_4/` (repo root is its parent). You do ONE part per session, then stop and report; do not start the next part.
`A4.ipynb` is both the assignment text and the report. It is NOT graded on content: it is the dossier for an individual oral exam where the students defend their own code and results with the notebook on screen. So everything you write must be correct, backed by the data, and quick to read aloud.

## Environment (Windows 11; use the Bash tool)
- gcc: `export PATH=/c/msys64/ucrt64/bin:$PATH` (the default PATH entry is broken).
- Python: `PY=/c/Users/20231118/Particle_based_simulation/Assignment_1/pbs-env/Scripts/python.exe` (3.13: numpy, scipy, pandas, matplotlib, nbformat, nbconvert). Set `PYTHONIOENCODING=utf-8`; use the Agg backend in scripts.
- The `ovito` Python module is needed for snapshots (install in Part C if missing: `$PY -m pip install ovito`).
- Run from `Assignment_4/`. Runs longer than ~1 min: `run_in_background`, do not poll, work on something else meanwhile. The machine has many cores: independent runs may go in parallel (one directory each).

## Hard rules
1. **Never edit a task cell** (metadata tag `task`; their text is checked against the release). Edit only `answer` cells; new cells you add get `metadata.tags = ["answer"]`. Edit the notebook with Python `nbformat` (read, change cells by index, write), never by hand-editing JSON. Cell map (after the restructure of `06_restructure.md`; find cells by their block tag, the indices are those of that restructure): 2 Group, 3 Results at a glance (`glance`), 8 Code map, 9 Sources, 10 Build & run log (wherever a rule below says "cell 9" for the log, it is cell 10 now); Part A: 11 task, 12 to 23 answer cells (`partA` setup, `A1`, `A2`, `Adt`, `Aphase`, `Awall`, `partA` Oral prep); B: 24 task, 25 to 34 (`partB`, `B0` to `B3`); C: 35 task, 36 to 59 (`partC`, `C1` to `C4`); D: 60 task, 61 to 90 (`partD`, `D1` to `D4`); 91 task and 92 Lessons learned (answer); 93 Hand-in. Each part is now blocks of text with its figures directly below (hidden code cells tagged `hide-input`, at most 8 printed lines each), the Oral prep is the last cell of each part, new answer cells get `["answer", "<block tag>"]`; where this contradicts "one text cell and one code cell per part" or "print the numbers quoted in the text" below, `06_restructure.md` wins. Group cell 2 is already filled; leave it.
   After every notebook change run `$PY agent_prompts/verify_tasks.py` (must print OK).
2. **Edit only what the part needs.** Code: only the stub functions of your part (keep the `\todo` comment, replace the `(void)` placeholder/`return false`, add a brief comment on how it is done; match the surrounding style and Doxygen `///` comments). Do NOT change macro defaults, constants, physics, `code/html/`, `Doxyfile`, other assignments, or the repo root. No binaries in `code/`; no second `main` in the top level of `code/` (the hand-in must build with `gcc -O3 *.c -o dem -lm`).
3. **Deliverables are not negotiable.** Re-read your task cell(s) from the notebook JSON yourself and tick every "Deliver" bullet and every "In the oral you should be able to" bullet. If a deliverable cannot be produced as specified (e.g. OVITO fails), STOP and tell the user; never substitute (no matplotlib instead of OVITO, no skipped run, no reduced table).
4. **Runs.** One directory per run: `runs/<name>/` (gitignored). Build and run, from `Assignment_4/`:
   `mkdir -p runs/<name> && gcc -O3 -Wall -Wextra <-D flags> code/*.c -o runs/<name>/dem.exe -lm && (cd runs/<name> && time ./dem.exe > run.dat)`
   Record the exact command and the measured wall-clock time of every run you keep. Trajectories stay in `runs/`; never copy them to `data/`. Anything needing a trajectory is reduced by a script in `code/` to a small file in `data/` (e.g. one frame via `np.savetxt`).
5. **Data and notebook code cells.** Hand-in files go in `data/` with prefix of the part (`A_`, `B_`, `C_`, `D_`). Notebook code cells only load from `data/` with relative paths (no absolute paths, never start a simulation), make the figure, and print the numbers quoted in the text. Figures: titled, labelled axes with units, readable fonts, legend where needed.
6. **Honesty.** Never tune parameters, windows or definitions to hit a literature value. If a result surprises you, check the code first, then report it plainly in the notebook and in the changelog under `!`. State every definition you use (onset, window, core region) in the notebook. Numbers in text come from printed output, not from memory; re-check them after the last edit.
7. **No git commands** (no add/commit/push/checkout). Do not run the checker's zip step (that is the finalize prompt).
8. **Ask only when blocked** by a decision that is the user's to make. Otherwise state your short plan (max 10 lines) and proceed.

## Answer style (the notebook is read aloud in the oral)
- Concise bullets with a bold lead-in, short sentences, plain words, no filler, no repetition of the task text. One sub-task = one short block headed `**A1 …**`.
- Order inside a Results cell: key numbers and figures first (with units and uncertainties), then the explanation (a few lines), then **`**Oral prep**`**: one `Q:` / `A:` pair per "In the oral you should be able to" bullet, answer in 1 to 3 lines, in the students' voice ("we", "our code"), citing the function or data file.
- Inline LaTeX with `$...$`. Link code with relative links (`[cylinder_wall](code/walls.c)`).
- Every statement must be something that can be defended: if you are not sure of a physical claim, test it or word it as what the data show.

## Self-check before you finish (all must pass; put the result in the changelog)
- `$PY agent_prompts/verify_tasks.py` -> OK.
- `gcc -O3 -Wall -Wextra code/*.c -o runs/_check/x.exe -lm` (mkdir first) compiles with no new warnings from your code.
- The notebook runs top to bottom in a fresh kernel: `$PY -m jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=900 A4.ipynb` (working dir `Assignment_4/`), then run verify_tasks again. No `…` placeholders left in the cells of your part.
- `git status` shows only intended files (look, do not fix with git).

## Bookkeeping (do last; the team reads this at the end)
a) **Build & run log (cell 9):** append, in the A3 style (machine line once: Windows 11, 20 cores, gcc 14.2 MSYS2 ucrt64), one row per production run or script whose output is in `data/`: exact compile+run command, simulated time, measured wall-clock, output file. Replace the cell's `- …`.
b) **Code map (cell 7):** replace the `…` of your rows with clickable relative links (`[fn](code/file.c)`, `[file](data/file)`); also list helper scripts. Every link must resolve.
c) **`code/README.md`:** Part A creates it (layout of `Assignment_3/code/README.md`: Build and run, Source files, What was added to the course code, How the data were made, Known limitations); later parts append their rows. Mention the starter `README.txt` stays.
d) **`CHANGELOG_A4.md`** (in `Assignment_4/`): append your block in the template format at the top of that file; if you used AI knowledge/literature that belongs in Sources, add it to the `Checks` line as `sources: ...`. Keep it to about 15 lines: the team reads it at the end.
