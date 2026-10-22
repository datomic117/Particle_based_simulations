Read agent_prompts/00_COMMON.md first, then do this notebook restructure and nothing else. Parts A-D and the Part C follow-up are done (CHANGELOG_A4.md, code/README.md; append, do not rewrite). This job changes the LAYOUT and the READABILITY of A4.ipynb. It changes no result: no new physics, no simulations, no new definitions, no tuning. Where this file contradicts 00_COMMON.md on notebook layout (one text cell + one big code cell per part; "print the numbers quoted in the text"), this file wins; everything else in COMMON (task cells read-only, honesty, no git, LF/CRLF care) still holds.

# Restructure: part by part, figures next to their text, hidden code, short prints, display formulas

## Why (the user's complaints) and decisions already taken
1. The answer text and the figures are far apart: each part has one long text cell followed by one big code cell with all figures and all printing. The reader must scroll back and forth.
2. The code cells are in the way of reading.
3. Several code cells print very long text (cell 12: 66 lines, cell 18: 97 lines, cell 21: 66 lines, lines up to 250 characters). It hurts readability and is mostly not quoted anywhere.
4. Formulas are almost all inline (Parts A, B, D have no display formula at all; C has one). Relations that matter belong on their own line.

The user has decided (do not ask again):
- **Scope:** all of Parts A to D. Cells 7 (code map) and 9 (build & run log) only get cosmetic fixes. Cell 8 (Sources), the Group cell and Lessons learned are not touched.
- **Code:** every code cell is collapsed by default in the notebook, and collapsed (expandable) in the HTML too.
- **Oral prep:** one block per part, at the end of the part, collapsed by default.
- **Numbers:** the text stays static markdown. Each code cell prints a small table or a few lines with exactly the numbers its text quotes (no live-generated text).

## Hard constraints
- Never edit a task cell (tag `task`); `$PY agent_prompts/verify_tasks.py` must print OK after every notebook write. Task cells keep their order; you may put answer cells anywhere between them.
- Every deliverable of the task cells stays (re-read the "Deliver" bullets and the "In the oral you should be able to" bullets of the task cells for Parts A to D; write the checklist "bullet -> block -> figure" into the changelog block).
- Every number in the text must still come from an output or a data file. Before you edit, copy the executed notebook to `runs/_restructure/A4_before.ipynb`. After the last edit run a number check: a small script (keep it in `runs/_restructure/`, not in code/) that lists the numeric tokens of the old answer text that no longer occur in the new text, and the other way round. Explain every difference in the changelog (they must all be deliberate).
- Edit the notebook with `nbformat` only (read, build the new cell list, write). The notebook is rewritten many times: write LaTeX-heavy markdown with the Write tool or with `chr(92)`, never via a Bash heredoc into Python (backslashes collapse, see the memory note), and after the last `python -m nbconvert --execute` normalise `A4.ipynb` back to LF (nbconvert writes CRLF here) and scan the answer cells for control characters (`ord(c) < 32` except newline).
- New cells get `metadata.tags = ["answer", "<block tag>"]`, where the block tag names the block (`partA`, `A1`, `A2`, `partB`, `B0`, `B1`, ..., `C1`..`C4`, `D1`..`D4`, `glance`). Keep the existing answer cells' tag `answer`.
- Do not change `code/`, `data/` contents, the run scripts, or `nb2html.py`. Cell indices will shift: update the "Cell map" line of `agent_prompts/00_COMMON.md` and the cell references of `agent_prompts/05_finalize.md` (use the block tags to find cells, and write the new map in the changelog). Do not touch other prompt files.

## Target layout

Order of cells (task cells unchanged, shown as T):

```
T  title, assessment, physics, model, code, runs          (cells 0-6, unchanged)
   Group (cell 2), Code map (7), Sources (8), Build & run log (9)   (kept; cosmetic fixes only)
   NEW  "Results at a glance"  [tag glance]  -- one table, one row per sub-task A1 ... D4:
        headline result with unit and error, and a link to its block (#anchor). At most 15 rows.
T  Part A task
   [hidden code, no output] part setup: imports, rcParams (figure style), data loading, helper functions
   markdown  ### A1 ...   text, formulas, key-number table, discussion
   [hidden code]          the figure(s) of A1 (+ at most 8 printed lines)
   markdown  caption + what to see (optional, 1 to 3 lines)
   ... same for A2
   markdown  Oral prep (collapsed)
T  Part B task   ... B0, B1, B2, B3 the same way ...
T  Part C task   ... C1, C2, C3, C4 ...
T  Part D task   ... D1, D2, D3, D4 ...
T  Lessons learned (untouched)   T  Hand-in
```

Per sub-task block, in this order (leave out what a sub-task does not have):
1. `### X1 <short title>` and one line "what we did" (2 to 4 lines at most).
2. The governing relations as display formulas (see below), each followed by a one-line "where ..." if symbols are new.
3. **Key numbers** as a small markdown table or 3 to 6 bullets (value, unit, error). No repetition of a number later in the same block.
4. The figure cell (hidden code), then a one-line caption `**Figure X1.** ...` that says what to look at.
5. The explanation: a few short bullets, hypotheses marked in italics as `*Hypothesis, untested:*`.
6. Pointers to code and data as relative links, as a single line `Files:` under the heading (not a long paragraph at the top of the part).

Placement rule: a figure sits directly under the sentence that discusses it, never more than one screen away. A code cell makes one figure (two only if they make one argument, as the two panels of one comparison). If a figure serves two sub-tasks, put it under the first and refer back with a link.

Each part's text is rearranged from the existing answer cells (11, 14, 17, 20): same content, same numbers, split by sub-task (A1/A2, B0 to B3, C1 to C4, D1 to D4). Light rewording to make the blocks self-contained is fine; new claims are not. Keep the existing Oral prep Q/A pairs word for word unless a number or cross-reference changes. Format:

```
<details>
<summary><b>Oral prep</b> (click to open)</summary>

- **Q:** ...
  **A:** ...

</details>
```
(blank lines inside `<details>` are needed for markdown and math to render; test it in the HTML.)

## Hidden code
- Every code cell: `metadata = {"jupyter": {"source_hidden": True}, "tags": ["answer", "hide-input", "<block tag>"]}`. JupyterLab and VS Code collapse `source_hidden` cells; the outputs stay visible.
- The HTML export does NOT honour `source_hidden` (tested: nbconvert 7.17 still prints the source). Use a custom template that wraps hidden cells in a closed `<details>`. Tested recipe (works with the installed nbconvert):
  - `html_template/collapse/conf.json`: `{"base_template": "lab", "mimetypes": {"text/html": true}}`
  - `html_template/collapse/index.html.j2`:
    ```
    {%- extends 'lab/index.html.j2' -%}
    {% block input_group -%}
    {%- if 'hide-input' in cell.metadata.get('tags', []) or cell.metadata.get('jupyter', {}).get('source_hidden') -%}
    <details class="hidecode"><summary>Show code</summary>{{ super() }}</details>
    {%- else -%}{{ super() }}{%- endif -%}
    {% endblock input_group %}
    ```
  - Command: `python -m nbconvert --to html --embed-images --template collapse --TemplateExporter.extra_template_basedirs=html_template A4.ipynb`
- Do not edit `nb2html.py` (course file). Write a new `make_html.py` in `Assignment_4/` that imports its functions (`strip_mathjax`, `rewrite_math`, `MATH_CSS`, `convert_text`, as in its `export`), runs the command above instead of its own, and then does what `export` does (MathML, CSS). Keep `html_template/` and `make_html.py` outside `code/` (they are not part of the zip). Add a few lines of CSS for `details.hidecode > summary` (small grey text, pointer cursor) and for the Oral prep summary.
- Checks of the HTML: no "could not convert" lines from the MathML step; `<math` present; `base64` figures present; the number of `<details class="hidecode">` equals the number of code cells; Oral prep blocks are closed and their formulas converted.
- The code is still the examinee's: say in the changelog that the code is one click away in both formats.

## Printing: small and exactly what is quoted
- Budget per code cell: at most 8 printed lines, each at most 100 characters, and nothing printed that the text does not quote. If more numbers are quoted, print one compact table instead of lines.
- Use a tiny helper in each part's setup cell that shows a markdown table (e.g. `IPython.display.Markdown` of a pipe table, LaTeX allowed in headers); no raw `print(df)` dumps, no 200-character lines.
- Typical cuts (look at each cell and decide; keep the computation, drop the print):
  - Part A (cell 12, 66 lines): print one table of e per case and the 2 or 3 numbers of the explanation; the rest stays in `data/A_*.csv`.
  - Part C (cell 18, 97 lines): `data/C_wall_check.txt` is not dumped; print its verdict line plus 3 numbers (tests, failures, max deviation) and link the file. The scans of phi against the margins go into the figure (and the 3 numbers quoted), the per-frame settling list becomes a figure panel or a 4-row table, the profile peaks/minima become markers on the figure plus one line, the five-seed table becomes one small table.
  - Part D (cell 21, 66 lines): one table for D, Q, error, one for the fit, the rest into figure annotations.
- Put numbers that a figure is about into the figure itself (title, legend, a text box): then the text can refer to the figure and the printed table can be short.
- Every number in the text appears either in a printed table of the same block, in a figure, or in a linked data file (say which, in the "Files:" line). Nothing is quoted from memory.

## Figures
- One style, set once per part in the setup cell: `figsize` about (8, 4.5) for a single panel and at most (12, 4.5) for two; font size 10 to 11; grid alpha 0.3; the same colour for the same thing across the notebook (edge on = one colour, edge off = another; use the colour-blind-safe palette the earlier parts already use).
- Titles, axis labels with units, legends where there is more than one curve, no overlapping labels (look at every figure after execution: extract the PNGs from the executed notebook and view them).
- OVITO PNGs (`data/*.png`) are shown with `IPython.display.Image(filename=..., width=...)` in their own hidden cell with a caption markdown cell.
- Do not enlarge the notebook: the executed notebook is 3.4 MB; keep dpi at the matplotlib default (100) and do not add figures that no sub-task needs.

## Formulas: display where they matter
Use `$$ ... $$` on separate lines for every relation the oral will ask about; keep inline `$...$` for symbols, short expressions and numbers with units. After each display formula a one-line "where ..." for symbols that are new. Check that `nb2html`'s converter (latex2mathml) can convert each one: avoid `\begin{aligned}`, `\tag`, `\label`, `\text` with odd characters; use one `$$` per line instead of an `aligned` block; the MathML step must report 0 failures. Candidates (confirm against the text of each part; do not invent new ones):
- **A:** restitution coefficient from the energies, contact time and time step (`t_c`, `\Delta t = t_c/20`), the pair and the wall reduced-mass relations if they are quoted.
- **B:** the sliding criterion $\tan\theta_c = \mu_w$, $\theta = \dot\theta\,(t - t_\mathrm{tilt})$, $\mathbf g = g(\sin\theta, 0, -\cos\theta)$, the rolling speed ratio argument ($v_\mathrm{top} = 2 v_c$).
- **C:** the wall contract ($\mathbf r_{iw} = \mathbf r - \mathbf r_w$, overlap $R - |\mathbf r_{iw}|$), the cylinder $\mathbf r_{iw}$, the edge $\mathbf r_{iw}$ and $r_{iw}^2$ (exists), the triangle-inequality step, $\phi = \sum_i V_i / V_\mathrm{core}$, $\rho_b = \phi\rho_p$, $Z = 2N_c/N$, the profile definition $n(r)V_p$ with the shell volume, the damped cosine of the fit, the ratio $Q_\mathrm{no\,edge}/Q_\mathrm{edge}$ with its mean and spread.
- **D:** $Q = -\mathrm{d}N/\mathrm{d}t$ and the fit window, the Beverloo law $Q = C\rho_b\sqrt g\,(D - k d_p)^{5/2}$, its linearisation $Q^{2/5} = s\,(D - k d_p)$ with $C = s^{5/2}/(\rho_b\sqrt g)$ and $k$ from the intercept, the error propagation $\sigma = 0.4\,Q^{-0.6}\sigma_Q$, the dimensional derivation (the three exponent equations on separate lines), the 2D exponent, the water-tank (Torricelli) comparison, the cost scaling.
Example of the form wanted:
```
The closest plate point to a centre over the hole is the rim point in its own azimuth, so
$$\mathbf r_{iw} = \big((r_c - R_\mathrm{orifice})\,\mathbf e_r,\ z - z_\mathrm{silo}\big), \qquad r_{iw}^2 = (R_\mathrm{orifice}-r_c)^2 + (z - z_\mathrm{silo})^2 .$$
```

## Smaller improvements to include (the user invited suggestions)
- **Results at a glance** (cell tagged `glance`, after the Group cell): one table of headline results with links to the blocks. The oral starts there.
- **Captions** under every figure ("what to see" in one line).
- **One convention for uncertainty and for hypotheses:** values as `x ± dx (stat)` / `(syst)`; every untested explanation starts with `*Hypothesis, untested:*`; every tested one cites the data file.
- **Cell 9 (build & run log, 1268 words):** keep the content; only make it easier to scan (one table per part, wall-clock column kept). No new runs.
- Do not rewrite the Lessons learned cell or the Sources cell (finalize does).

## Steps
1. Read the task cells 10, 13, 16, 19 and the four answer text cells and four code cells (11/12, 14/15, 17/18, 20/21); copy the executed notebook to `runs/_restructure/A4_before.ipynb`; list for every part the sub-tasks, their figures and the numbers each text quotes (write this list into the changelog block as the checklist).
2. Build the `make_html.py` + template first and test it on the current notebook (code collapsed, formulas converted), so the end state can be checked.
3. Restructure part by part, A, B, C, D, in this order. After each part: write, `verify_tasks`, execute the notebook (`$PY -m nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=900 A4.ipynb`, NOT `python -m jupyter nbconvert`), look at every figure of the part, count the printed lines per cell (budget above), fix.
4. Add the "Results at a glance" cell last (it quotes the final numbers).
5. Final: LF normalisation, control-character scan, `verify_tasks`, `gcc -O3 -Wall -Wextra code/*.c -o runs/_check/x.exe -lm` (code is unchanged; compile only as the COMMON self-check), number check script, `python make_html.py A4.ipynb` (writes `A4.html`), the HTML checks above, `git status` shows only intended files.
6. Update the cell map in `00_COMMON.md` and the references in `05_finalize.md`.

## Deliverables checklist
Each part goes sub-task by sub-task with its figures directly below | all code cells hidden by default in the notebook and collapsed in the HTML | at most 8 short printed lines per code cell, every quoted number still traceable | display formulas for the relations listed, 0 MathML failures | Oral prep collapsed at the end of each part, Q/A unchanged | "Results at a glance" cell | all task-cell deliverables still present (checklist in the changelog) | `make_html.py` and `html_template/` | new cell map in `00_COMMON.md` and `05_finalize.md`.

## Bookkeeping
- `code/README.md`: nothing to add (no code changed); if you add helper scripts under `code/`, there should be none.
- Code map (cell 7): unchanged unless a row points to something that moved.
- Build & run log (cell 9): add one row for `make_html.py` (command, what it does).
- Changelog block "Restructure" (about 12 lines): the new cell map, the print counts before and after per code cell, the number-check result, what moved where, the HTML checks (counts), the formulas now in display form, and `!` for anything the team must look at (for example a figure that had to be split, or a number that the check flagged).
- Self-check from COMMON: verify_tasks, compile, the notebook executes top to bottom in a fresh kernel, no `…` placeholders in the cells of your parts (the Lessons learned cell is not yours), `git status`.
