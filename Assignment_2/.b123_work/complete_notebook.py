from pathlib import Path
import copy
import json
import time
import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter

root = Path(__file__).resolve().parents[1]
path = root / 'A2.ipynb'
original = nbformat.read(path, as_version=4)
nb = copy.deepcopy(original)
nb.cells[6].source = nb.cells[6].source.replace('| B1 I/O | … | … |',
    '| B1 I/O | [`main`](code/main.c) | [CSV check](data/b1_output_check.csv) |').replace('| B2 types + LJ | … | … |',
    '| B2 types + LJ | [`initialise_types`, `initialise_velocities`](code/initialise.c), [`update_velocities_half_dt`](code/dynamics.c), [`calculate_forces_nb`](code/forces.c) | [force check](data/b3_forcecheck.txt), [NVE 1 fs](data/b3_nve_1fs.csv), [NVE 0.5 fs](data/b3_nve_05fs.csv), [NVT continuation](data/b3_equil_long.csv) |')
nb.cells[7].source = '''## Sources

**AI tools used:** OpenAI Codex (B1–B3 assistance).
**Used for:** reviewing the existing type-dependent masses and LJ implementation,
analysing saved CSV and force-check output, and drafting B1–B3 explanations and plots.
**Other external sources:** course-supplied MD starter code and assignment text.
Add any other sources or earlier AI assistance used by the group.
**Parts with substantial AI-assisted content:** B1–B3 report text and
[`analyse_b123.py`](code/analyse_b123.py), also included in the notebook below.
**Verification performed in this session:** compiled the C code and executed the
analysis against the saved data, including numerical CSV consistency checks.
The group still needs to review these explanations and complete the equilibrium test.
'''
nb.cells[8].source = '''## Build & run log

- Review on 2026-09-29, Windows with MSYS2 GCC: from the notebook directory,
  `gcc -O3 -Wall -Wextra code/*.c -o md_b123_check.exe -lm` succeeded
  (about 2 s). Warnings include unfinished bonded routines and the existing
  `%lu`/`size_t` format in the XYZ writer.
- The analysis below was executed using the existing `../pbs-env` Python
  environment. All plotted simulation data already existed before this review;
  no new production trajectory was generated in this session.
- Saved NVE output spans 50 ps: 50,000 steps at 1 fs and 100,000 at 0.5 fs,
  sampled every 10 and 20 steps respectively. These step counts and intervals
  are verified from the CSV files. The exact historical build commands, wall
  times and parameter snapshots were not saved here; they must not be invented.
- For future runs: edit `code/setparameters.c`, then from `code/` run
  `gcc -O3 *.c -o md -lm` and `./md > ../data/NEW_RUN.csv`
  (PowerShell: `./md.exe | Out-File -Encoding utf8 ../data/NEW_RUN.csv`).
  Use a new output name and archive the parameter settings for each run.
  Restart and trajectory paths in the current configuration are relative to `code/`.
'''
nb.cells[14].source = r'''Diagnostic output is written by [`main`](code/main.c) as CSV on standard
output, every `num_dt_output` steps, after both velocity-Verlet half updates.
Thus $E_\mathrm{pot}$ and $E_\mathrm{kin}$ refer to the same full timestep.
For NVT, kinetic energy is recalculated after velocity scaling.

The current columns are step, internal time, potential energy, kinetic energy,
total energy, total temperature, and the two site-type temperatures. The older
[`b1_output_check.csv`](data/b1_output_check.csv) contains the first six columns.
Its 10 rows are separated by 10 steps. The analysis below verifies finite
values, $E_\mathrm{tot}=E_\mathrm{pot}+E_\mathrm{kin}$ to $10^{-8}$ internal energy
units and the temperature identity to $10^{-12}$ K, and prints the actual residuals.

For $N=2000$ flexible sites with centre-of-mass momentum removed,

$$T=\frac{2E_\mathrm{kin}}{3N-3}.$$

Energy is in $k_B\,\mathrm{K}$, distance in Å and time in
$t_0=1.096687535$ ps; the numerical temperature is in Kelvin.
The older B1 file is only a starter-model I/O test, not a 293 K pentane run.
All plots below read saved files from `data/`.
'''
nb.cells[16].source = r'''### B2: particle types, masses and pair forces

`initialise_types` assigns each block of five sites the sequence
CH3–CH2–CH2–CH2–CH3. For 2000 sites this gives 800 CH3 sites
($m=15.035$ amu) and 1200 CH2 sites ($m=14.027$ amu).
The same type selects the mass in both velocity half updates,
$\Delta\mathbf v_i=\mathbf F_i\Delta t/(2m_i)$, and in
$E_\mathrm{kin}=\sum_i m_i|\mathbf v_i|^2/2$.
Initial Gaussian velocity components have variance $T/m_i$; subtracting the
mass-weighted centre-of-mass velocity removes total momentum.

The pair parameters used by `calculate_forces_nb` are:

| Pair | $\epsilon/k_B$ (K) | $\sigma$ (Å) |
|---|---:|---:|
| CH3–CH3 | 98 | 3.75 |
| CH3–CH2 | 67.1416 | 3.85 |
| CH2–CH2 | 46 | 3.95 |

With $\mathbf r_{ij}=\mathbf r_i-\mathbf r_j$ under the minimum-image convention,

$$\mathbf F_i=\frac{24\epsilon_{ij}}{r^2}
\left[2(\sigma_{ij}/r)^{12}-(\sigma_{ij}/r)^6\right]\mathbf r_{ij},
\qquad \mathbf F_j=-\mathbf F_i.$$

The pair-specific energy at $r_c=14$ Å is subtracted from each interacting
pair's energy. This constant does not change the force inside the cut-off.
The virial accumulator adds $\mathbf r_{ij}\cdot\mathbf F_i$ per pair before
division by $3V$. The existing mass and LJ expressions agree with these formulas.

**Scope:** there are currently no bonds in `initialise_bond_connectivity`,
and the bonded force routines are unfinished. These B3 runs test an **unbonded
binary LJ system**, with pentane's 2:3 type ratio and mass density, not liquid
pentane yet. Molecular exclusions become active only after connectivity is
implemented. The total mass gives $L=42.46121$ Å at
$\rho=626\ \mathrm{kg\,m^{-3}}$, satisfying $r_c<L/2$.

### B3: force and virial checks

The saved log tests all 2000 particles (none reported below resolution):

| Check | Finite-difference step | Largest reported relative error | Accepted tolerance |
|---|---|---:|---:|
| Force vector, maximum over sites | $\delta=10^{-6}$ Å | $7.39444\times10^{-6}$ | $10^{-4}$ |
| Global virial | scaling $h=10^{-6}$ | $4.57505\times10^{-11}$ | $10^{-8}$ |

The virial is one global check, not a per-particle maximum. The checker divides
the force error norm by the larger analytical/FD force norm; its virial
normalisation is $\max(|W_\mathrm{FD}|,|U_+|)$.

**Tolerance estimate.** Double precision has $\epsilon_\mathrm{mach}\simeq
2.22\times10^{-16}$. Using the saved early NVT endpoint as the energy scale,
$|U|\simeq1.26\times10^6$, the checker's estimate
$16\epsilon_\mathrm{mach}|U|/\delta$ is about $4.48\times10^{-3}$ internal
force units. The smallest logged force norm is about 90.2, giving a relative
resolution estimate of $5.0\times10^{-5}$; $10^{-4}$ is a conservative tolerance
just above this estimate. This is an order-of-magnitude round-off estimate,
not a rigorous bound on accumulation error. For the virial,
$16\epsilon_\mathrm{mach}|U|/(h|W|)\simeq1.1\times10^{-10}$ with
$|W|\simeq4.14\times10^7$; $10^{-8}$ allows additional summation/truncation error.
The old log does not identify its restart or parameter snapshot, so it cannot
establish the provenance or equilibration of the tested state. Repeat and
record that information for the final verification.

### B3: NVE energy conservation and timestep dependence

The following analysis fits $E(t)=a+bt$ over each complete 50 ps trace and
defines relative drift as $b/|\overline E|$, per ps. It also reports the RMS
residual **after subtracting the fitted slope**, so oscillation and drift
are not confused.

| Timestep | Relative drift (ps$^{-1}$) | Detrended RMS ($k_B$ K) |
|---|---:|---:|
| 1 fs | $-4.8970\times10^{-8}$ | 4.15875 |
| 0.5 fs | $-9.5157\times10^{-9}$ | 1.01997 |

For a smooth potential, velocity-Verlet has bounded energy error of order
$\Delta t^2$, so halving the timestep should reduce its amplitude about fourfold.
The measured RMS ratio is **4.077**, consistent with this expectation.
The slope magnitude drops by **5.146**, but a fitted slope from two finite,
chaotic trajectories is not proof of a universal drift scaling. Velocity-Verlet's
characteristic error is the bounded oscillation, not secular energy loss.
Here the energy-shifted LJ potential also has a force discontinuity at the
cut-off; crossings can add errors that are not described by the smooth-potential
argument. The plots zoom in on total energy as well as displaying all three energies.

### B3: mass check and the actual problem with these runs

For each type $s$, the diagnostic is
$T_s=m_s\langle |\mathbf v|^2\rangle_s/3$.
Means over the last 25 ps are:

| Run | CH3 (K) | CH2 (K) | Set point (K) |
|---|---:|---:|---:|
| NVE, 1 fs | 312.214 | 311.271 | 293 |
| NVE, 0.5 fs | 310.532 | 310.340 | 293 |
| Longer NVT continuation, 1 fs | 292.988 | 292.667 | 293 |

The two masses produce similar temperatures; both NVE types heat together.
Halving the timestep does not remove this heating, whereas the longer NVT
continuation brings both types near the set point. Combined with the code review,
this supports correct mass handling and points to **incomplete equilibration**.
NVE has no active temperature set point: falling potential energy becomes
kinetic energy while the total stays nearly constant.

The per-type diagnostic uses $3N_s$ in its denominator. Removing total momentum
slightly changes its finite-system expectation to
$T_s=T(1-m_s/M_\mathrm{total})$, about 292.847 K for CH3 and 292.858 K for CH2
at $T=293$ K. This roughly 0.15 K correction is smaller than the observed
sampling variations. The code below reports descriptive block standard errors;
these are not equilibrium confidence intervals for a drifting trajectory.

**Equilibration has not yet been demonstrated.** The mean potential energy in
successive 10 ps blocks of the longer NVT run falls from approximately 964261
to 914478, 864777, 859820 and 852644 internal units. A thermostat keeping $T$
near 293 K does not establish structural stationarity. These saved runs give
useful non-bonded verification, but do not yet meet B3's requirement to start
the final tests from an equilibrated configuration.

To finish that requirement, continue NVT from `b3_equil_long_restart.dat`
into **new** output files until potential-energy and both type-temperature
block means show no systematic trend over several blocks. Save this state,
run `force_test = 1` from it, and start matched-duration NVE runs at 1 fs and
0.5 fs from that same restart (`is_NVT = 0`), keeping the density fixed.
Record the parameters, commands and runtime for each run; update these tables
with the new measurements. Do not just run the current settings unchanged:
they restart from the older `b3_equil_extra_restart.dat`.

An initial lattice can make forces cancel by symmetry and samples only a small
set of separations, concealing errors that appear in irregular configurations.
A thermalised configuration tests nonzero forces and diverse local environments;
the finite-difference check still cannot detect wrong masses or a potential and
force that consistently implement the same wrong model.
'''
nb.cells[17].source = (root/'code/analyse_b123.py').read_text(encoding='utf-8')
# Execute only the completed analysis; leave later assignment placeholders alone.
test_nb = nbformat.v4.new_notebook(cells=[copy.deepcopy(nb.cells[17])])
start = time.perf_counter()
NotebookClient(test_nb, timeout=120, kernel_name='python3',
               resources={'metadata': {'path': str(root)}}).execute()
nb.cells[17] = test_nb.cells[0]
for i, cell in enumerate(original.cells):
    if 'answer' not in cell.metadata.get('tags', []):
        assert nb.cells[i].source == cell.source, f'Task cell {i} changed'
nbformat.write(nb, path)
html, _ = HTMLExporter().from_notebook_node(nb)
(root/'A2.html').write_text(html, encoding='utf-8')
print(f'Updated notebook and HTML; analysis executed in {time.perf_counter()-start:.1f} s')
for out in nb.cells[17].outputs:
    if out.output_type == 'stream': print(out.text)
