# 6EMA02 Oral exam – study guide (built from YOUR four reports + code)

Written for: you (Tomas), with very little time. Everything below is taken from your own notebooks, code and the files in this folder.
Numbers are the ones in your reports. Items marked **[check]** are my own reasoning, not from your reports: verify before you rely on them.

---------------------------------------------------------------------
## 0. How to spend your time

| Time you have | Do this |
|---|---|
| **~2 h** | Read §1 (exam), then for each assignment only "30-second story", "numbers to know", and the sample-question answers. Read §8 (weak spots). |
| **~half a day** | Add: open each assignment's code pointers (§A.3) and be able to say what every function does. Do the predict drills (§7). Practise the showcase aloud. |
| **a full day** | Add: V&V matrix (§6.1), methods landscape (§6.2), and a mock oral (ask me; I examine you, you answer out loud). |

**The one fact that decides pass/fail:** D1 (can you explain your OWN code) must be ≥ 2, and no dimension may be 0. Your A3 Part C/D and almost all of A4 were written with Claude (declared in Sources). "We didn't write that part" is explicitly not an answer. So the most valuable hours are spent reading those functions until you can narrate them (list in §8.2).

---------------------------------------------------------------------
## 1. The exam in one page

25 min total, 20 examined, recorded on Teams, you share your screen (notebooks + code open).

| Segment | Time | What happens | What to do |
|---|---|---|---|
| Showcase | ~3 min | You pick ONE result and explain it. Follow-up "why does it come out like that?" is not preparable. | Pick something you can explain to L3 (see §2). |
| Deep dive | ~8 min | Examiner picks 1 assignment (you don't know which). "Show me where your code computes ..., why that way?" then a sample question about your numbers. | Know code locations (§A.3 etc.). |
| V&V | ~6 min | Other assignments: "How do you know this is right? What does that test prove? What bug could slip through?" | §6.1 matrix. |
| Extension | ~3 min | Adaptive, pushes you to L4. Weak answers get a fresh question. | Think "limits of method / design a new test". |

Levels: L1 describe → L2 explain why → L3 predict & connect (quantitative, your numbers) → L4 critique & design. **Grade = ceiling you reach.**
Dimensions (0–4 each): D1 code command, D2 V&V (what a test can and cannot catch), D3 physics of your results, D4 methods/transfer (MD vs DPD vs DEM).
Grade from sum S: 6–7 → 6, 8–9 → 7, 10–12 → 8, 13–14 → 9, 15–16 → 10. Fail if any dimension 0, D1 ≤ 1, or S ≤ 5.

Always in the exam: walk through any function/plot; "predict before you reason" (commit to an answer, then justify); a bug you hit; what you'd do differently.

**Habits that raise your score:**
1. Commit to a prediction in one sentence ("it goes up, by about 4×") *before* explaining.
2. For every test, say both halves: "it proves X, it can never catch Y".
3. Quote one of your own numbers in every answer.
4. When unsure: reason from the physics aloud, say what you'd measure. Silence/"don't know" gets you a fresh question but no points.
5. Declared untested hypotheses are fine: say "that is a hypothesis, I did not test it, here is how I would".

---------------------------------------------------------------------
## 2. Showcase: pick one (60-90 s script each)

Pick a result where you can survive "why?" down to L3.

1. **A4 Beverloo / flow independent of height (recommended: clean story, lots of numbers).**
   "A silo of 7500 grains empties at a constant rate Q = 2804 ± 18 particles/s: first half of the window 2842/s, second half 2800/s (−1.5 %), while a water tank would drop 18 %. Fitting Q^(2/5) vs D gives k = 1.59 ± 0.05 (literature 1.5) and C = 0.518 (literature 0.58; 4σ low, I report that honestly)."
   Likely follow-up: why constant? (Janssen: wall friction carries the weight, pressure saturates, flow set by free-fall zone of size ~D.)
2. **A3 surface tension vs Groot–Warren.** σ_s from pressure-tensor anisotropy matches Eq. 36 within 1–2 % for Δa ≥ 14, control run at a_AB = 25 gives 0.009 ± 0.017 (zero). Follow-up: why must the sign be positive; why factor ½ L_x (two interfaces).
3. **A1 Apophis.** Miss distance 40,202 → 38,548 → 38,158 → 38,057 km for dt 2 h → 15 min, ratios 4.2, 3.9 (second order), Richardson 38,024 km vs Horizons ~38,000 km. Follow-up: why second order; what limits accuracy.
4. **A2 pentane effect.** g+g+ is 17.5 ± 1.5 times more likely than g+g−, ΔF ≈ 838 K ≈ 2.9 kT. Follow-up: which interaction (1–5 end-group LJ), and "did you test it?" (No: the ε=0 control was not run; be ready to say so).

---------------------------------------------------------------------
## A1. Planetary motion (N-body gravity, Euler vs velocity-Verlet)

### A.1 30-second story
12 bodies from a JPL Horizons snapshot (2026-09-01), SI units, Vec3D in C. Pairwise gravity with Newton's 3rd law (66 pairs). Compare Euler-forward and velocity-Verlet by conservation laws (E, L, COM), validate against the Horizons snapshot one year later, run outer planets for 548 years, predict Apophis 2029 flyby.

### A.2 Numbers to know
- 12 bodies, N(N−1)/2 = **66** pair evaluations per step.
- C3 (1 year, dt = 1, 0.5, 0.25, 0.125 day): max |ΔE/E0| Verlet 2.2e-6 → 3.5e-8, **slope 1.99**; Euler 8.7e-3 → 2.0e-3, **slope 0.71** (not yet asymptotic ~1). Euler L drift 6.7e-4 → 9.4e-5 (∝ dt); Verlet L conserved to round-off at every dt. COM velocity "drift" 1e-6–1e-5 is identical for both and not shrinking with dt: it is finite-difference noise from differentiating the stored R_com, not physics.
- C1: fitted Earth year 364.39 d (−0.24 % vs 365.256); Moon 27.357 d on 3 months (+0.13 %), 27.324 d on 1 year (+0.009 %).
- C2 (1 year vs Horizons, dt = 24 h / 6 h / 1 h): Mercury 6.4e9 / 1.8e8 / **4.7e6 m**, Earth 2.6e9 / 6.1e6 / 1.2e5, Jupiter 1.1e9 / 7.8e3 / 4.3e2, Pluto 4.5e8 / 138 / **138 (plateau)**, Sun 2.3e2 plateau. Inner bodies still dt-limited at 1 h; outer bodies hit a floor that dt cannot fix.
- C4: 200,000 d (~548 y), Sun + 5 outer planets: |ΔE/E0| ≲ 7e-8, bounded; COM moves in a straight line at ~0.16 m/s (system has non-zero total momentum); perihelion/aphelion within ~0.5 % of literature.
- C5: closest approach 2029-04-13 ~21:46 TDB. Miss distance 40,202 (2 h), 38,548 (1 h), 38,158 (30 min), 38,057 km (15 min). Difference ratios 4.23, 3.88 → O(dt²). Richardson 38,024 km, quadratic fit 38,017 km, Horizons ≈ 38,000 km. No Moon: 2.15 million km (a day earlier).

### A.3 Code pointers ("show me where you compute ...")
- Forces + potential, pair loop j = i+1, 3rd law: [forces.c:4-49](../Assignment_1/Submission/code/forces.c#L4-L49). Acceleration uses only m[j] (no division by own mass) → zero-mass Apophis needs no special case.
- Euler: [integrate.c:8](../Assignment_1/Submission/code/integrate.c#L8). Verlet: [integrate.c:28](../Assignment_1/Submission/code/integrate.c#L28): saves a_old, position update r += v dt + ½ a dt², **recompute force at new positions**, then v += ½ (a_old + a_new) dt. (Force in the middle: if you evaluate it elsewhere you use stale acceleration and lose 2nd order/symplectic structure.) Integrator switch `integrate_step` at [integrate.c:69](../Assignment_1/Submission/code/integrate.c#L69).
- main loop, CLI parsing, mallocs: [main.c](../Assignment_1/Submission/code/main.c). Vec3D helpers: [vec3d.h](../Assignment_1/Submission/code/vec3d.h). Reader: [io.c](../Assignment_1/Submission/code/io.c).

### A.4 Sample questions with answers
**(L2) Euler energy error grows steadily, Verlet stays bounded. Why?**
Euler is not symplectic. For an oscillation each step multiplies the amplitude by sqrt(1 + (ω dt)²) > 1, so the phase-space area expands and E grows secularly (error O(dt) per unit time). Velocity-Verlet is time-reversible and symplectic: it exactly conserves a "shadow Hamiltonian" H̃ = H + O(dt²), so the true E just oscillates in a band of size O(dt²) around a constant. Evidence: Verlet slope 1.99 and flat band vs Euler monotone drift (C3 figure).

**(L2) Zero-mass bodies?** a_i = Σ_j G m_j (r_j − r_i)/r³ does not contain m_i, so Apophis (m = 0) is accelerated by others and pulls on no one (its m_j = 0 contributes nothing). F/m would be 0/0 = NaN, which then propagates into v and r.

**(L3) Conservation tests were clean. What does the Horizons comparison check that they cannot?**
Conservation tests are internal consistency (verification): any central, pairwise-antisymmetric force law, even the wrong one (say 1/r^2.5, or wrong G), conserves L and COM exactly, and E if U is consistent with it. Horizons is independent truth: it validates equations, G·m values, initial data and which bodies are included (validation). Example from your data: barycentre vs planet-centre states for outer planets.

**(L3) Force magnitude bug, direction right: which diagnostics catch it?**
L and COM stay perfectly clean (they only need a central force and F_ij = −F_ji). Energy catches it **only if U is coded independently of a** (your U is a separate line, so a wrong magnitude in the acceleration makes E = K + U drift/oscillate, because a ≠ −∇U). If the same wrong factor sits in both a and U, even E is conserved and only Horizons catches it.

**(L3) What sets the largest usable dt, and what changes with a very eccentric orbit?**
The shortest dynamical time among bodies: Moon around Earth (27 d), Mercury (88 d), and the Apophis flyby (hours: needed 15 min). Rule of thumb dt ≲ T/100 and judge by an error budget. For a very eccentric body the relevant scale is the periapsis passage time ~ sqrt(r_p³/GM): much shorter than the period, so you need a much smaller dt, or an adaptive step (but a variable step breaks symplecticity unless done carefully).

### A.5 Extra likely questions
- Why is the Euler L drift O(dt)? One Euler step adds dt² Σ m_i v_i × a_i (not cancelled by Newton's 3rd law), Verlet's symmetric kick-drift-kick cancels it. (You corrected your own A4 table after C3: tell this as your "bug/learning" story.)
- Why can a symplectic integrator keep E bounded yet have phase error? Shadow Hamiltonian has slightly different orbital frequency; frequency error adds up linearly with number of orbits (C4 text), so the planet is at the wrong point on a correct ellipse.
- Why barycentre for Jupiter–Pluto? Model has no moons; the system barycentre moves smoothly under the Sun, planet centre wobbles around it. **[check]** order of magnitude: Pluto's centre orbits the Pluto–Charon barycentre at ~2100 km, speed ~24 m/s; a velocity error of ~10 m/s over a year gives ~3e8 m, vs your 138 m. Jupiter: ~1 m/s wobble → ~1e7 m. No dt can fix a wrong initial condition.
- Why does Richardson work? Error ∝ dt² so D(h) − D(2h) ≈ −3 c h² → D0 = D(h) + (D(h) − D(2h))/3.
- Why Apophis harder than planets? Close encounter amplifies errors (miss distance very sensitive to Earth's position and Apophis's velocity before the flyby); the encounter lasts hours so dt must resolve it.
- L4 test design: in 2-body Kepler check the Laplace–Runge–Lenz vector: it is conserved only for exactly 1/r² forces, so a wrong exponent shows up as perihelion precession that E, L, COM would never see.

---------------------------------------------------------------------
## A2. MD of liquid n-pentane (united atom)

### A.1 30-second story
400 pentane molecules (2000 sites, CH3/CH2 types), L = 42.46 Å at 626 kg/m³, T = 293 K, internal units (energy k_B·K, Å, amu; time unit 1.0967 ps; 1 fs = 9.118e-4). LJ (Lorentz–Berthelot, r_cut 14 Å, shifted) + bond/angle/torsion (RB). Verification: FD force tests, NVE energy conservation, per-type temperature (mass check). Berendsen thermostat. Production: conformer statistics (pentane effect) and MSD → D.

### A.2 Numbers to know
- Units: t0 = sqrt(M0 L0² / E0) = 1.0967 ps; v0 = 91.18 m/s; v_rms(CH2, 293 K) = 721.8 m/s.
- FD test (equilibrated config, δ = 1e-6 Å): force rel. error **5.7e-6**, virial 1.4e-9, tolerance 5e-5, noise floor F_noise ≈ 16 ε_mach |U|/δ = 4.3e-3 vs smallest force 136 → floor 3.2e-5.
- NVE (1 fs vs 0.5 fs): bounded oscillation ratio **3.89** (expected 4); linear drift −6.4e-8 vs −1.0e-8 per time unit (ratio 6.1, not 4: you honestly don't claim a scaling).
- Per-type temperatures T_CH3 = 292.88, T_CH2 = 292.85 K (diff 0.04 K) → masses typed correctly.
- Bonded FD errors: bond 1.0e-9, angle 1.2e-9, torsion 2.1e-8; torsion ΣF ~1e-15, torque ~1e-13.
- B6: closest pair 3.24 Å, r/σ = 0.82. B7 NVE from lattice: T falls from 293 to ~203 K (lattice melts, kinetic → potential), drift −1.15e-5/ps.
- B9 thermostat: NVE 204 K; τ = 0.01 ps → 292.99 ± 1.62 K; τ = 1 ps → 291.67 ± 3.52 K.
- C1 (4 replicas × 200 ps used): tt 0.473, tg (4 states) 0.483, g+g+ & g−g− 0.0413, g+g− & g−g+ 0.00236. g+g+/g+g− = **17.5 ± 1.5**, ΔF = 838 K = 2.86 kT. Gauche well 476 K above trans, barrier 1668 K = 5.7 kT. Residence times: trans 36 ps, gauche 15 ps, relaxation 10.4 ps.
- 1–5 distance in ideal geometry: tt 5.17, tg 4.56, g+g+ 3.76, g+g− 2.87 Å (σ33 = 3.75).
- C2: D = **(6.81 ± 0.07)e-9 m²/s** (fit window 10–50 ps), Yeh–Hummer correction +0.6e-9 → 7.4e-9; experiment ~5.6e-9 (unverified source!), so +22 % (+32 % corrected). MSD is diffusive after ~2 ps; ballistic limit 3kT/M t² checked without fit.

### A.3 Code pointers
- LJ non-bonded (types, mixing, shift, exclusion factor): [forces.c:447](../Assignment_2/code/forces.c#L447). Exclusion factors f12 = f13 = f14 = 0 in [setparameters.c:138](../Assignment_2/code/setparameters.c#L138), applied in [nbrlist.c:242](../Assignment_2/code/nbrlist.c#L242) (the 1–5 pair is a normal pair).
- Bond [forces.c:44](../Assignment_2/code/forces.c#L44), angle [forces.c:105](../Assignment_2/code/forces.c#L105), dihedral [forces.c:249](../Assignment_2/code/forces.c#L249) (works with cos φ, no arccos).
- Velocity half step with type masses [dynamics.c:38](../Assignment_2/code/dynamics.c#L38), Berendsen [dynamics.c:103](../Assignment_2/code/dynamics.c#L103): λ = sqrt(1 + dt/τ (T0/T − 1)), T = 2 E_kin/(3N − 3).
- Analysis: dihedral angle [analysis.c:13](../Assignment_2/code/analysis.c#L13), sampling [analysis.c:90](../Assignment_2/code/analysis.c#L90), MSD with unwrapped COM and many time origins [analysis.c:130](../Assignment_2/code/analysis.c#L130).

### A.4 Sample questions with answers
**(L2) Why no LJ for 1–2 and 1–3 pairs; what would you see if you forgot?**
Their geometry is already set by the bond and angle terms; counting LJ too is double counting, and LJ is wrong at such short range. Numbers: 1–2 pair r = 1.54 Å vs σ ≈ 3.85 Å: (σ/r)^12 = 2.5^12 ≈ 6e4, energy ~ 1e7 K. 1–3 pair r ≈ 2.58 Å: (σ/r)^12 ≈ 120, energy ~ 3e4 K ≈ 100 kT. **[check arithmetic]** You would see bonds stretched and angles opened (equilibrium geometry shifted), huge potential energy, wrong conformer statistics, possibly blow-up at 1 fs.

**(L3) Why must the FD force test be run from an equilibrated configuration?**
Fresh chain sits exactly at the minimum of bond/angle/torsion terms: bonded forces are ~0, so relative error = round-off/round-off (reports big error on correct code, would also pass wrong code). On a lattice, symmetry cancels forces or makes them tiny, and few pairs are close. Equilibrated state has generic non-zero forces from every term and a realistic distribution of pair distances. Also the FD floor is absolute (16 ε|U|/δ), so it hurts the smallest forces most; equilibrated config gives a safe minimum force (136) so the floor is 3e-5 → tolerance 5e-5.

**(L2) Why NVE for energy conservation, and how should the drift scale with dt?**
Berendsen is non-Hamiltonian: it adds/removes energy by design, so only NVE isolates integrator error. The bounded oscillation of E_tot is the integrator's own error and scales ∝ dt² (you measured 3.89 for halving). A systematic slope is NOT purely the integrator: sources are round-off, neighbour-list misses, and the force discontinuity at r_cut (you shift the potential, not the force), so it need not scale as dt² (you measured 6.1). Separate them (you fit a line, then take the residual).

**(L3) Raise T by 100 K: which way does the trans fraction move?**
Down (more gauche): Boltzmann weights exp(−ΔE/kT) with the gauche well ~476 K above trans (free-energy measurement). **[check]** From your 3×3 table, P(trans) per dihedral = 0.718 at 293 K (g/t = 0.39). Assuming ΔF ≈ 476 K constant: g/t at 393 K ≈ 0.39·exp(476·(1/293−1/393)) = 0.59, P(trans) ≈ 0.63. Also g+g− is penalised by ~2.9 kT at 293 K; this penalty shrinks in kT units at higher T, so g+g− grows relatively faster. D also increases.

**(L4) What does the FD test prove, which bug class can it never catch, and what covers the gap?**
Proves: coded force = −∇(coded energy), i.e. internal consistency (and virial consistency). Never catches: a wrong potential that is consistently wrong in both U and F: wrong ε/σ/mixing, wrong masses, wrong exclusions/1–5 handling, wrong units, a neighbour list that misses pairs (missing in both), wrong integrator or thermostat. Covered by: independent checks: per-type temperature (masses), N² brute-force pair count vs neighbour list, hand-computed energy of a small cluster, NVE conservation (catches dt/force-inconsistency only), and validation against physics (conformer ratio from force field, D vs experiment, density).

### A.5 Extra likely questions
- Why 3N − 3 in T? Total momentum fixed at zero removes 3 degrees of freedom. Why mass-weighted CM removal? Momentum is Σ m v.
- Why σ/ε sanity: Maxwell–Boltzmann width sqrt(kT/m) per type.
- Why does T fall to 203 K in B7? Lattice start is not an equilibrium liquid; melting converts kinetic into potential energy (E_pot rises). Your lessons-learned says: equilibrate the packing with NVT before B7–B9.
- Berendsen: what does short τ buy/cost? Fast equilibration, small T fluctuation (1.62 K) but strongly perturbs the dynamics; **does not sample the canonical ensemble** (suppresses fluctuations; kinetic energy fluctuations wrong). Why NVE for diffusion: thermostat rescaling alters dynamics.
- MSD: why unwrapped coordinates (wrapped MSD saturates ~L²/4), why 6 (= 2d for d = 3), why COM (removes vibration), why start fit at 10 ps (slope 1.1 at 2 ps; need pure diffusive), why end at 50 ps (fewer, overlapping time origins: noise), error from 4 replicas.
- Why 4 replicas × 300 ps and discard 100 ps? T_CH3 − T_CH2 was still 3.4 K in the first 100 ps of production: ends/interior equilibrate slower than conformers.
- Finite-size: D_PBC = D_∞ − 2.837 kT/(6πηL): diffusion in a periodic box is slower by hydrodynamic backflow.
- Why is D too high? United-atom is smoother (no H), fitted to thermodynamics not transport. (Your Lessons: not tested.)
- Free-energy view: −kT ln P(φ) vs torsion potential coincide (single dihedral = Boltzmann of U_tors); pentane effect only in the pair.
- Predict: set the 1–5 LJ (ε33) to zero → g+g+ and g+g− become about equally likely (not run!).

---------------------------------------------------------------------
## A3. Dissipative particle dynamics (Groot & Warren)

### A.1 30-second story
Converted your MD code to DPD: soft beads, kT = 1, m = 1, r_c = 1, ρ = 3, a_AA = a_BB = 25. Three pair forces: conservative F^C = a(1 − r)r̂, dissipative F^D = −γ w² (r̂·v_ij) r̂, random F^R = σ w ζ r̂/√Δt. Thermostat = dissipative + random (pairwise, momentum conserving); Berendsen removed. Tests with parts switched off (energy conservation, ideal gas, velocity variance), then validation vs Groot–Warren: g(r), demixing, χ(Δa), surface tension, chains.

### A.2 Numbers to know
- γ = 4.5, **σ = sqrt(2γkT) = 3**, Δt = 0.04, r_shell 0.4 (run script) / 1 in committed file, a = 25 = 75kT/ρ (water compressibility).
- B1 (conservative only, t = 8): max rel. energy deviation 1.1e-3, 1.7e-4, 5e-5, 8e-6 for dt 0.04 … 0.005, scaling **dt^2.33**.
- B2: g(r) flat = 1.00 ± 0.01–0.02 (ideal gas). B3: measured kT = 1.030–1.031 vs 1.0 (3 % too hot), for conservative on and off; Groot–Warren report ~2 % at dt 0.04 for λ = ½.
- C1: g(r) peak 1.15 at r ≈ 0.88, block scatter 0.002.
- C2/C3 (a_AB 35/39/43/47): φ_A = 0.106 / 0.0196 / 0.0066 / 0.0017; χ = 2.71 / 4.07 / 5.08 / 6.39. Slope χ/Δa through origin 0.286–0.288 vs Eq. (24) **0.286**; free fit over χ > 3: 0.268. Steady after ~20,000 steps (bulk φ_A scatter 0.003).
- C4: σ_s = (L_x/2)[P_xx − ½(P_yy + P_zz)]. Control a_AB = 25: 0.009 ± 0.017. a_AB 35, 37, 39, 43, 47: 0.23, 0.57, 0.86, 1.37, 1.79; measured/Eq. 36 = 1.08, 1.05, 1.02, 0.99, 1.00.
- D1: ⟨R²⟩ = 6.75 vs ideal 3(N−1)kT/C = 6 (N = 5, C = 2).
- D3: χN/Δa ≈ 0.84–0.97 (N=3), 1.42–1.48 (N=5), 2.3–2.8 (N=10) vs 0.306 N; fitted slope 0.266.

### A.3 Code pointers
- DPD pair forces: [forces.c:112](../Assignment_3/code/forces.c#L112): `wC = 1 − r/r_cut`, `fmag = a wC`, `Epot += 0.5 a r_cut wC²`, `wR = 1 − r`, random `zetadt = gauss()/sqrt(dt)`, dissipative `−γ wR² (r̂·v_ij)`, applied with +ft to i and −ft to j (so ζ_ij = ζ_ji by construction: one draw per pair).
- σ in [setparameters.c:43](../Assignment_3/code/setparameters.c#L43). Switches `conservative_force_on`, `dissipative_random_force_on` [setparameters.c:27](../Assignment_3/code/setparameters.c#L27).
- `gauss()` in [random.c](../Assignment_3/code/random.c) is actually **uniform on [−√3, √3]** (variance 1), as recommended by the paper; name is misleading, generator is `rand()`.
- Spring (r0 = 0, k = C = 2): [forces.c:48](../Assignment_3/code/forces.c#L48). Verlet loop (no thermostat) [main.c](../Assignment_3/code/main.c): half-kick, drift, PBC, nbr list, forces, half-kick.
- Thermostat removed: `thermostat` still exists at [dynamics.c:108](../Assignment_3/code/dynamics.c#L108) but is not called.
- Virial per axis and pressure tensor: [forces.c:23](../Assignment_3/code/forces.c#L23) and [dynamics.c:41](../Assignment_3/code/dynamics.c#L41) (kinetic part Σ m v_α²/V).
- Analysis (g(r), profiles, pressure, chain size): [analysis.c](../Assignment_3/code/analysis.c) (`sample_rdf`, `sample_profile`, `sample_pressure`, `sample_chain_size`).

### A.4 Sample questions with answers
**(L3) The 1/√Δt in the random force. Where from, what breaks without?**
Fluctuation–dissipation: ⟨F^R(t)F^R(t')⟩ = 2γkT w² δ(t − t'). On a discrete grid δ(t) → 1/Δt, so the amplitude ∝ 1/√Δt: the velocity kick per step F^R Δt = σ w ζ √Δt has variance ∝ Δt (Wiener increment), exactly like friction removes energy ∝ Δt. Without it (F^R = σwζ) the kick ∝ Δt has variance ∝ Δt²: noise vanishes as Δt → 0 relative to friction, so the temperature drops far below kT and depends on Δt (system freezes); too large amplitude would heat. Your B3 test (T from velocity variance = 1.03) is what catches it.

**(L2) Why are σ and γ not independent?**
FDT σ² = 2γ kT (with w^R² = w^D). The balance of random heating and frictional cooling gives equilibrium temperature kT = σ²/2γ. Wrong σ² by factor f → T off by f (σ off by f → T off by f²). With γ = 4.5 you derive σ = 3.

**(L3) Ideal-gas test (conservative off): what must come out, which bugs does it isolate?**
g(r) = 1 flat (no structure), Maxwellian velocities with variance kT/m = 1 (and p = ρkT = 3 if you check pressure). Only dissipative + random remain, they are not derivable from a potential (U = 0), so positions are uncorrelated. Isolates: g(r) normalisation (your first versions were not flat), neighbour list/PBC pair completeness, a conservative force that is "not really off" (would show a peak at r ≈ 0.9 and depletion at small r), and the FDT balance via T. Does NOT test: sign/magnitude of the conservative force (B1 energy conservation and C1 g(r) do), or momentum conservation per pair.

**(L2) Why did Berendsen have to go?**
In DPD the pairwise dissipative + random forces ARE the thermostat. Extra global velocity rescaling fights the FDT balance (the equilibrium T is set by σ, γ; rescaling forces a second target and suppresses fluctuations), and is non-local so it destroys the hydrodynamics that DPD's local momentum conservation gives.

**(L3) For what physics does DPD momentum conservation matter?**
Anything with flow: viscosity, shear, hydrodynamic interactions between colloids/polymers, interfacial dynamics, blood cells/micelle transport, domain growth (Navier–Stokes at large scales). Langevin/Berendsen thermostats (not momentum conserving) give Brownian dynamics without hydrodynamics. In static equilibrium properties (g(r), χ, σ_s, phase separation) both thermostats would give the same ensemble, and indeed your Part C/D results are mostly equilibrium.

### A.5 Extra likely questions
- Derive U^C: F = a(1 − r/rc) r̂ → U = (a rc/2)(1 − r/rc)², U(0) = a rc/2 = 12.5. (**Your notebook B1 cell typed (1 − r/(2 rc))²: typo; the code `0.5*a*r_cut*wC*wC` is right.**)
- Why does B1 need smaller dt than production? With thermostat on, friction/noise absorb integration errors; with only conservative forces nothing removes them. Largest step (0.04) still conserves to 0.1 %.
- Why 3 % too hot at dt = 0.04? Velocity-Verlet with dissipative force evaluated at half-step velocity (λ = ½ in G–W); integration error, expected to shrink with dt (you did not check).
- Why bonded beads also repel (f12 = f13 = f14 = 1)? Opposite of A2: DPD repulsion is soft and a chain bead is a monomer bead; exclusion would give ideal chains (⟨R²⟩ = 6 not 6.75) and break the χ ↔ Δa mapping.
- Measuring χ: Flory–Huggins, equal chemical potential of A in both phases → χ = ln[(1−φ)/φ]/(1−2φ) for N = 1; requires flat plateaus away from interfaces, steady state, symmetric coexistence (you checked plateau slope, trend, mirror).
- Surface tension sign/trend/limit: σ_s > 0 (interface costs free energy; negative → phases would mix/interface grow); increases with Δa; → 0 at critical point χ ≈ 2.36 as (1 − 2.36/χ)^(3/2).
- Why longer chains demix more easily: χ per bead depends on Δa only; mixing entropy per chain is ~kT whatever N, but contact penalty ∝ N: χN = 0.306 N Δa, critical χN ≈ 2.
- Why does a_AB − 25 map to χ (Δa)? χ ≈ 0.286 Δa at ρ = 3 (G–W Eq. 24). a_AA = 25 fixed by compressibility.
- Equilibrated ↔ discarded first 20,000 steps because B-phase bulk φ_A still rising.
- L4 design: measure total momentum (should be exactly conserved to round-off, the check that separates DPD from Langevin); viscosity from shear/Poiseuille flow; Schmidt number is ~1 (too low vs real liquids ~1000) – a known limitation of DPD; dt-dependence of T.

---------------------------------------------------------------------
## A4. Discrete element method: silo discharge and Beverloo

### A.1 30-second story
Linear spring-dashpot (normal) + Cundall–Strack tangential spring with stored displacement and Coulomb cap, spheres d_p = 2 mm, ρ_p = 2500. Parts: A restitution tests (e = 0.85 wall, 0.92 pair), B tilted bed (rotating gravity), C silo walls, fill, packing, D discharge, Beverloo, blocked small orifice.

### A.2 Numbers to know
- e_n pp 0.92, pw 0.85; e_t 0.33; μ pp 0.50, μ pw 0.45. k_n = 1000 N/m **chosen**; t_c = **0.219 ms** follows; Δt = t_c/20 = 1.095e-5 s; wall k = 2× pair k; wall η 3.9× pair η.
- A: e measured 0.848 / 0.849 / 0.848 (wall), 0.919 / 0.920 / 0.919 (pair) at 0.1/0.3/1 m/s. 10 steps/contact gives <1 % at all three speeds, 5 does not. Error not smooth in dt: force jumps on/off at full dashpot strength, so error depends on **phase of contact onset on the time grid** (single collision up to +1.5 %, mean right).
- B: onset θ_c ≈ **12°** (11.7° at 20°/s, 12.1° at 10°/s) vs block-sliding prediction arctan(μ_w) = **24.2°**; bottom layer speed ratio 0.57 (rolling 0.5, sliding 1) → bottom layer rolls.
- C: 7500 grains, R_cyl = 17.7 mm (17.7 d_p wide), bed ~57 mm; E_kin at t_open = 5.6e-7 of max; φ = **0.597 ± 0.005** (−0.013 syst), ρ_b = **1493 ± 12** kg/m³, Z = 4.51; wall layering period 0.86 d_p, range ~3.5 d_p. Edge of the hole: without it discharge **16 ± 2 %** faster (5 seeds). wall check 685,455 random particles, 0 failures.
- D: Q(7 d_p) = **2804 ± 18 /s** = 29.33 g/s; halves 2842 / 2800 (−1.5 %) vs water tank −18 %. Beverloo fit of Q^(2/5) vs D (6,7,8,9 d_p): **k = 1.59 ± 0.05**, **C = 0.518 ± 0.015** (4.2σ below 0.58), χ² 1.39/2 dof. D = 3 d_p: blocked, 8 of 7500 grains in 0.8 s (fitted law would give ~78). Cost of a 2×-wide silo: 3.4 h (×13.7).
- Check: C ρ_b √g (D − k d_p)^(5/2) at D = 7 d_p: 0.518 × 1493 × 3.13 × (10.8 mm)^(5/2) ≈ 29 g/s ✓.

### A.3 Code pointers
- Contact forces: [forces.c:60](../Assignment_4/code/forces.c#L60) (pp), [forces.c:189](../Assignment_4/code/forces.c#L189) (pw): normal spring k·overlap, dashpot −η (v_rel·n̂)n̂ scaled by mass_factor = sqrt(m_eff/m_ref), tangential spring −k_t·t_ij (**stored contact displacement tij**), tangential dashpot on surface velocity incl. spin, Coulomb check `dft.sq >= fric² dfn.sq` → rescale and reset stored displacement to −F_t/k_t (sliding).
- Constants from (e, k_n): [setparameters.c:162-175](../Assignment_4/code/setparameters.c#L162-L175): t_c = sqrt(½ m_ref (π² + ln²e)/k_n), η = −2 ln e · m_eff/t_c, k_t = (m/7)(π² + ln²e_t)/t_c² (the 2/7 comes from rotational inertia of a sphere), Δt = DT_FRACTION · t_c ([setparameters.c:185](../Assignment_4/code/setparameters.c#L185)).
- Rotating gravity: RUN_SLOPE branch [setparameters.c:203-212](../Assignment_4/code/setparameters.c#L203-L212): g = g(sinθ, 0, −cosθ).
- Walls: cylinder [walls.c:134](../Assignment_4/code/walls.c#L134), silo plate + hole edge [walls.c:168](../Assignment_4/code/walls.c#L168), fill [initialise.c:125](../Assignment_4/code/initialise.c#L125), collision test [initialise.c:49](../Assignment_4/code/initialise.c#L49). Analysis in Python: [analyse_a.py](../Assignment_4/code/analyse_a.py), [analyse_c.py](../Assignment_4/code/analyse_c.py), [analyse_d.py](../Assignment_4/code/analyse_d.py).
- Your notebook already has "Oral prep" boxes after A, B, C, D: read those (they are good).

### A.4 Sample questions with answers
**(L3) Constants from restitution, and why that way round?**
Normal contact = damped oscillator m_eff x'' + η x' + k x = 0 with m_eff = m/2 (pair) or m (wall). Solution x = (v/ω_d) e^(−βt) sin ω_d t: separation at t_c = π/ω_d; rebound speed v e^(−βπ/ω_d) → e = exp(−βπ/ω_d). Hence t_c = sqrt(m_eff (π² + ln²e)/k), η = −2 ln e · m_eff/t_c. Reason: e is the measurable, material-level quantity (and flow is insensitive to the exact spring value); a physical k (Hertz, steel/glass) is far too stiff (tiny t_c, impossible dt), so you pick a soft k, keep e exact, and t_c sets dt (contact must be resolved by ≥ 10–20 steps). Note: in your code you pick k_n = 1000 N/m and t_c is derived; the examiner may phrase it "set e and t_c": equivalent, one of (k, t_c) is free.

**(L3) Why does tangential force need memory?**
Static friction is a spring on the accumulated tangential displacement since the contact formed; Coulomb cap μ|F_n| decides stick vs slide. Depends on history, not on current positions (MD forces are functions of current r only). Without stored displacement friction would act only on sliding speed: zero for a bed at rest → a bed (or heap) could not hold at any slope, and your tilt test would creep immediately. In code: `tijs[k]` (pp) and `tiw` (pw), reset to −F_t/k_t when sliding.

**(L2) Is energy conserved? What is the diagnostic still worth?**
No: dashpots dissipate (e < 1) and friction sliding dissipates. E_kin + E_pot decreases. Still useful: (i) E_kin → ~0 shows the fill has settled before t_open (5.6e-7 of max), (ii) E_kin onset defines the tilt angle, (iii) sudden unphysical increase flags bugs (bad overlap, too-large dt), (iv) in the A tests e² = E_kin,after/E_kin,before.

**(L3) Constant discharge while fill drops; why different from water?**
Janssen: grain contacts form force chains that push sideways, wall friction carries most of the weight, so pressure at the bottom saturates beyond ~1 diameter of fill; flow set by free-fall region ~D above orifice, Q = Cρ_b√g(D − k d_p)^(5/2), independent of head. Water: p = ρ g h, Torricelli v = sqrt(2gh), Q ∝ √h; your halves would differ by sqrt(3927/5799) = 0.82 (−18 %); you got −1.5 % (1.2× its error).

**(L2) Why is measuring e from an impact a genuine verification?**
e is an INPUT to the constants; the run measures it independently from E_kin before/after (e² = ratio), so recovering 0.85/0.92 (at three speeds) checks the entire chain: constants, force routine, integrator, dt, overlap detection. Speed independence (exact linear model) is a second test: a spread over speeds = time-integration error. Limit: normal direction only (no spin): tangential spring/friction not tested by Part A.

### A.5 Extra likely questions
- Why is t_c/20 enough? 10 steps already <1 % at all speeds; but phase effect: single collision up to +1.5 %; mean correct.
- Why wall ≠ pair constants? m_eff differs (m vs m/2): wall k = 2×; e differs (0.85 vs 0.92) → η 3.9×.
- Why does the bed fail at 12°, not 24.2°? Floor contacts roll rather than slide (speed ratio 0.57 ≈ 0.5), so μ_w does not set θ_c. **Hypothesis, untested**: did not vary μ_w; no spins in output.
- Edge of hole: closest point on rim at same azimuth, r_iw = ((r_c − R_o) e_r, z − z_s); proof by triangle inequality; checked vs brute force 685,455 particles.
- Why 16 % and not ~50 %? Hole effectively wider by ~0.3 d_p not 1 d_p; extra grains cross in the slow outer ring (84 % at r > 5 mm). Hypothesis (rim = slow zone) untested.
- Why does Beverloo have exponent 5/2? Q = ρ_b^a g^b D^c: mass a = 1, time b = ½, length −3a + b + c = 0 → c = 5/2 (2D: 3/2). Physical: free fall speed √(gD) × area D² × ρ_b.
- Which constant do you trust, k or C? k (ratio of fit parameters; ρ_b, m, g cancel; within 1.8σ). C depends on ρ_b of flowing material; 4σ low; candidates: flowing density < settled density, silo only 2 orifices wide at D = 9 d_p, exponent fixed at 5/2.
- D = 3 d_p: arch of mutually supporting grains blocks; friction + geometry; random event (one run); Beverloo has no arching mechanism so valid only for D ≳ 6 d_p (fitted range).
- Packing: φ = 0.597 (random loose packing between sc 0.52 and fcc 0.74); wall layers: oscillation period 0.86 d_p, decays over ~3.5 d_p (four layers), so silo must be ≫ 7 d_p wide for D = 9 d_p.
- L4: untested gaps: tangential/friction model not verified in isolation (design: oblique impact with spin, or measure angle of repose vs μ); vary μ_w to test the rolling explanation; multiple seeds for all D; 2×-wide silo for k, C systematics; add rolling friction.

---------------------------------------------------------------------
## 6. Cross-cutting (V&V segment and D4 live here)

### 6.1 V&V matrix: what each test proves and what it cannot catch

| Test | Proves | Cannot catch |
|---|---|---|
| A1 E, L, COM conservation | Internal consistency; integrator symplectic; central pairwise force | Wrong force law/G/masses/initial data; any force that is central + antisymmetric (e.g. wrong exponent) |
| A1 Euler vs Verlet order (slopes 0.71 / 1.99) | Order of method | Physics errors |
| A1 Horizons comparison | Model + data + implementation vs independent truth (validation); where error plateaus (model, not dt) | Which of several possible causes; needs ablation (no-Moon run) |
| A1 Apophis dt-convergence + Richardson | Numerical convergence | Model error (Yarkovsky etc.): converged ≠ correct |
| A2 FD force test | Force = −∇(coded U) | Wrong parameters/exclusions/masses/units; neighbour-list misses |
| A2 per-type temperature | Masses used per type | Everything else |
| A2 NVE drift ∝ dt² | Integrator + consistent forces | Wrong potential; cutoff force jump gives non-dt² slope |
| A2 torsion ΣF = 0, torque = 0 | Translation/rotation invariance of torsion | Correct magnitude of torsion force (FD test does that) |
| A2 MSD ballistic limit (3kT/M t²) | MSD code and units | Force field adequacy |
| A2 D vs experiment | Model validation (22–32 % high: model limitation, not bug: no factor 2, slope 1) | Source of discrepancy (not tested; experimental value unverified) |
| A3 B1 conservative-only energy (dt^2.33) | Conservative force/potential consistent, integrator | Dissipative/random part |
| A3 B2 ideal-gas g(r) = 1 | g(r) normalisation; dissipative+random don't create structure; conservative really off; pair finding | Conservative-force correctness |
| A3 B3 velocity variance (T = 1.03) | FDT σ² = 2γkT, 1/√Δt, integrator heating (+3 %) | Momentum conservation, conservative force |
| A3 C1 g(r) vs paper (digitised by eye) | Conservative force in a fluid | Unlike pairs, precision (paper points read by eye ±0.05–0.1) |
| A3 C3 χ slope 0.286 | Phase separation + χ measure vs G–W Eq. 24 (validation vs literature) | Same code-wide systematic errors if G–W had them too |
| A3 C4 σ_s control = 0 and ratio to Eq. 36 | No x/yz bias in pressure tensor; magnitude and factor ½ | Control alone can't catch a wrong overall factor, hence the magnitude check |
| A4 A restitution recovery (3 speeds × 4 dt) | Normal contact constants, force routine, dt | Tangential/friction part; non-smooth dt behaviour (phase effect) |
| A4 C wall check vs brute force | Wall geometry (closest point, flags) | Dynamics of contact |
| A4 B slope onset | Static friction memory (sticks until ~12°) | Why it fails there (rolling is a hypothesis; μ_w not varied) |
| A4 D2 halves ratio, D3 Beverloo fit | Constant Q; Beverloo scaling in 6–9 d_p; k within 1.8σ | C 4σ low unexplained; one seed; silo only 2 orifices wide at 9 d_p |

Generic sentence for any test: "It verifies that ___ is implemented as intended; it cannot detect ___, which I cover with ___."
Verification = "do I solve my equations right?" (conservation, FD, convergence, controls). Validation = "are these the right equations/parameters?" (Horizons, experiment, G–W, Beverloo literature).

### 6.2 Methods landscape (D4)

| | MD (A2) | DPD (A3) | DEM (A4) |
|---|---|---|---|
| Particle | atom/united atom | bead = cluster of molecules (coarse-grained) | real grain (mm) |
| Forces | conservative from potential (LJ, bonds…) | soft conservative + pairwise dissipative + random (thermostat) | contact only: spring, dashpot, friction with history |
| Energy | conserved (NVE) | not conserved; kT set by FDT | dissipated (e < 1), no temperature |
| Momentum | conserved globally | conserved locally per pair → hydrodynamics | conserved in collisions; gravity/walls external |
| Time step | fs (set by bond vibration/LJ) | 0.04 time units (soft potential allows large dt) | t_c/20 (resolve contact) |
| Strength | chemistry, conformers, transport of molecules | mesoscale: phase separation, polymers, surfactants, with hydrodynamics | granular flow, packing, silos |
| Limits | tiny length/time scales | no chemistry; Schmidt number ~1 (liquids ~10³); soft, compressible | no thermal motion; contact law models; cost grows with N and stiffness |

How to choose: scales of length and time; are thermal fluctuations relevant (Péclet = convection/diffusion; Pe ≲ 1–30 → Brownian effects matter → mesoscale/Langevin/DPD/SRD; Pe ≫ 1 → deterministic, continuum CFD/DEM–CFD); Reynolds number (Stokes vs inertial); need hydrodynamics (momentum-conserving methods: DPD, SRD, LB; not Langevin); is there a solvent to resolve.
One-liners from old exams: **Langevin dynamics** = implicit solvent, friction + random force, no momentum conservation, no hydrodynamic interactions, cheap for dilute colloids. **DPD** as above. **SRD/MPCD** = solvent particles streaming and periodically rotating velocities in cells about a random axis, conserves momentum and energy locally, gives hydrodynamics, colloids/polymers; parameters set by dimensionless numbers (viscosity, Schmidt, Pe, Re) not by mapping to real SI units (naive water mapping needs absurdly small cells).
Worked transfer example (2018 exam, blood centrifuge) **[check numbers]**: v = (2/9)(ρc − ρp) ω²L R²/μ ≈ 1.2e-3 m/s; D = kT/(6πμR) ≈ 5.6e-14 m²/s; Pe = vR/D ≈ 6e4, Re = ρvR/μ ≈ 3e-3. Re ≪ 1 (Stokes) but Pe ≫ 1: sedimentation dominates over Brownian motion, so a fluctuation-resolving mesoscale method would be wasteful and hard to map (matching a huge Pe is expensive); use continuum Stokes/Lattice-Boltzmann or a deterministic coarse-grained/Stokesian dynamics model.
Your own Lessons learned (A4) says you were unsure what "which method first" meant: prepare a real answer: "Start from the scale of the question: atomistic detail → MD; hydrodynamics + thermal fluctuations at 10 nm–µm → DPD/SRD; grains and no thermal motion → DEM; always begin with the simplest version, validate against a known limit, then add physics."

### 6.3 Bug stories you can tell (from your reports)
- A1: after C3 you found Euler's L drifts O(dt) (Δt² Σ m v × a term), and corrected your A4 table that said L drift always signals a force-loop bug. The COM-velocity drift turned out to be finite-difference round-off.
- A2: B7 temperature dropped to 203 K → lattice not melted; every later run started from it; you'd equilibrate with long NVT first. PowerShell redirect wrote UTF-16 CSV. Didn't log early runs.
- A3: g(r) not flat in the first ideal-gas runs (normalisation); B1 energy "wandered by 1 % independent of dt": noise force was still on; found because the dt-scaling exposed it; regenerated with `run_b1.sh`.
- A4: first wall check "failed" because the brute-force reference could not resolve azimuth near the axis, not the code; one seed gave 17 % vs 16 ± 2 % over five seeds; empty output files from a missing timing tool.

### 6.4 What would you do differently (prepare a 30 s answer)
Record every run (command, parameters) as a script from the start; equilibrate properly (A2); several seeds/replicas for every number (A4, A3); run the control experiments you only argued (A2: ε33 = 0 control; A4: vary μ_w); verify literature reference values against primary sources (A2 D_exp).

---------------------------------------------------------------------
## 7. Predict drills (say the answer in one sentence before reasoning)

A1
- Double dt: Verlet energy error → ×4; Euler → ×~1.6–2; Verlet L error stays ~0, Euler L ×2.
- Remove the Moon from the 2029 run: Earth's initial velocity (12 m/s wobble) is no longer balanced → Earth ends ~10^6 km off → flyby destroyed (your run: 2.15 million km). **[check]**
- Use planet-centre states for Pluto: ~10^8 m error after a year; no dt helps.
A2
- Halve dt: E_tot oscillation ÷4. Switch to NVT: E not conserved. Raise T by 100 K: trans ↓ (≈ 0.72 → 0.63 **[check]**), D ↑. Turn off 1–5 LJ: g+g+ ≈ g+g−. Forget 1–2/1–3 exclusion: huge repulsion, bonds stretched. Wrong mass for CH3: T_CH3 ≠ T_CH2 (the one test that sees it). Density ↑ → D ↓.
A3
- γ doubled with σ fixed: kT halves. σ doubled with γ fixed: kT ×4. Remove 1/√Δt: kT collapses (≪ 1) and depends on Δt. Δa ↑ → χ ↑ ∝ 0.286 Δa, σ_s ↑, φ_A in the wrong phase ↓. N ↑ at fixed Δa → χN ↑ ∝ N, demix more easily. Switch conservative off → g(r) = 1. Add Berendsen on top → suppressed fluctuations, FDT broken.
A4
- Orifice 7 → 8 d_p: Q ratio = ((8 − 1.59)/(7 − 1.59))^2.5 ≈ 1.5. Fill twice as high: Q unchanged. Silo twice as wide (if already ≫ D): Q unchanged (but at D = 9 d_p wall effects change it). g reduced to 1/6: Q × 0.41. Grain size ↑ at fixed D: (D − k d_p) shrinks → Q ↓. Friction ↑: arches more likely, min orifice ↑. Δt ×2 (t_c/10): e error ~±1 %, phase-dependent. e ↓: η ↑, rebound lower, bed settles faster. Tilt rate ×2: onset angle slightly later (0.4° in your data).

---------------------------------------------------------------------
## 8. Honest weak spots in YOUR reports (the examiner can see these)

### 8.1 Things that are wrong, thin or unfinished
1. **A1 "Lessons learned" is empty** ("*Lessons learned:* …"). A standing question in the oral; prepare it (see §6.4).
2. **A1 C5, no-Moon explanation.** Report says "chaos/keyhole". The task hint: the effect is larger than the Moon's own pull on Apophis, so "what dominates?" **[check]** Better, and consistent with your own C1 note: the snapshot's Earth velocity includes the ~12 m/s monthly wobble about the Earth–Moon barycentre; without the Moon this becomes an error in Earth's heliocentric orbit: 12 m/s × 2.6 yr ≈ 1e9 m ≈ 10^6 km, same order as your 2.1e6 km. The flyby is extremely sensitive to Earth's position.
3. **A1 C1** "Earth–Moon envelope period ≈ 200 days" is stated without derivation; don't defend a number you can't derive. Safe version: the Sun perturbs the Moon's orbit (solar tide), periodic modulation, not a bug since E is flat at 1e-9.
4. **A1 C2** Pluto/Sun errors plateau (138 m, 231 m): the report blames "initial conditions / barycentre representation" vaguely. Safe statement: dt is no longer the limiting factor; remaining error is model/data (physics not in my model e.g. GR, asteroids; representation), not integration. Say how you'd test (add a body/GR).
5. **A2 "g+g+ vs g+g−"**: causal claim (1–5 LJ) is argued from geometry, **control with ε33 = 0 was not run**; D_exp from an unchecked web source; 4-replica error on g+g− probably too small (they say ±10 % in the oral).
6. **A2 "why not lattice for FD test"**: answer is thin; use the fuller one in §A2.4.
7. **A3 notebook typo:** U^C written as (a rc/2)(1 − r/(2 rc))² in B1; correct is (1 − r/rc)². Code is right.
8. **A3 A1 notes** don't mention 1/√Δt at all; **A3 pressure claim [check]**: you wrote the control pressure 23.7 "as the EOS p = ρkT + 0.101 a ρ² predicts", but that EOS gives 3 + 0.101·25·9 = **25.7**, 8 % higher. Either explain (kinetic part is 3.09 at T = 1.03; α fit from G–W; you used ρ = 3) or don't say they agree exactly. Be ready to either defend or admit.
9. **A3** `gauss()` is uniform (variance 1) and uses `rand()`; B3 5 % temperature error not tested vs dt; r_c = 1 hard-coded in the weight (wR = 1 − r_ij) and r_shell = 1 in the committed `setparameters.c` vs 0.4 in your runs.
10. **A4**: normal restitution is verified, but the tangential model (e_t = 0.33, friction) is not verified on its own; "rolling" explanation and Beverloo C deviation are explicitly untested hypotheses; D3 error from only 4 window pieces (~30 % uncertain).

### 8.2 Read-and-narrate list (AI-written parts you MUST be able to walk through)
Priority order (D1 is the fail-gate):
1. A4: `cylinder_wall` and `silo_bottom_wall` (walls.c:134, 168), `initialise_positions_silo` (initialise.c:125), `count_in_silo`, `initialise_collision`, RUN_SLOPE branch of `update_parameters`, and what `analyse_a.py/c.py/d.py` do (restitution from energies, φ core, Q window fit).
2. A3: pressure-tensor accumulation (forces.c virial_xyz; dynamics.c kinetic part), `sample_profile` / `sample_pressure` / `sample_rdf` in analysis.c, run scripts.
3. A2: `analysis_update`, `sample_dihedrals`, `sample_msd` (analysis.c), `run_c.sh`; speed-up in `calculate_forces_nb` (pair-parameter table).
4. A1: main.c CLI parsing (declared AI) and the notebook cells for C1–C5.
Method: open the function, say aloud (a) inputs, (b) what it computes, (c) which test exercises it, (d) one thing that could go wrong. If you can't, ask me to explain it line by line.

---------------------------------------------------------------------
## 9. Last-hour checklist
- [ ] Four notebooks open in tabs (A1.ipynb, A2.ipynb, A3.ipynb, A4.ipynb or their HTML), code folder open in the editor, binaries built (A1 main.exe, A2/A3 md/dpd, A4 dem).
- [ ] Showcase chosen and rehearsed out loud (≤ 3 min), with its follow-up.
- [ ] For each assignment: 30-second story, 3 numbers, 1 bug story, the V&V "proves/can't catch" sentence.
- [ ] Code locations: forces, integrator, thermostat/DPD forces, contact forces, wall functions.
- [ ] Prepared lessons-learned (what differently, what surprised, what still bothers).
- [ ] D4 answer: MD vs DPD vs DEM, when each, what transfers.
- [ ] Be able to run one thing live (e.g. re-plot a figure, compile and run a short A1 case, or the A4 test collision) in case of the "live action".
