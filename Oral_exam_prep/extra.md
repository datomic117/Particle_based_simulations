These are examples of the kind of questions asked in the oral examination, at the published depth levels (L1 describe → L2 explain → L3 predict & connect → L4 critique & design). In the exam, most questions are asked about your own code and results — the best preparation is to be able to answer these for your submission, not in the abstract.
A1 — Planetary motion

    (L2) Euler-forward's energy error grows steadily while velocity-Verlet's stays in a bounded band. Why?
    (L2) How does your code handle zero-mass bodies, and why can't you compute their acceleration as F/m?
    (L3) Your energy, angular momentum and center-of-mass diagnostics were all clean. What does the comparison against the JPL Horizons snapshot check that those tests cannot?
    (L3) Suppose your force magnitude had a bug but the direction was right. Which of your diagnostics would catch it — and which would stay perfectly clean, and why?
    (L3) What sets the largest usable timestep for your set of bodies, and what would you change if a body with a very eccentric orbit were added?

A2 — Molecular dynamics of n-pentane

    (L2) Why do 1-2 and 1-3 bonded pairs get no Lennard-Jones interaction? What would you observe if you forgot that?
    (L3) Why must the finite-difference force test be run from an equilibrated configuration rather than from the initial lattice?
    (L2) Why is energy conservation tested in NVE rather than NVT, and how should the residual drift scale with the timestep?
    (L3) You raise the temperature by 100 K: which way does your trans fraction move, and why?
    (L4) What does the finite-difference force test prove — and name a bug class it can never catch. What covers that gap?

A3 — Dissipative particle dynamics

    (L3) There is a 1/√Δt in the random force. Where does it come from, and what goes wrong without it?
    (L2) Why are σ and γ not independent parameters in your code?
    (L3) In the ideal-gas test (conservative force off), what exactly must come out, and which bug classes does this test isolate?
    (L2) Why did the Berendsen thermostat have to be removed when you converted your MD code to DPD?
    (L3) Both the DPD thermostat and Berendsen hold the temperature. For what kind of physics does DPD's momentum conservation matter?

A4 — DEM: silo discharge and the Beverloo law

    (L3) You never set the spring and dashpot constants directly — you set restitution coefficients and a contact time. Explain how the constants follow, and why we do it that way around.
    (L3) Why does the tangential force need memory (the accumulated tangential displacement), when nothing in MD ever needed contact history?
    (L2) Is energy conserved in your DEM simulation? What is the energy diagnostic still worth?
    (L3) Your discharge rate is constant while the fill height drops. Why — and why is a draining water tank different?
    (L2) Why is measuring the restitution coefficient from an impact simulation a genuine verification of the contact model?

And in every exam

    Walk me through your own code: any function you submitted, any plot in your report (“we didn't write that part” is not an answer).
    Predict before you reason: “if we change X, what happens to Y — commit to an answer, then justify it.”
    Tell me about a bug you hit, what the diagnostics showed, and how you found it.
    What would you do differently if you started over?


Individual Oral Examination
How you are assessed

Your grade for 6EMA02 is determined in a 25-minute individual oral examination in weeks 9–10. The four assignments and their report notebooks are mandatory, but their result is determined in the oral — nobody grades the text of a report; the oral is what the report is about.

Each report carries a weight of 0.25: hand in all four and your final grade is simply your oral grade. A report you never hand in scores 0 for its weight (final grade = oral × 0.75 with one report missing). Handed-in reports stay valid for 3 years.
The 20 minutes

    Showcase (~3 min): you open with one result from your four assignments that you are proud of.
    Deep dive (~8 min): the examiner picks one of your four assignments and goes into your code and results. You don't know which one in advance — prepare all four.
    Verification & validation (~6 min): questions across the other assignments about how you verified and validated your simulations.
    Extension (~3 min): a final stretch to find how far your understanding reaches.

What we look for

The difference between having run a simulation and understanding it:

    predict before running — what happens if we double the timestep?
    explain why a test proves correctness — what does the force test catch, and what can it never catch?
    connect what you measured to the physics — why this conformer ratio?
    transfer — when MD, when DPD, when DEM?

The rubric

Four dimensions, each scored 0–4:
Dimension	What it measures
D1 Implementation command	Can you walk through and reason about your own code?
D2 Verification & validation	Do you know what each test proves, and what it cannot?
D3 Physical interpretation	Can you connect your results to the underlying physics, quantitatively?
D4 Methods & transfer	Which method for which problem — and why?

Score anchors per dimension: 0 = cannot engage · 1 = describes procedure only · 2 = explains own work correctly with prompting · 3 = explains unprompted, predicts correctly, connects · 4 = critiques, designs, transfers beyond the material.

Grade mapping (S = sum over the four dimensions): any dimension at 0, or D1 ≤ 1 → fail. S 6–7 → 6 · S 8–9 → 7 · S 10–12 → 8 · S 13–14 → 9 · S 15–16 → 10. A student who cannot explain their own submitted work fails, regardless of report quality.
Question depth levels

Questions escalate through four levels; your grade comes from the ceiling you reach, not from time survived:

    L1 describe — what did you do, what does this plot show;
    L2 explain — why is this correct, what does this test verify;
    L3 predict & connect — what changes if …, why does the physics come out this way;
    L4 critique & design — limits of the method, design a new test.

Practical

The exam is on campus at a table, and runs through a Teams meeting so it can be recorded (the required record for a single-examiner oral; the recording starts after the ID check). Join the meeting from your laptop and share your screen. Have ready: your four report notebooks open in tabs, your code open in an editor, your binaries built. Slot scheduling by availability poll in week 8.

Rehearsal: an in-class demo oral in week 4 and a peer-oral exercise in week 6 — the format will not be new to you at the graded moment. See also the page with sample questions.
AI and sources

You may use AI assistants freely in this course — for code, analysis, and text. Declare the use in the Sources section of each report notebook, together with any other external sources: anything in your submission your group did not produce itself. Declared sources are never penalised. Everything in your submission is fair game in the oral, and “we didn't write that part” is not an answer. Presenting undeclared external work as your own is academic misconduct.

The best preparation for this exam is the verification work the assignments already ask of you: if you understand why your tests pass, you understand your simulation.
