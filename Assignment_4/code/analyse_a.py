"""Part A: restitution coefficient from the runs of run_a.sh and run_a_scan.sh.

Run from Assignment_4/:   python code/analyse_a.py
Reads runs/A_<test>_v<v>_f<frac>/run.dat and runs/A_scan_<test>_f<frac>_v<v>/run.dat
(columns: step, time, Z, Epot, Ekin, Etot, N_silo) and writes
  data/A_restitution.csv   the 24 runs of the table
  data/A_phase_scan.csv    e against the contact phase (scan over the impact speed)
  data/A_A1_v0.3.dat, data/A_A2_v0.3.dat   table of the 0.3 m/s, dt = t_c/20 runs (figure)

Definitions (Epot is the energy in the normal spring, > 0 only while the
particles overlap):
  Ekin_before = Ekin at the last step before the first step with Epot > 0
  Ekin_after  = Ekin at the first step after the last step with Epot > 0
                (Epot is 0 at every later step, checked)
  e = sqrt(Ekin_after / Ekin_before)   (head-on, no rotation, no gravity)
  phase = fraction of a time step by which the first grid time with overlap
          lies after the moment the surfaces first touch. The particles start
          0.05 R from touching (A1: wall, speed v; A2: gap 0.1 R, closing speed
          2 v), so touching happens at t_touch = 0.05 R / v, the first overlap is
          at step n1 = ceil(t_touch / dt) and phase = n1 - t_touch / dt, in
          (0, 1]. The first overlap there is phase * v * dt. Checked against the
          first step with Epot > 0 in every run.
"""
import glob
import os
import shutil
import numpy as np

R_MIN = 0.975e-3                        # setparameters.c
E_TAB = {1: 0.85, 2: 0.92}              # particle-wall, particle-particle
SPEEDS = ["0.1", "0.3", "1.0"]
FRACS = ["0.2", "0.1", "0.05", "0.025"]


def analyse(path, test, v):
    """Return Ekin_before, Ekin_after, e, phase and the number of steps with overlap."""
    with open(path) as f:
        dt = float(f.readline().split()[-1])        # '# run ... dt <value>'
    d = np.loadtxt(path)
    epot, ekin = d[:, 3], d[:, 4]
    touching = np.flatnonzero(epot > 0)
    first, last = touching[0], touching[-1]
    assert first >= 1, f"{path}: contact already at the first printed step"
    assert last + 1 < len(epot), f"{path}: contact not finished"
    assert np.all(epot[last + 1:] == 0), path
    t_touch = 0.05 * R_MIN / v
    n1 = int(np.ceil(t_touch / dt))
    assert n1 == d[first, 0], f"{path}: first overlap at step {d[first, 0]}, expected {n1}"
    phase = n1 - t_touch / dt
    ek_before, ek_after = ekin[first - 1], ekin[last + 1]
    return ek_before, ek_after, np.sqrt(ek_after / ek_before), phase, last - first + 1


os.makedirs("data", exist_ok=True)

# the 24 runs of the table
with open("data/A_restitution.csv", "w") as f:
    f.write("num_test,v_impact,dt_fraction,steps_per_contact,Ekin_before,Ekin_after,"
            "e,e_tabulated,rel_error,phase,steps_with_overlap\n")
    for test in (1, 2):
        for v in SPEEDS:
            for frac in FRACS:
                eb, ea, e, ph, n = analyse(f"runs/A_{test}_v{v}_f{frac}/run.dat", test, float(v))
                f.write("%d,%s,%s,%d,%.9e,%.9e,%.6f,%.2f,%.5f,%.4f,%d\n"
                        % (test, v, frac, round(1 / float(frac)), eb, ea, e, E_TAB[test],
                           e / E_TAB[test] - 1, ph, n))
                print("A%d v=%-4s f=%-5s e=%.4f (tab %.2f) rel.err %+.2f %%  phase %.3f  overlap steps %d"
                      % (test, v, frac, e, E_TAB[test], 100 * (e / E_TAB[test] - 1), ph, n))

for test in (1, 2):
    shutil.copy(f"runs/A_{test}_v0.3_f0.05/run.dat", f"data/A_A{test}_v0.3.dat")

# the scan over the impact speed (if run_a_scan.sh was run)
paths = sorted(glob.glob("runs/A_scan_*/run.dat"))
if paths:
    with open("data/A_phase_scan.csv", "w") as f:
        f.write("num_test,dt_fraction,v_impact,e,e_tabulated,rel_error,phase,steps_with_overlap\n")
        for test in (1, 2):
            for frac in FRACS:
                for p in sorted(glob.glob(f"runs/A_scan_{test}_f{frac}_v*/run.dat")):
                    v = p.split("_v")[-1].split("/")[0].split("\\")[0]
                    eb, ea, e, ph, n = analyse(p, test, float(v))
                    f.write("%d,%s,%s,%.6f,%.2f,%.5f,%.4f,%d\n"
                            % (test, frac, v, e, E_TAB[test], e / E_TAB[test] - 1, ph, n))
    print("wrote data/A_phase_scan.csv from", len(paths), "runs")
print("wrote data/A_restitution.csv, data/A_A1_v0.3.dat, data/A_A2_v0.3.dat")
