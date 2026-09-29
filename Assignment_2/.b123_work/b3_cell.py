# B3: analyse saved C simulation results. Run this notebook from Assignment_2.
# The C runs produce the data; this cell loads them and makes every result below.
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress, t as student_t
from IPython.display import display

time_unit_ps = np.sqrt(1.66053906660e-27 * 1e-20 / 1.380649e-23) / 1e-12
N = 2000
target = 293.0
total_mass = 800*15.035 + 1200*14.027

b1 = pd.read_csv("data/b1_output_check.csv")
assert np.isfinite(b1.to_numpy()).all() and np.all(np.diff(b1.step) == 10)
print("B1 max |Etotal-Epotential-Ekinetic|:",
      np.max(np.abs(b1.Etot_internal-b1.Epot_internal-b1.Ekin_internal)))
print("B1 max temperature identity error:",
      np.max(np.abs(b1.T_internal-2*b1.Ekin_internal/(3*N-3))))

# These literal paths also let the submission checker find every input file.
equil_files = __EQUIL_FILES__
nve_1fs = pd.read_csv("data/b3_final_nve_1fs.csv")
nve_half = pd.read_csv("data/b3_final_nve_05fs.csv")
force_logs = [Path("data/b3_final_fd_1e-06.txt").read_text(),
              Path("data/b3_final_fd_1e-05.txt").read_text()]
meta_files = __META_FILES__
metadata = [json.loads(Path(name).read_text()) for name in meta_files]
restart = Path("data/b3_final_equil___STAGE___restart.dat")
restart_hash = hashlib.sha256(restart.read_bytes()).hexdigest()
for run in metadata:
    if run["protocol"] != 0:
        assert run["input_sha256"] == restart_hash, "Tests did not use the same restart!"
    assert run["returncode"] == 0

equil_parts, offset = [], 0.0
for name in equil_files:
    frame = pd.read_csv(name)
    frame["time_ps"] = frame.time_internal*time_unit_ps + offset
    offset = frame.time_ps.iloc[-1]
    equil_parts.append(frame)
equil = pd.concat(equil_parts, ignore_index=True)
for frame in [equil, nve_1fs, nve_half]:
    assert np.isfinite(frame.to_numpy()).all()
    assert np.allclose(frame.Etot_internal, frame.Epot_internal+frame.Ekin_internal,
                       rtol=1e-12, atol=1e-8)
    assert np.allclose(frame.T_internal, 2*frame.Ekin_internal/(3*N-3), rtol=1e-12)
    assert np.allclose(800*frame.T_CH3_K+1200*frame.T_CH2_K,
                       (N-1)*frame.T_internal, rtol=1e-12)

def blocks_and_trend(time_ps, values, number=10):
    """Contiguous block means reduce correlation before estimating trend error."""
    indices = np.array_split(np.arange(len(values)), number)
    bt = np.array([np.mean(np.asarray(time_ps)[part]) for part in indices])
    bv = np.array([np.mean(np.asarray(values)[part]) for part in indices])
    fit = linregress(bt, bv)
    half_ci = student_t.ppf(0.975, number-2)*fit.stderr
    return bt, bv, fit.slope, half_ci

# Stationarity: final 50 ps of preparation plus the separate 50 ps NVT
# validation branch, using twenty 5 ps blocks. The shared restart is at 100 ps.
window_start = equil.time_ps.iloc[-1]-100.0
late_equil = equil.loc[equil.time_ps > window_start+1e-8].copy()
assert len(late_equil) == 1000
trend_rows = []
for field, label in [("Epot_internal", "Potential energy"),
                     ("T_CH3_K", "CH3 temperature"), ("T_CH2_K", "CH2 temperature")]:
    bt, bv, slope, half_ci = blocks_and_trend(late_equil.time_ps, late_equil[field], 20)
    trend_rows.append({"quantity": label, "mean": bv.mean(),
                       "block standard error": bv.std(ddof=1)/np.sqrt(len(bv)),
                       "slope per ps": slope, "95% slope CI half-width": half_ci,
                       "slope CI contains zero": abs(slope) <= half_ci})
trend_table = pd.DataFrame(trend_rows)
print(f"Equilibration window: {window_start:.0f}–{equil.time_ps.iloc[-1]:.0f} ps of the new continuation.")
display(trend_table)

fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, constrained_layout=True)
for field, label in [("T_CH3_K", "CH3"), ("T_CH2_K", "CH2")]:
    axes[0].plot(equil.time_ps, equil[field], lw=0.5, alpha=0.5, label=label)
axes[0].axhline(target, color="black", ls="--", label="293 K")
axes[0].set_ylabel("Site-type temperature (K)")
axes[0].legend()
axes[1].plot(equil.time_ps, equil.Epot_internal, lw=0.6, alpha=0.7, label="Potential energy")
for i, part in enumerate(equil_parts):
    bt, bv, _, _ = blocks_and_trend(part.time_ps, part.Epot_internal, len(part)//50)
    axes[1].plot(bt, bv, "o-", color="black", ms=3, label="5 ps means" if i == 0 else None)
axes[1].set(xlabel="Additional NVT time (ps)", ylabel="Potential energy / (kB K)")
axes[1].legend()
for ax in axes:
    ax.axvspan(window_start, equil.time_ps.iloc[-1], color="green", alpha=0.1)
    ax.axvline(100, color="purple", ls=":", lw=1)
    ax.grid(alpha=0.2)
fig.suptitle("Equilibration at 626 kg/m³; shaded region used for stationarity and mass checks")
plt.show()

# Finite differences: preserve the supplied checker's error definitions.
fd_rows, fd_details = [], []
for text in force_logs:
    header = re.search(r"# U=(\S+) delta=(\S+) h=(\S+) N=(\d+)", text)
    energy, delta, h = map(float, header.groups()[:3])
    errors = [float(x) for x in re.findall(r"Particle .*?rel.err=(\S+)", text)]
    ids = [int(x) for x in re.findall(r"Particle (\d+):", text)]
    assert ids == list(range(N)) and len(errors) == N and "not tested" not in text
    force_norms = [np.linalg.norm([float(x) for x in xyz.split(',')])
                   for xyz in re.findall(r"Force analytical: \(([^)]+)\)", text)]
    virial = re.search(r"Virial test: analytical W = (\S+), finite-difference W = (\S+), rel.err=(\S+)", text)
    W, Wfd, Werror = map(float, virial.groups())
    force_floor = 16*np.finfo(float).eps*abs(energy)/delta
    relative_floor = force_floor/min(force_norms)
    fd_rows.append({"check": "maximum force, 2000 sites", "FD step": f"{delta:g} Å",
                    "relative error": max(errors), "tolerance": 1e-4, "pass": max(errors) < 1e-4})
    fd_rows.append({"check": "global virial", "FD step": f"{h:g} (dimensionless)",
                    "relative error": Werror, "tolerance": 1e-8, "pass": Werror < 1e-8})
    fd_details.append({"delta": delta, "h": h, "energy": energy, "min_force": min(force_norms),
                       "absolute_floor": force_floor, "relative_floor": relative_floor,
                       "max_force_error": max(errors), "virial_error": Werror,
                       "virial_floor": 16*np.finfo(float).eps*abs(energy)/(h*max(abs(Wfd),abs(energy)))})
display(pd.DataFrame(fd_rows))
display(pd.DataFrame(fd_details).drop(columns=["max_force_error", "virial_error"]))

# NVE: fit the secular component, then measure the fluctuations about that fit.
fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex="col", constrained_layout=True)
drift_rows, temperature_rows, nve_trends = [], [], []
for col, (label, frame) in enumerate([("1 fs", nve_1fs), ("0.5 fs", nve_half)]):
    time_ps = frame.time_internal.to_numpy()*time_unit_ps
    energy = frame.Etot_internal.to_numpy()
    slope, intercept = np.polyfit(time_ps, energy, 1)
    fit = slope*time_ps+intercept
    residual = energy-fit
    drift_rows.append({"timestep": label, "duration (ps)": time_ps[-1],
                       "energy slope (internal/ps)": slope,
                       "relative drift (1/ps)": slope/abs(energy.mean()),
                       "detrended energy RMS": residual.std(ddof=0),
                       "relative detrended RMS": residual.std(ddof=0)/abs(energy.mean())})
    for field, legend in [("Epot_internal", "Potential"), ("Ekin_internal", "Kinetic"), ("Etot_internal", "Total")]:
        axes[0,col].plot(time_ps, frame[field], label=legend, lw=0.7)
    axes[0,col].set(title=f"NVE, {label}", ylabel="Energy / (kB K)")
    axes[0,col].legend()
    axes[1,col].plot(time_ps, energy-energy.mean(), lw=0.6, label="Total energy minus mean")
    axes[1,col].plot(time_ps, fit-energy.mean(), "k--", label="Linear fit")
    axes[1,col].set_ylabel("Energy deviation / (kB K)")
    axes[1,col].legend()
    axes[2,col].plot(time_ps, residual, lw=0.6)
    axes[2,col].axhline(0, color="black", lw=0.6)
    axes[2,col].set(xlabel="Time (ps)", ylabel="Detrended energy / (kB K)")
    for field, species, mass in [("T_CH3_K", "CH3", 15.035), ("T_CH2_K", "CH2", 14.027)]:
        _, bv, b, ci = blocks_and_trend(time_ps, frame[field])
        temperature_rows.append({"run": f"NVE {label}", "site": species,
                                 "mean temperature (K)": bv.mean(), "block SE (K)": bv.std(ddof=1)/np.sqrt(10),
                                 "set point (K)": target,
                                 "COM-corrected target (K)": target*(1-mass/total_mass)})
        nve_trends.append({"run": label, "site": species, "slope K/ps": b, "95% CI half-width": ci})
for ax in axes.flat:
    ax.grid(alpha=0.2)
fig.suptitle("Matched 50 ps NVE runs from the same equilibrated restart")
plt.show()
drift_table = pd.DataFrame(drift_rows).set_index("timestep")
display(drift_table)
fluctuation_ratio = drift_rows[0]["detrended energy RMS"]/drift_rows[1]["detrended energy RMS"]
drift_ratio = abs(drift_rows[0]["relative drift (1/ps)"]/drift_rows[1]["relative drift (1/ps)"])
print(f"Halving dt: detrended RMS reduction = {fluctuation_ratio:.4f}; drift magnitude reduction = {drift_ratio:.4f}.")
for field, species, mass in [("T_CH3_K", "CH3", 15.035), ("T_CH2_K", "CH2", 14.027)]:
    _, bv, _, _ = blocks_and_trend(late_equil.time_ps, late_equil[field], 20)
    temperature_rows.insert(0, {"run": "NVT stationary 100 ps", "site": species,
                                "mean temperature (K)": bv.mean(), "block SE (K)": bv.std(ddof=1)/np.sqrt(len(bv)),
                                "set point (K)": target, "COM-corrected target (K)": target*(1-mass/total_mass)})
display(pd.DataFrame(temperature_rows))
display(pd.DataFrame(nve_trends))

fig, axes = plt.subplots(1, 2, figsize=(11, 3.5), constrained_layout=True)
for ax, (label, frame) in zip(axes, [("1 fs", nve_1fs), ("0.5 fs", nve_half)]):
    time_ps = frame.time_internal.to_numpy()*time_unit_ps
    for field, name in [("T_CH3_K", "CH3"), ("T_CH2_K", "CH2")]:
        bt, bv, _, _ = blocks_and_trend(time_ps, frame[field])
        ax.plot(bt, bv, "o-", label=name)
    ax.axhline(target, color="black", ls="--", label="293 K")
    ax.set(title=f"NVE {label}", xlabel="Time (ps)", ylabel="5 ps mean temperature (K)")
    ax.legend()
    ax.grid(alpha=0.2)
plt.show()

display(pd.DataFrame([{"protocol": m["protocol"], "stage": m["stage"],
                       "FD delta": m["fd_delta"], "wall time (s)": m["wall_seconds"],
                       "input restart": m["input_restart"]} for m in metadata]))
print("Independent typing, mass-update and all-pairs check:")
print(Path("data/b3_independent_check_final.txt").read_text(encoding="utf-8-sig"))

# Used to insert the measured numbers into the accompanying explanation.
b3_summary = dict(equilibration=trend_rows, fd=fd_details, drift=drift_rows,
                  temperature=temperature_rows, nve_trends=nve_trends,
                  fluctuation_ratio=fluctuation_ratio, drift_ratio=drift_ratio,
                  equilibration_ps=float(equil.time_ps.iloc[-1]), window_start=float(window_start),
                  common_restart_sha256=restart_hash)
