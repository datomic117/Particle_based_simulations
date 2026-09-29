"""B1-B3 analysis. Run from the notebook directory; all inputs are in data/."""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display

time_unit_ps = np.sqrt(1.66053906660e-27 * 1e-20 / 1.380649e-23) / 1e-12
b1 = pd.read_csv("data/b1_output_check.csv")
nve1 = pd.read_csv("data/b3_nve_1fs.csv")
nvehalf = pd.read_csv("data/b3_nve_05fs.csv")
equil = pd.read_csv("data/b3_equil_long.csv")
force_log = Path("data/b3_forcecheck.txt").read_text(encoding="utf-8")

# B1: check the actual saved output, rather than just describing its format.
energy_residual = np.max(np.abs(b1.Etot_internal - b1.Epot_internal - b1.Ekin_internal))
temperature_residual = np.max(np.abs(b1.T_internal - 2*b1.Ekin_internal/(3*2000-3)))
assert np.isfinite(b1.to_numpy()).all()
assert np.all(np.diff(b1.step) == 10)
assert energy_residual < 1e-8 and temperature_residual < 1e-12
print(f"B1: {len(b1)} finite rows, spacing 10 steps; max |E-K-U| = {energy_residual:.3g}; "
      f"max temperature identity error = {temperature_residual:.3g} K")

# The checker reports a vector-relative force error for each particle.
force_errors = np.array([float(v) for v in re.findall(r"Particle .*?rel.err=([^\s]+)", force_log)])
force_norms = np.array([np.linalg.norm([float(v) for v in xyz.split(',')])
                        for xyz in re.findall(r"Force analytical: \(([^)]+)\)", force_log)])
virial_error = float(re.search(r"Virial test:.*rel.err=([^\s]+)", force_log)[1])
assert len(force_errors) == 2000 and "not tested" not in force_log
# Scale estimate from the saved 20 ps NVT endpoint; the old log lacks run metadata.
earlier_equil = pd.read_csv("data/b3_equil_retry.csv")
energy_scale = abs(earlier_equil.Epot_internal.iloc[-1])
delta = 1e-6
absolute_floor = 16*np.finfo(float).eps*energy_scale/delta
relative_floor = absolute_floor/force_norms.min()
display(pd.DataFrame({
    "check": ["max force error, 2000 sites", "global virial error"],
    "FD step": ["delta = 1e-6 angstrom", "h = 1e-6 (dimensionless)"],
    "relative error": [force_errors.max(), virial_error],
    "accepted tolerance": [1e-4, 1e-8],
    "within tolerance": [force_errors.max() < 1e-4, virial_error < 1e-8],
}))
print(f"Force noise estimate: {absolute_floor:.3g} internal force units; "
      f"smallest |F| = {force_norms.min():.3g}; relative floor ~ {relative_floor:.3g}.")

fig, axes = plt.subplots(2, 2, figsize=(12, 7), constrained_layout=True)
drift_rows = []
for column, (label, frame) in enumerate([(1.0, nve1), (0.5, nvehalf)]):
    t = frame.time_internal.to_numpy()*time_unit_ps
    energy = frame.Etot_internal.to_numpy()
    slope, intercept = np.polyfit(t, energy, 1)
    residual = energy - (slope*t + intercept)
    for field, legend in [("Etot_internal", "total"), ("Ekin_internal", "kinetic"), ("Epot_internal", "potential")]:
        axes[0, column].plot(t, frame[field], label=legend, linewidth=0.8)
    axes[0, column].set(title=f"NVE, timestep {label:g} fs", ylabel="Energy / (kB K)")
    axes[0, column].legend()
    axes[1, column].plot(t, energy-energy.mean(), lw=0.7, label="total energy minus mean")
    axes[1, column].plot(t, slope*t+intercept-energy.mean(), color="black", label="linear fit")
    axes[1, column].set(xlabel="Time (ps)", ylabel="Energy deviation / (kB K)")
    axes[1, column].legend()
    drift_rows.append({"dt (fs)": label, "duration (ps)": t[-1], "slope (internal energy/ps)": slope,
                       "relative drift (1/ps)": slope/abs(energy.mean()),
                       "detrended RMS (internal energy)": np.std(residual)})
for ax in axes.flat:
    ax.grid(alpha=0.2)
plt.show()
drift = pd.DataFrame(drift_rows)
display(drift)
print("1 fs / 0.5 fs ratios: |slope| = "
      f"{abs(drift_rows[0]['slope (internal energy/ps)']/drift_rows[1]['slope (internal energy/ps)']):.3f}; "
      "detrended RMS = "
      f"{drift_rows[0]['detrended RMS (internal energy)']/drift_rows[1]['detrended RMS (internal energy)']:.3f}")

# Five contiguous blocks expose slow relaxation; block SEs are descriptive,
# not reliable equilibrium uncertainties when those block means keep drifting.
temperature_rows = []
fig, axes = plt.subplots(1, 3, figsize=(14, 3.8), constrained_layout=True)
for ax, (label, frame) in zip(axes, [("NVE 1 fs", nve1), ("NVE 0.5 fs", nvehalf), ("NVT, longer continuation", equil)]):
    blocks = np.array_split(np.arange(len(frame)), 5)
    times = [frame.time_internal.iloc[idx].mean()*time_unit_ps for idx in blocks]
    for field, kind in [("T_CH3_K", "CH3"), ("T_CH2_K", "CH2")]:
        ax.plot(times, [frame[field].iloc[idx].mean() for idx in blocks], "o-", label=kind)
        late = frame[field].to_numpy()[len(frame)//2:]
        means = np.array([part.mean() for part in np.array_split(late, 5)])
        temperature_rows.append({"run": label, "type": kind, "mean, last half (K)": late.mean(),
                                 "5-block SE (K)": means.std(ddof=1)/np.sqrt(5), "set point (K)": 293.0})
    ax.axhline(293, color="black", ls="--", label="293 K")
    ax.set(title=label, xlabel="Time within run (ps)", ylabel="Block mean temperature (K)")
    ax.legend()
    ax.grid(alpha=0.2)
plt.show()
display(pd.DataFrame(temperature_rows))
display(pd.DataFrame({"NVT time block (ps)": ["0-10", "10-20", "20-30", "30-40", "40-50"],
                      "mean potential energy (internal)": [part.mean() for part in np.array_split(equil.Epot_internal.to_numpy(), 5)]}))
print("Conclusion: good force/energy checks, but continued structural relaxation; equilibrium is not yet demonstrated.")
