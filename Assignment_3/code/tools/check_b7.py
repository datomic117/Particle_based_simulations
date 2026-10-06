import pandas as pd
import numpy as np

filename = "../../data/b7_nve.csv"

try:
    df = pd.read_csv(filename)
except UnicodeError:
    df = pd.read_csv(filename, encoding="utf-16")

# 1 internal time unit = 1.0967 ps
time_ps = df["time"] * 1.0967

E = df["Etot"].to_numpy()
T = df["T"].to_numpy()

mean_E = np.mean(E)
std_E = np.std(E, ddof=1)

E_range = np.max(E) - np.min(E)
relative_range = E_range / abs(mean_E)

slope, intercept = np.polyfit(time_ps, E, 1)
relative_drift_per_ps = slope / abs(mean_E)

print("Number of samples:", len(df))
print("Simulation time:", time_ps.iloc[-1], "ps")
print()

print("TOTAL ENERGY")
print("Mean Etot =", mean_E)
print("Std Etot =", std_E)
print("Minimum Etot =", np.min(E))
print("Maximum Etot =", np.max(E))
print("Relative full range =", relative_range)
print("Linear slope =", slope, "energy units / ps")
print("Relative drift =", relative_drift_per_ps, "per ps")

print()
print("TEMPERATURE")
print("First saved T =", T[0], "K")
print("Final T =", T[-1], "K")
print("Mean T first 100 samples =", np.mean(T[:100]), "K")
print("Mean T last 100 samples =", np.mean(T[-100:]), "K")