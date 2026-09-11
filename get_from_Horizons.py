#!/usr/bin/env python3
"""Fetch a solar-system snapshot from JPL Horizons (6EMA02, assignment A1).

Writes one line per body: id, name, mass (kg), position (m), velocity (m/s),
solar-system barycentric, ecliptic J2000 — the format of the provided
bodies_<epoch>.dat files. Masses are derived from Horizons' GM values as
m = GM / G with G = 6.6743015e-11; use the same G in your simulation and
your G*m products reproduce Horizons' GM exactly. Bodies without a GM in
Horizons (e.g. the asteroid Apophis) get mass 0.

Usage:
    python3 get_from_Horizons.py 2026-09-01 bodies_2026-09-01.dat
    python3 get_from_Horizons.py 2029-04-13 encounter.dat 399 301 "99942;"

Bodies are Horizons COMMAND codes: 10 = Sun, 199/299/399/499 = planet
centers, 301 = Moon, single digits 5..9 = planetary-system barycenters
(used for the moon-bearing outer planets, with their system masses — see
SYSTEM_GM below); a number with a trailing ';' looks up a small body by
its designation (99942; = Apophis). Default: the 12 bodies of the course
hand-out. Needs the `requests` package and internet access.
"""
import sys
import re
import datetime
import requests

G = 6.6743015e-11
API = "https://ssd.jpl.nasa.gov/api/horizons.api"
DEFAULT_BODIES = ["10", "199", "299", "399", "301", "499",
                  "5", "6", "7", "8", "9", "99942;"]

# For the moon-bearing outer planets the snapshot uses the planetary-system
# BARYCENTER (COMMAND 5..9), not the planet center: a model that contains no
# moons propagates exactly that point. Horizons returns no GM for barycenters,
# so their system GM (planet + moons, km^3/s^2) comes from the DE440/DE441
# constants (Park et al. 2021, AJ 161:105) — the same ephemeris Horizons uses.
SYSTEM_GM = {
    "5": 126712764.10,      # Jupiter system
    "6": 37940584.8418,     # Saturn system
    "7": 5794556.40,        # Uranus system
    "8": 6836527.10058,     # Neptune system
    "9": 975.50,            # Pluto system
}


def query_horizons(command, start, stop):
    """One Horizons vector query; returns the 'result' text."""
    params = {
        "format": "json", "EPHEM_TYPE": "VECTORS", "CENTER": "@0",
        "OUT_UNITS": "KM-S", "REF_PLANE": "ECLIPTIC", "REF_SYSTEM": "J2000",
        "VEC_CORR": "NONE", "VEC_TABLE": "2", "VEC_LABELS": "YES",
        "OBJ_DATA": "YES", "STEP_SIZE": "2d",
        "COMMAND": f"'{command}'",
        "START_TIME": f"'{start}'", "STOP_TIME": f"'{stop}'",
    }
    r = requests.get(API, params=params, timeout=60)
    r.raise_for_status()
    data = r.json()
    if "result" not in data:
        raise RuntimeError(f"Horizons returned no result for {command}: {data}")
    return data["result"]


def parse_name(txt):
    m = re.search(r"Target body name:\s*([^(\n]+)", txt)
    if not m:
        return "unknown"
    name = m.group(1).strip()
    if name.endswith(" Barycenter"):            # 'Jupiter Barycenter' -> Jupiter
        name = name[:-len(" Barycenter")]
    return name.split()[-1]


def parse_mass(txt):
    """Mass in kg from the object-data header, m = GM/G. GM lines look like
    'GM (km^3/s^2) = 398600.435436' or 'GM, 10^11 km^3/s^2 = 1.32712440018'.
    Returns 0.0 when Horizons lists no GM (test particles like Apophis)."""
    for line in txt.splitlines():
        # e.g. 'GM, km^3/s^2 = 4902.800066', 'GM (km^3/s^2) = 398600.435436',
        # the Sun's 'GM, 10^11 km^3/s^2 = 1.32712440018'. The unit must follow
        # GM directly ('GM 1-sigma' lines do not match), and the value is read
        # from the '=' belonging to GM — some lines hold two columns.
        m = re.search(r"GM[,(]?\s*(?:10\^(\d+)\s*)?\(?km\^3/s\^2\)?\s*[=~]\s*([\d.]+)",
                      line)
        if m:
            value = float(m.group(2))
            if m.group(1):                        # the Sun's '10^11' prefix
                value *= 10.0 ** int(m.group(1))
            return value * 1e9 / G                # km^3/s^2 -> m^3/s^2 -> kg
    m = re.search(r"Mass,?\s*x?\s*10\^(\d+)\s*\(?kg\)?\s*[=~]\s*([\d.]+)", txt)
    if m:
        return float(m.group(2)) * 10.0 ** int(m.group(1))
    return 0.0


def parse_state(txt):
    """First state vector between $$SOE and $$EOE, converted to m and m/s."""
    m = re.search(r"\$\$SOE(.*?)\$\$EOE", txt, re.S)
    if not m:
        raise RuntimeError("no $$SOE/$$EOE block in the Horizons result")
    numbers = re.findall(r"[XYZ][XYZ ]?=\s*(-?[\d.]+E[+-]\d+)", m.group(1))
    if len(numbers) < 6:
        raise RuntimeError("could not parse the state vector")
    x, y, z, vx, vy, vz = (float(v) * 1e3 for v in numbers[:6])
    return x, y, z, vx, vy, vz


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    epoch, outfile = argv[1], argv[2]
    bodies = argv[3:] if len(argv) > 3 else DEFAULT_BODIES
    start = datetime.date.fromisoformat(epoch)
    stop = start + datetime.timedelta(days=1)

    with open(outfile, "w") as f:
        print(f"# 6EMA02 solar-system snapshot, epoch {epoch} 00:00 (TDB), "
              "solar-system barycentric, ecliptic J2000; outer planets = "
              "planetary-system barycenters with system masses", file=f)
        print("# id  name  mass_kg  x_m  y_m  z_m  vx_m/s  vy_m/s  vz_m/s   "
              f"(masses = GM/G with G = {G:.7e})", file=f)
        for body in bodies:
            txt = query_horizons(body, start, stop)
            name = parse_name(txt)
            if body in SYSTEM_GM:
                mass = SYSTEM_GM[body] * 1e9 / G
            else:
                mass = parse_mass(txt)
            x, y, z, vx, vy, vz = parse_state(txt)
            ident = int(body.rstrip(";"))
            print(f"{ident:8d} {name:>9s} " +
                  " ".join(f"{v:20.12e}" for v in (mass, x, y, z, vx, vy, vz)),
                  file=f)
            print(f"fetched {name} (mass {mass:.6e} kg)")
    print(f"wrote {outfile}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
