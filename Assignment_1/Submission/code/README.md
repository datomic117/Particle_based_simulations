# Building and running

## Build

From the `code/` directory:

```
gcc -O2 -o main main.c io.c forces.c integrate.c -lm
```

This produces `main` (or `main.exe` on Windows).

## Run

```
./main [--integrator euler|verlet] [--dt SECONDS] [--steps N]
       [--trajectory-every K] [--diagnostics-every K]
       [--bodies NAME1,NAME2,...]
```

Defaults: `--integrator verlet --dt 3600 --steps 10 --trajectory-every 2
--diagnostics-every 1`, all bodies in `../data/bodies_2026-09-01.dat`.

- `--integrator euler|verlet` selects Euler-forward or velocity-Verlet.
- `--dt` sets the time step in seconds.
- `--steps` sets the number of steps to integrate.
- `--trajectory-every` writes a trajectory frame every K-th step.
- `--diagnostics-every` writes a diagnostics row every K-th step. This
  only thins out the CSV; forces and accelerations are still computed
  every step regardless, so it does not affect the simulated physics.
  Useful for long runs (C4, C5) where writing every step would produce
  an unwieldy `diagnostics.csv`.
- `--bodies` restricts the run to a comma-separated subset of body names
  (case-insensitive, matched against the `name` column of the input
  file), e.g. `--bodies Sun,Earth,Moon`. Omit it to use every body in
  the file.

The input file is currently hard-coded in `main.c` as
`../data/bodies_2026-09-01.dat`.

## Output

Two CSV files are written to `../data/`:

- `diagnostics.csv`: one row every `--diagnostics-every` steps (plus
  step 0) with `step,t,K,U,E,Lx,Ly,Lz,Rcom_x,Rcom_y,Rcom_z`.
- `trajectory.csv`: one row per body per saved frame with
  `step,t,body,x,y,z`.

## Files

- `vec3d.h`: the `Vec3D` type and inline vector helpers (B2).
- `io.c` / `io.h`: reading the Horizons snapshot files and writing the
  diagnostics/trajectory CSVs (B3, B6, B7).
- `forces.c` / `forces.h`: pairwise gravitational accelerations and
  total potential energy (B4).
- `integrate.c` / `integrate.h`: Euler-forward and velocity-Verlet
  steppers, selected at runtime (B5).
- `main.c`: command-line parsing, allocation, and the main time loop.
