# Freesurf: free-surface wave simulation in 2D and 3D

The project solves the incompressible Navier–Stokes equations for water and air together, tracking the free surface with a geometric volume-of-fluid method. It has two independent solver apps in C and one Next.js front end that runs and visualises both.

```
CFD-freesurface-flow/
├── freesurf-2d/   2D solver (vertical x–z slice)            → freesurf-2d/README.md
├── freesurf-3d/   3D solver (rectangular tank, OpenMP)      → freesurf-3d/README.md
├── web/           Next.js app: start runs, scrub results    → web/README.md
├── runs/          runs started from the web app
├── Makefile       builds and tests both solvers
└── PLAN.md        plan for the 2D/3D split (the original 2D plan is freesurf-2d/PLAN.md)
```

## Quick start

```
make && make test                       # both solvers, built with GNU GCC (brew install gcc)
cd web && npm install && npm run dev    # http://localhost:3000
```

In the browser, pick **2D** or **3D**, set the tank, grid and wave, and press **Run**. The run opens immediately and fills in as the solver works.
- **Time bar:** drag to scrub, press play (0.25×–4× real time), or step with ← →.
- **Fields:** liquid fraction, pressure, dynamic pressure, speed and vorticity.
- **3D runs** open on an interactive surface view (rotate and zoom). They also have x–z, y–z and x–y slices with a position slider.

From the command line:
```
make run2d      # freesurf-2d/build/freesurf freesurf-2d/cases/standing_wave.cfg
make run3d      # freesurf-3d/build/freesurf3d freesurf-3d/cases/standing_wave_3d.cfg
```

## The two solvers

| | freesurf-2d | freesurf-3d |
|---|---|---|
| Domain | vertical slice x–z | tank x–y–z |
| Initial surface | h + A cos(mπx/L) | h + A cos(mπx/Lx) cos(nπy/Ly) |
| Interface | PLIC (line), MYC normals | PLIC (plane), 3D MYC normals |
| Pressure solver | MIC(0)-PCG | block-Jacobi MIC(0)-PCG, OpenMP |
| Frame storage | float32 | float16 + float32 surface map |
| Default run | 128×64, 10 s, ~20 s | 64×32×32, 10 s, ~1.5 min |
| Period vs linear theory | +0.06 % | +0.26 % (mode 1,1) |

The two code bases agree to 6×10⁻⁶ in period on the same 2D wave (3D mode (1,0)). See each app's README for the numerics, validation, options and output formats.

## Web app

| Endpoint | |
|---|---|
| `GET /api/runs` | list runs (2D and 3D) |
| `POST /api/runs` | start a run; the JSON body holds `dim` (2 or 3) and the parameters of `web/lib/params.ts` |
| `GET /api/runs/:id`, `DELETE /api/runs/:id` | status and manifest; stop and delete |
| `GET /api/runs/:id/frames/:n` | 2D frame |
| `GET /api/runs/:id/frames/:n?view=surface` | 3D surface height map |
| `GET /api/runs/:id/frames/:n?view=xz\|yz\|xy&at=k` | 3D slice, in the 2D frame layout |

- At most 2 solver processes run at once, of which at most 1 is 3D.
- Solvers are started with an argument array (no shell), and run ids are validated before the filesystem is touched.

## Limitations

- No surface tension or contact-angle model.
- Explicit viscosity, so fine grids take small time steps.
- Momentum advection is not mass-consistent. This is fine for sloshing, but strongly breaking waves would need a momentum-conserving VOF.
- One scenario per solver (standing waves). New cases go in each app's `src/init.c`.
