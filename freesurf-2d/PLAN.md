# Plan: 2D two-phase VOF Navier–Stokes wave solver in C, with a Next.js viewer

## Context
The project directory `/Users/wgong/projects-code-base/personal/CFD-freesurface-flow` is empty. The user wants a C application that solves the incompressible Navier–Stokes equations for waves moving on a free liquid surface. They also want a Next.js front end that shows how the flow changes over time, with a time bar they can drag. Their choices:
- **2D** vertical (x–z) slice of a tank, with a **two-phase VOF** method (water and air, sharp PLIC interface)
- A **standing wave / sloshing** scenario, which can be checked against linear theory
- **VTK output for ParaView** (kept as an option)
- **Web app:** start runs from the browser (a Next.js API route runs the C solver), view finished runs and runs still in progress, and use a **field selector** (C, p, |u|, vorticity) with a colour legend and a time bar

Environment: Apple clang 21, make, Node 26 and npm 11 are installed. There is no cmake or OpenMP. The solver is plain C11 plus libm with no outside dependencies.

## Part A: C solver

### Equations
- Momentum: ρ(∂u/∂t + u·∇u) = −∇p + ∇·[μ(∇u+∇uᵀ)] + ρg
- Mass: ∇·u = 0
- Interface: ∂C/∂t + ∇·(uC) = 0
- Fluid properties: ρ = Cρ_l+(1−C)ρ_g and μ = Cμ_l+(1−C)μ_g. Surface tension is left out because it is negligible at tank scale.

### Numerics
- **Grid:** staggered MAC grid. p, C, ρ, and μ sit at cell centres, u/w on the faces, with 2 ghost layers. Double precision.
- **Each step:**
  1. Compute dt = min(CFL, viscous 0.25·h²·ρ_min/μ_max, gravity √(2h/g), dt_max), clipped so output times are hit exactly.
  2. **VOF advection** with the Weymouth–Yue (2010) operator-split PLIC scheme. It conserves mass to machine precision and keeps C between 0 and 1. The x/z sweep order alternates each step. Interface normals come from Youngs' method, and line↔volume conversions use the Scardovelli–Zaleski analytic formulas.
  3. Update ρ and μ. Face densities are arithmetic averages; corner μ is a 4-cell average.
  4. **Predictor u\*:** ENO2 upwind advection, the full variable-μ viscous term ∇·(2μD), and gravity added at the faces.
  5. **Pressure Poisson equation** ∇·(1/ρ_f ∇p) = ∇·u\*/dt. It is a 5-point SPD matrix with Neumann conditions at the walls. One air cell in a top corner has its pressure fixed so the system has a unique solution. It is solved with MIC(0)-preconditioned conjugate gradient, warm-started from the previous p.
  6. **Correction:** u = u\* − dt/ρ_f ∇p, then apply boundary conditions.
- This ordering behaves like symplectic Euler, so wave energy does not drift. A flat interface is in exact discrete hydrostatic balance.
- **Boundary conditions:** closed tank. Normal velocity is 0. Tangential velocity is free-slip by default, with `wall_bc=noslip` available. C has zero gradient at walls.

### CLI and config
- Usage: `build/freesurf <case.cfg> [key=value ...]`. Overrides are applied after the file is read, and unknown or invalid keys are errors.
- Default `cases/standing_wave.cfg`:
  - Tank Lx=1.0 m, Lz=0.5 m. Water depth h=0.25, amplitude A=0.01, mode n=1.
  - Fluids: ρ 1000/1.2, μ 1e-3/1.8e-5, g=9.81.
  - Grid 128×64, cfl=0.3, t_end=10, output_dt=0.025, `vtk_output=1`.
  - Analytic period: T = 2π/√(gk·tanh(kh)) ≈ 1.3975 s.

### Outputs (everything goes in `output_dir`)
| File | Purpose |
|---|---|
| `manifest.json` | Read by the web app. Contains: status (`running`/`done`), progress, params, nx, nz, Lx, Lz, analytic period, field names, global min/max per field, and the frame list `[{i,t}]`. Rewritten atomically after each output (write to a temp file, then rename). |
| `frames/frame_NNNNNN.bin` | One frame. 5 Float32 little-endian fields, cell-centred, row-major with the bottom row first: C, p, p_dyn (p minus the flat-surface hydrostatic profile, so the wave signal is visible), \|u\|, ω. Written as temp file + rename. |
| `vtk/*.vti`, `*.vtp`, `*.pvd` | ParaView output, only written when `vtk_output=1`. The `.vti` files (ImageData, appended raw Float32) hold C, p, ρ, velocity, and ω. The `.vtp` files hold the PLIC interface segments. Each `.pvd` series lists all outputs with their simulation times. |
| `probes.csv` | Columns: t, η_left, η_center, η_right, volume, KE, PE, divmax, poisson_iters |
| stdout | Progress lines. At the end, a summary: measured period (from zero crossings) vs analytic, volume drift, and amplitude ratio. |

### C files (all new)
```
Makefile                  make → build/freesurf; make debug (ASan+UBSan); make test; make run; make clean
cases/standing_wave.cfg
src/main.c                CLI + time loop
src/config.{c,h}          Params, defaults, parser, overrides, validation
src/sim.{c,h}             Sim struct, alloc/free, ghost-aware index macros
src/bc.{c,h}              ghost fill for u, w, C
src/vof.{c,h}             Youngs normals, PLIC α↔V, WY split advection, property update
src/momentum.{c,h}        ENO2 advection, viscous stress, gravity → u*
src/poisson.{c,h}         coefficients, MIC(0), PCG
src/init.{c,h}            standing wave: exact cell fractions by quadrature of η(x)=h+A·cos(kx)
src/diagnostics.{c,h}     probes, volume, energies, div, period measurement, summary
src/vtk_io.{c,h}          .vti/.vtp/.pvd writers
src/web_io.{c,h}          frame .bin + manifest.json writers (atomic), running field min/max
tests/test_plic.c         α(V)→V(α) round-trip ≤1e-12
tests/test_poisson.c      manufactured solution, ~2nd-order convergence, 1000:1 coefficient jump
tests/test_vof.c          solid-body rotation of a disk: mass ≤1e-12, 0≤C≤1
```

## Part B: Next.js front end (`web/`)
Scaffold with `npx create-next-app@latest web` (TypeScript, App Router, Tailwind, ESLint, npm). Before writing route handlers, read the installed version's own docs (`web/node_modules/next/dist/docs/`), for example on async `params` in route handlers. The only extra dependencies are those from create-next-app; colormaps and canvas rendering are written by hand.

### Server (Node runtime, `dynamic = 'force-dynamic'`)
- `web/lib/paths.ts` sets up these paths:
  - `ROOT` defaults to `..` and can be overridden with the `FREESURF_ROOT` env var.
  - `RUNS_DIR` is `ROOT/runs`, and `BIN` is `ROOT/build/freesurf`.
  - Run ids must match `^[a-z0-9-]+$` before any filesystem access, which prevents path traversal.
- `web/lib/params.ts` holds one shared schema used by both the form and the server:
  - Fields: nx (16–512), nz (8–256), depth, amplitude, mode (1–8), t_end (≤60), output_dt (≥0.005, at most 2000 frames), wall_bc, vtk_output.
  - Cross-checks: 0 < depth−A and depth+A < Lz.
- `web/lib/runner.ts` starts and tracks solver processes:
  - `spawn(BIN, [cfg, ...'key=value', 'output_dir=…'])` with an args array and no shell.
  - At most 2 runs at once, held in a module-level Map.
  - stdout/stderr go to `runs/<id>/log.txt`.
  - On a non-zero exit it writes `runs/<id>/status.json` as `{status:'failed', tail}`.
- API routes:
  - `GET /api/runs`: list runs, newest first, with id, status, progress, frame count, and key params. Merges `manifest.json` and `status.json`.
  - `POST /api/runs`: validate, create an id (timestamp + random suffix), spawn, return `{id}`. Returns 409 if the run limit is reached and 400 if validation fails.
  - `GET /api/runs/[id]`: the manifest, or 404.
  - `DELETE /api/runs/[id]`: kill the process if it is running, then remove the directory.
  - `GET /api/runs/[id]/frames/[n]`: send the frame bytes as `application/octet-stream`, with an immutable cache header.

### UI
- `app/layout.tsx` has a sidebar with:
  - a **New simulation** form (`components/RunForm.tsx`). It shows the analytic period live as you edit, and opens the new run once it is created.
  - a **Run list** (`components/RunList.tsx`) with status badges and progress bars. It polls every 1 s while any run is in progress and has a delete button.
- `app/page.tsx` shows an empty state. `app/runs/[id]/page.tsx` shows the viewer.
- `components/Viewer.tsx` contains:
  - **Field selector:** C, p, p_dyn, |u|, ω.
  - **Canvas heatmap** (`components/FieldCanvas.tsx`): an offscreen nx×nz ImageData built from a 256-entry colormap lookup table, scaled to the tank's aspect ratio with z pointing up, plus axis ticks in metres and a hover readout of (x, z, value).
  - **Colour legend** (`components/ColorLegend.tsx`). Its range is the run's global min/max for that field, so colours stay consistent over time. ω uses a diverging map symmetric about 0; the other fields use a sequential map.
  - **Time bar** (`components/TimeBar.tsx`): a range slider over the frame index, the current t / t_end, play/pause, ±1 frame buttons, a speed setting (0.25×–4×, where 1× is real time based on output_dt), ←/→ keys, and a "follow live" toggle for runs still in progress.
- `web/lib/frames.ts` handles frame loading:
  - It parses a frame `ArrayBuffer` into Float32Array views, one per field, after checking the byte length against nx·nz·nfields·4.
  - Frames are cached with an LRU limit of about 300 MB.
  - The frame being viewed is fetched first and stale requests are cancelled; other frames are prefetched in the background, 4 at a time.
  - "Loaded N/M" is shown.
- `web/lib/colormaps.ts` holds the viridis and blue-white-red lookup tables.

## Docs
`README.md` covers the equations and numerics, how to build and run the solver, the config keys, ParaView tips, and the web app (`make && cd web && npm install && npm run dev`, then open http://localhost:3000).

## Verification
1. `make && make test`: the PLIC, Poisson order, and VOF mass tests pass.
2. `make debug` plus a tiny run (`nx=32 nz=16 t_end=0.5`): sanitizers report no errors.
3. Still water (`amplitude=0 t_end=1`): max|u| should be about 1e-10 or less.
4. `make run`: period within about 2% of 1.3975 s, volume drift about 1e-12, and the period error should shrink at 256×128. Check that `manifest.json` is valid JSON and that every frame's size is 5·nx·nz·4 bytes.
5. `cd web && npm run lint && npm run build` both pass.
6. With `npm run dev`:
   - `curl -X POST /api/runs` with a small case should return an id.
   - Polling `GET /api/runs/<id>` should show progress and then `done`.
   - Fetching frame 0 should return the right number of bytes.
   - Invalid params should return 400, and an id like `../x` should return 400 or 404.
7. Manual browser check at http://localhost:3000: start a run, watch it fill in live, drag the time bar, play and pause, switch fields, and check the legend and hover readout. The surface should slosh like a cosine with a node at x=L/2.
