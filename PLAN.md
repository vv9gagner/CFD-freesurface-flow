# Plan: split the 2D solver into its own app, add a separate 3D solver, one shared front end

## Context
The project currently has one 2D two-phase VOF Navier–Stokes solver in C (`src/`, `tests/`, `cases/`, root `Makefile`) and a Next.js viewer (`web/`). That viewer starts runs, polls their progress, and has a time bar plus a field selector. The 2D solver is validated: at 128×64 the period error is +0.06%, volume drift is about 1e-13, and still water stays exactly at rest.

The user wants three things:
1. **The current 2D code moved into its own folder as a standalone app.**
2. **3D simulation added as a separate code base.** The two solvers share no C code.
3. **The same front end** for both, with this plan **saved as MD files** in the repo.

User choices:
- 3D viewer: **slice views** (x–z, y–z, x–y, with a position slider) **plus an interactive 3D free-surface view (three.js)**, all driven by the existing time bar.
- 3D scenario: **standing wave with modes (m, n)**.

User feedback on the first draft: **use GNU C**. Both solvers will be built with **GNU GCC** (`gcc-16` 16.2 from Homebrew is installed) in GNU C mode (`-std=gnu11`) instead of Apple clang.

Facts about the environment:
- The repo has `git init` but no commits, and all files are untracked, so moving them with plain `mv` loses nothing.
- `gcc-16 -fopenmp` compiles and links with libgomp, so OpenMP needs no extra library. The ASan and UBSan runtimes are installed.
- Node 26 has `Float16Array`.
- **Only 31 GB of disk is free**, so 3D frames must be compact.
- `runs/` holds two finished 2D runs from the user. They must keep working.

## Compiler (both apps)
- Each Makefile picks GNU GCC when `CC` is left at its default: the first of `gcc-16 gcc-15 gcc-14 gcc-13` found on PATH. If none is found, it stops with "GNU GCC not found: `brew install gcc`". `CC=…` can still override it.
- Flags: `-std=gnu11 -O2 -Wall -Wextra -Wpedantic -Wshadow`. The 3D release build adds `-fopenmp`.
- Debug: `-O0 -g -fsanitize=address,undefined`. If GCC's ASan does not run on this macOS/arm64 setup, the debug target falls back to UBSan only, and the README says so.
- Changing compiler alters floating-point rounding (FMA contraction), so 2D results are expected to match the clang build to about 4 significant digits, not bit for bit.

## Target layout
```
CFD-freesurface-flow/
├── Makefile              # delegates: all / test / clean / run2d / run3d
├── README.md             # overview of both apps + web
├── PLAN.md               # this plan
├── freesurf-2d/          # moved as-is: Makefile, README.md, PLAN.md (old plan), src/, tests/, cases/
├── freesurf-3d/          # new: Makefile, README.md, src/, tests/, cases/standing_wave_3d.cfg
├── web/                  # same Next.js app, now handles 2D and 3D runs
└── runs/                 # unchanged
```

## Step 0: Save the plans
- This plan is saved as `PLAN-3D.md` at the repo root (done while planning). The original 2D plan in `PLAN.md` is left untouched until Step 1.
- In Step 1, move `PLAN.md` to `freesurf-2d/PLAN.md` and rename `PLAN-3D.md` to `PLAN.md`.

## Step 1: Move the 2D solver into `freesurf-2d/` (mechanical)
- `mv src tests cases Makefile freesurf-2d/`. Delete the root `build/`, which will be rebuilt as `freesurf-2d/build/`.
- Move the solver-specific sections of the README (method, validation, CLI, output, ParaView) to `freesurf-2d/README.md`.
- Only code change: `src/output.c` `output_manifest()` writes `"dim": 2`. The web app treats manifests without `dim` as 2D, which keeps old runs working.
- The 2D `Makefile` switches to GNU GCC (see "Compiler" above).
- The new root `Makefile` runs `$(MAKE) -C freesurf-2d` and `-C freesurf-3d` for all, test and clean, and adds `run2d` and `run3d`.
- `.gitignore` ignores `*/build/`, `output/` and `runs/`.

## Step 2: 3D solver `freesurf-3d/` (new code base, GNU C11 + OpenMP)
It mirrors the 2D file structure: `util`, `config`, `sim`, `bc`, `vof`, `momentum`, `poisson`, `init`, `solver`, `diagnostics`, `output`, `vtk_io`, `main`, plus `tests/` and `cases/standing_wave_3d.cfg`. The binary is `build/freesurf3d`. Code is copied from the 2D files and generalised; nothing is shared.

| Part | 3D design |
|---|---|
| Grid | Staggered MAC grid with u, v, w on faces and cell-centred p, C, ρ, μ. `IX(s,i,j,k)` uses one shared array layout `(nx+1+2NG)(ny+1+2NG)(nz+1+2NG)`. z is vertical (gravity −z) and y is the new horizontal axis |
| BCs | 6 closed walls: normal velocity 0, slip or no-slip tangential velocity, zero-gradient C |
| PLIC | `plane_volume(n, α)`: analytic plane/cube volume (Scardovelli & Zaleski). `plane_alpha(c, n)`: safeguarded Newton with bisection fallback, since volume is monotonic in α. Normals: 3D Youngs, with a centred-columns (height-function) candidate along the dominant axis chosen by the MYC rule |
| VOF advection | Weymouth & Yue with 3 sweeps per step and a cyclic order (xyz, yzx, zxy). Slab fluxes use the same scaled-normal trick as in 2D |
| Momentum | ENO2 advection of u, v, w; full ∇·(2μD) with edge-averaged μ (4 cells) for the shear terms; **reduced gravity** [ρ]g(z_int−h)∇C/ρ_f, where z_int is the PLIC plane height at the cell's column centre |
| Pressure | 7-point variable-coefficient Poisson; 3D MIC(0)-PCG; p pinned in one gas cell. The matvec, dot products and updates are OpenMP-parallel; the preconditioner stays sequential |
| Time step | Same limits as 2D (CFL ≤ 0.5, explicit viscosity, gravity, dt_max), applied to all 3 directions |
| OpenMP | `#pragma omp parallel for collapse(2)` on the stencil loops. GCC `-fopenmp` with libgomp; `OPENMP=0` gives a serial build. Debug builds are serial |
| Init | η = h + A cos(mπx/Lx) cos(nπy/Ly); cell fractions from 16×16 sub-samples in (x, y) of the exact column fill |
| Diagnostics | η(x,y) as column sums; probes at corner (0,0), centre and opposite corner; volume; KE and PE (column-centroid formula); max \|∇·u\|; period from corner-probe zero crossings; wave-energy ratio. Analytic period uses the two-layer formula with k = π√((m/Lx)² + (n/Ly)²) |

### 3D config keys
- `Ly`, `ny`, `mode_x`, `mode_y`.
- Everything else matches 2D.
- Defaults:
  - tank Lx=1, Ly=0.5, Lz=0.5, grid 64×32×32 (cubic cells),
  - depth h=0.25, amplitude A=0.01, modes (1,1),
  - t_end=10, output_dt=0.05, dt_max=0.005,
  - `vtk_output=0`.

### 3D output (disk-aware)
| File | Content |
|---|---|
| `manifest.json` | Format `"freesurf3d-web-1"`, `dim: 3`. Also: nx, ny, nz, Lx, Ly, Lz; stored fields `["C","p_dyn","u","v","w"]` (float16); `p_hydro[nz]`; display fields with global ranges (C, p, p_dyn, speed, vorticity ±max\|ω\|); surface η range; frame_times; byte sizes; summary. Written atomically |
| `frames/frame_N.bin` | 5 × nx·ny·nz **float16** values, i fastest then j then k, cell-centred. The default grid gives 655 KB per frame and about 130 MB per run |
| `surface/surface_N.bin` | nx·ny float32 η(x,y) |
| `vtk/` (optional) | 3D `.vti` files (C, p, velocity, vorticity) and a `.pvd` series |

## Step 3: Web app (same front end)
- **`lib/paths.ts`**: `SOLVERS = {2: {bin, case}, 3: {bin, case}}`, pointing at `freesurf-2d/build/freesurf` and `freesurf-3d/build/freesurf3d`.
- **`lib/params.ts`**:
  - Add `dim: 2 | 3`, `Ly`, `ny`, `mode_y`, with defaults and limits per dimension. 3D limits: up to 128×128×64 cells and 512k cells in total, at most 600 frames.
  - `analyticPeriod` handles (m, n).
  - Size and cost estimates per dimension (3D frames are float16).
  - `solverArgs` maps `mode` and `mode_y` to the 3D keys `mode_x` and `mode_y`.
- **`lib/types.ts`**: `Manifest` becomes a 2D/3D union, with a missing `dim` read as 2.
- **`lib/runner.ts`**:
  - Picks the solver by `dim`.
  - At most 2 runs at once, and at most 1 of them 3D, since a 3D run uses every core.
  - `run.json` stores `dim`.
- **Frame API** (`app/api/runs/[id]/frames/[n]/route.ts`):
  - 2D behaves as before.
  - 3D adds `?view=xz|yz|xy&at=<k>`. The server reads the float16 frame and **returns a slice in the same 5-field float32 layout as a 2D frame**: C, p = p_dyn + p_hydro, p_dyn, speed, and ω normal to the slice from in-plane central differences. The existing viewer pipeline (FrameStore, FieldCanvas, ranges, CSV) therefore works unchanged on slices.
  - `?view=surface` returns η(x,y).
  - A bad `view` or `at` returns 400.
  - The extraction code lives in `lib/slices.ts`: relative imports only and pure functions, so it can be tested.
- **`lib/frames.ts`**: `FrameStore` takes a URL builder and an expected byte size, so slices and surfaces get their own cache.
- **`components/Viewer.tsx`**:
  - 2D runs work as before.
  - 3D runs get view tabs: **3D surface | Slice x–z | Slice y–z | Slice x–y**.
  - A slice-position slider shows the index and its coordinate. Vertical slices default to the mid-plane; the horizontal slice defaults to the cell just below the mean surface.
  - Setup and results panels show Ly, ny and the mode (m, n).
- **`components/FieldCanvas.tsx`**: axis names become props (x/z, y/z, x/y), and the surface line appears only on vertical slices.
- **`components/Surface3D.tsx`** (new, client-only via `next/dynamic` with `ssr: false`; the Next 16 lazy-loading docs are checked first):
  - three.js mesh of nx×ny vertices, displaced by η−h with a vertical-exaggeration slider.
  - Coloured with the existing diverging blue–red scale from `lib/colormaps.ts`.
  - Tank wireframe box, OrbitControls, theme-aware background.
  - Resizes with the container and frees GPU resources on unmount.
  - Dependency: `three` with `@types/three`.
- **`components/RunForm.tsx`**:
  - A 2D/3D segmented toggle at the top. In 3D mode it also shows Ly, ny and mode_y.
  - Disk and time estimates for 3D.
- **`components/RunList.tsx`**: a 2D/3D tag, the grid `64×32×32`, and the modes `(m,n)`.
- **`npm test`**: `node --test lib/*.test.ts` (Node 26 strips TypeScript types natively). `slices.test.ts` tests slice shapes and values on synthetic data, and checks that the vorticity of a solid-body rotation equals 2Ω.

## Verification
1. **2D after the move (now built with gcc-16):**
   - The build has no warnings and `make -C freesurf-2d test` passes.
   - The default run gives the same summary as the clang build to about 4 significant digits (period error +0.06%, energy ratio 0.877, volume drift about 1e-13).
   - Still water stays exactly at rest.
   - The two existing runs still open in the viewer.
2. **3D unit tests** (`make -C freesurf-3d test`):
   - `plane_volume` against numerical integration of the 2D area formula over z (random planes, tolerance 1e-8).
   - α↔V round trip and slab-split conservation (1e-12).
   - Poisson: second-order convergence, and an exact solve across a 1000:1 jump.
   - Solid-body rotation of a sphere: volume ≤ 1e-12, 0 ≤ C ≤ 1, shape error < 5%.
3. **3D physics:**
   - Still water stays exactly at rest.
   - **Cross-check with 2D:** mode (1,0) with ny=4 must match the 2D period at the same nx, nz to within 0.1%.
   - Default mode (1,1): period within about 1% of linear theory, volume drift about 1e-12.
   - Serial and OpenMP builds agree; record the speedup and wall time in the README.
4. **Sanitizers:** `make -C freesurf-3d debug` plus a tiny run with ASan/UBSan reports nothing.
5. **Web:**
   - `npm run lint`, `tsc`, `npm run build`, and `npm test` pass.
   - curl:
     - create a 3D run and poll it until done,
     - fetch each slice type and the surface, checking byte lengths,
     - bad view/at returns 400,
     - a second concurrent 3D run returns 409,
     - old 2D runs still load.
6. **UI** (headless Chrome over DevTools, with swiftshader WebGL):
   - 2D page unchanged.
   - 3D page: switch views, move the slice slider, the 3D surface renders and rotates, the time bar scrubs and plays.
   - Screenshots in light and dark mode.

## Implementation notes (done, 2026-10-26)
All steps were implemented and verified. Where the result differs from the plan above:
- **Pressure preconditioner (3D).** Profiling showed the sequential MIC(0) solves taking 57% of the run time, so OpenMP alone gave almost no speedup. I replaced it with a **block-Jacobi MIC(0)**: one slab per thread along the longest axis, plain MIC(0) when serial. The CG vector operations are fused into fewer parallel loops. 4 threads were the fastest (1.9× over serial). Beyond that, the extra CG iterations and libgomp's kernel-level barriers on macOS cost more than they save, so `threads=0` uses at most 4 performance cores.
- **`plane_volume`.** It integrates the exact 2D section area (Simpson's rule between the kinks, which is exact) instead of using the closed-form cubic. The cubic loses precision when one normal component is tiny, as in quasi-2D waves.
- **Results.**
  - 3D mode (1,1): period +0.26%, volume drift 3e-13, 86 s for the default case.
  - 3D mode (1,0) against the 2D app at the same grid: periods agree to 6e-6.
  - Both solvers build warning-free with gcc-16 (`-std=gnu11`) and pass their unit tests and sanitizer runs.
- **Web.** Lint, typecheck, `npm test` and `next build` pass. The UI was driven in headless Chrome (2D page, 3D surface/slices, light and dark) with no console errors.
