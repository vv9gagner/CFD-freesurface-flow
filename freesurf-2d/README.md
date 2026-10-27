# freesurf-2d: 2D two-phase Navier–Stokes solver for free-surface waves

This is the 2D app of the project. It solves water and air together in a vertical (x–z) slice of a closed tank and tracks the free surface with a geometric volume-of-fluid (VOF) method. The code is GNU C11 with libm only, built with Homebrew GCC (`gcc-16`).

```
make            # build/freesurf (GNU GCC; override with make CC=...)
make test       # PLIC, Poisson and VOF unit tests
make run        # standing wave, 128x64, 10 s  (~20 s)
make debug      # build/freesurf_debug with AddressSanitizer + UBSan
```

## Physics

```
rho (du/dt + u.grad u) = -grad p + div( mu (grad u + grad u^T) ) + rho g
div u = 0
dC/dt + div(u C) = 0                    C = liquid volume fraction
rho = C rho_l + (1-C) rho_g,   mu = C mu_l + (1-C) mu_g
```

The default fluids are water and air, with g = 9.81 m/s². Surface tension is ignored. The tank is closed, with free-slip or no-slip walls.

**Case: standing wave.** The surface starts as `eta(x) = h + A cos(m pi x / L)`, at rest. Linear two-layer theory gives
`omega^2 = g k (rho_l - rho_g) / (rho_l coth(k h) + rho_g coth(k (Lz - h)))` with `k = m pi / L`.

## Numerical method

| Part | Method |
|---|---|
| Grid | Staggered MAC grid with 2 ghost layers |
| Interface | PLIC reconstruction with Mixed-Youngs-Centred normals and analytic line/area formulas |
| VOF advection | Weymouth & Yue (2010) split scheme. It conserves mass to round-off and keeps 0 ≤ C ≤ 1 for CFL ≤ ½ |
| Momentum | ENO2 upwind advection and the full variable-viscosity stress ∇·(2μD), both explicit |
| Gravity | Balanced-force "reduced gravity": the solver uses p_rgh = p + ρg(z−h), and gravity acts as an interfacial force [ρ]g(z_int−h)∇C/ρ at the PLIC interface height. This removes spurious air jets at a tilted interface |
| Pressure | Variable-coefficient Poisson equation, solved by conjugate gradient with a MIC(0) preconditioner |
| Time step | Advect C, update ρ and μ, predict u\*, then project (a symplectic-Euler splitting) |

## Validation

`make test`:
- **PLIC:** volume ↔ line round trip exact to 7×10⁻¹⁶.
- **Poisson:** second-order convergence, and an exact solve across a 1000:1 density jump to 6×10⁻¹¹.
- **VOF:** solid-body rotation conserves volume to 1×10⁻¹⁵, with a shape error of 0.3%.

Standing wave (h = 0.25 m, A = 0.01 m, L = 1 m, mode 1, 10 s, analytic T = 1.3993 s):

| Grid | Measured period | Error | Volume drift | Wave energy kept after 7 periods | Run time |
|---|---|---|---|---|---|
| 64×32   | 1.4003 s | +0.075 % | 5×10⁻¹³ | 71 % | 3 s |
| 128×64  | 1.4002 s | +0.063 % | 5×10⁻¹³ | 88 % | 20 s |
| 256×128 | 1.3996 s | +0.025 % | 9×10⁻¹⁴ | 97 % | 5.5 min |

Still water (`amplitude=0`) stays exactly at rest. The energy loss is numerical damping, and it drops by about 4× per refinement. The 3D app reproduces the 64×32 row to 6×10⁻⁶ when run as mode (1,0) (see `../freesurf-3d/README.md`).

## Command line

```
build/freesurf <case.cfg> [key=value ...]
build/freesurf cases/standing_wave.cfg nx=256 nz=128 amplitude=0.03 output_dir=output/big
```

| Key | Default | Meaning |
|---|---|---|
| `Lx`, `Lz` | 1.0, 0.5 | tank length and height [m] |
| `nx`, `nz` | 128, 64 | cells (keep dx ≈ dz) |
| `depth`, `amplitude`, `mode` | 0.25, 0.01, 1 | initial surface h + A cos(mπx/L) |
| `rho_liquid`, `rho_gas`, `mu_liquid`, `mu_gas`, `gravity` | water/air | fluid properties |
| `wall_bc` | `slip` | `slip` or `noslip` |
| `t_end`, `cfl`, `dt_max`, `output_dt` | 10, 0.3, 0.005, 0.025 | time stepping [s] |
| `poisson_tol`, `poisson_maxit` | 1e-9, 5000 | pressure solver |
| `output_dir`, `vtk_output`, `print_every` | output/standing_wave, 1, 100 | output |

Exit codes: 0 done, 1 bad input, 2 diverged or output failure, 3 cancelled (SIGINT/SIGTERM).

## Output (in `output_dir`)

| File | Content |
|---|---|
| `manifest.json` | `"format": "freesurf-web-1", "dim": 2`. Contains metadata, parameters, global field ranges, frame times and the summary. Written atomically |
| `frames/frame_NNNNNN.bin` | 5 float32 fields (C, p, p_dyn, \|u\|, ω), each nx·nz values, x fastest, bottom row first |
| `probes.csv` | Surface heights at the left wall, centre and right wall; volume; KE and PE; max \|∇·u\|; Poisson iterations |
| `vtk/fields.pvd`, `vtk/interface.pvd` | ParaView time series (cell fields; PLIC interface segments) |

**ParaView:** open both `.pvd` files and use the 2D view (the tank's z is VTK's y). Colour by `C`, or add *Contour* at C = 0.5.

## Source layout

```
src/main.c        CLI and time loop          src/solver.c       time step, dt limits
src/config.c      parameters and parser      src/sim.c          grid and arrays
src/vof.c         PLIC and VOF advection     src/momentum.c     predictor and reduced gravity
src/poisson.c     MIC(0)-PCG and projection  src/bc.c           wall boundary conditions
src/init.c        initial conditions         src/diagnostics.c  probes, energy, period
src/output.c      web frames and manifest    src/vtk_io.c       ParaView files
tests/            unit tests                 cases/             case files
```
