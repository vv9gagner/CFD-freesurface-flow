# freesurf-3d: 3D two-phase Navier–Stokes solver for free-surface waves

This is the 3D app of the project. It solves water and air together in a closed rectangular tank and tracks the free surface with a geometric volume-of-fluid (VOF) method. It is a separate code base from `../freesurf-2d`: the same design, written independently for three dimensions. The code is GNU C11 with OpenMP (libgomp), built with Homebrew GCC (`gcc-16`).

```
make            # build/freesurf3d (OpenMP; make OPENMP=0 for a serial build)
make test       # 3D PLIC, Poisson and VOF unit tests
make run        # standing wave mode (1,1), 64x32x32, 10 s  (~1.5 min)
make debug      # serial build/freesurf3d_debug with AddressSanitizer + UBSan
```

## Physics

These are the same equations as in 2D, with velocity (u, v, w). z is vertical and y is the second horizontal axis.

**Case: standing wave with modes (m, n).** The surface starts at rest as
`eta(x, y) = h + A cos(m pi x / Lx) cos(n pi y / Ly)`. Linear two-layer theory gives the same dispersion relation as in 2D, with
`k = pi sqrt((m/Lx)^2 + (n/Ly)^2)`. Mode (m, 0) is a 2D wave, which allows a cross-check against the 2D app.

## Numerical method

| Part | 3D design |
|---|---|
| Grid | Staggered MAC grid: u, v, w on faces and p, C, ρ, μ at cell centres; 2 ghost layers |
| PLIC | `plane_volume` integrates the exact 2D section area over the smallest normal component. This is exact, because the section area is piecewise quadratic and Simpson's rule is applied between its kinks. It also stays accurate for nearly 2D normals, where the closed-form cubic loses precision. `plane_alpha` uses Newton's method with a bisection safeguard |
| Normals | 3D Mixed-Youngs-Centred: centred-column (height-function) candidates along x, y and z, plus Youngs, with the MYC selection rule. Planar interfaces are reconstructed exactly |
| VOF advection | Weymouth & Yue (2010) with three split sweeps per step, in a cyclic order (xyz, yzx, zxy) |
| Momentum | ENO2 advection of u, v, w; full ∇·(2μD) with viscosity on cell edges; explicit |
| Gravity | Balanced-force reduced gravity [ρ]g(z_int−h)∇C/ρ at the PLIC plane height |
| Pressure | 7-point variable-coefficient Poisson equation. PCG with a block-Jacobi MIC(0) preconditioner: one slab per thread along the longest axis, with plain MIC(0) when serial. Vector operations are fused into a few OpenMP loops per iteration |

### Performance

The test case is the default case on an Apple M-series chip (6 performance and 4 efficiency cores), 64×32×32 cells for 2000 steps. About 92% of the time goes to the pressure solve.

| Threads | 1 | 2 | 4 | 6 | 10 |
|---|---|---|---|---|---|
| Wall time per simulated second | 15.9 s | 10.2 s | **8.6 s** | 9.2 s | 15.1 s |
| CG iterations per step | 69 | 90 | 114 | 131 | 156 |

More threads mean more preconditioner blocks, so CG needs more iterations, and libgomp barriers on macOS go through the kernel. `threads=0`, the default, therefore uses at most 4 performance cores. A multigrid preconditioner would be the next step for larger grids.

## Validation

`make test`:
- **PLIC:** `plane_volume` matches independent numerical integration to 2×10⁻¹³. The α ↔ V round trip and slab splitting are exact to 1×10⁻¹⁵, including degenerate (2D-like) normals. Exact normals are reconstructed for planar interfaces.
- **Poisson:** second-order convergence (error ratio 2.00 at 16³ → 32³ → 64³), and an exact solve across a 1000:1 density jump to 3×10⁻¹¹.
- **VOF:** a sphere makes one revolution about y, then one about x, on a 40³ grid. Volume drift is 3×10⁻¹⁵, C stays within [0, 1], and the shape error is 5%.

Physics checks:

| Case | Measured period | Analytic | Error | Wave energy kept | Volume drift | Run time |
|---|---|---|---|---|---|---|
| Mode (1,1), 64×32×32, 10 s | 0.7828 s | 0.7807 s | +0.26 % | 46 % after 13 periods | 3×10⁻¹³ | 86 s |
| Mode (1,0), 64×4×32, 10 s | 1.400344 s | 1.399290 s | +0.075 % | 70.6 % | 7×10⁻¹⁴ | 26 s |
| 2D app, 64×32, 10 s (reference) | 1.400335 s | 1.399290 s | +0.075 % | 70.6 % | 5×10⁻¹³ | 3 s |

- **Cross-check with 2D.** The two independent code bases agree to 6×10⁻⁶ in period and 3×10⁻⁴ in energy for the same 2D wave.
- **Still water.** With `amplitude=0` the tank stays exactly at rest: kinetic energy 0, max |∇·u| 0, volume drift 0.
- **Energy loss in mode (1,1).** The wavelength is only about 57 cells, so numerical damping is stronger than in the 2D runs. Refining the grid reduces it, as the 2D convergence table shows.

## Command line

```
build/freesurf3d <case.cfg> [key=value ...]
build/freesurf3d cases/standing_wave_3d.cfg mode_y=0 ny=4 Ly=0.0625     # quasi-2D cross-check
```

| Key | Default | Meaning |
|---|---|---|
| `Lx`, `Ly`, `Lz` | 1.0, 0.5, 0.5 | tank size [m] (z vertical) |
| `nx`, `ny`, `nz` | 64, 32, 32 | cells (keep them roughly cubic) |
| `depth`, `amplitude`, `mode_x`, `mode_y` | 0.25, 0.01, 1, 1 | initial surface |
| `rho_liquid`, `rho_gas`, `mu_liquid`, `mu_gas`, `gravity`, `wall_bc` | water/air, `slip` | fluids and walls |
| `t_end`, `cfl`, `dt_max`, `output_dt` | 10, 0.3, 0.005, 0.05 | time stepping [s] |
| `poisson_tol`, `poisson_maxit` | 1e-9, 5000 | pressure solver |
| `output_dir`, `vtk_output`, `print_every`, `threads` | output/standing_wave_3d, 0, 100, 0 | output, OpenMP threads (0 = automatic) |

Exit codes are the same as in 2D: 0 done, 1 bad input, 2 diverged or output failure, 3 cancelled.

## Output (in `output_dir`)

The output is sized for limited disk space. The default run writes 201 frames, about 130 MB.

| File | Content |
|---|---|
| `manifest.json` | `"format": "freesurf3d-web-1", "dim": 3`. Contains the grid, `storage` description, `p_hydro` (flat-surface hydrostatic pressure per layer), display-field ranges for slices (C, p, p_dyn, speed, vorticity), `surface_range`, frame times and the summary |
| `frames/frame_NNNNNN.bin` | 5 **float16** cell-centred fields (C, p_dyn, u, v, w), field-major, x fastest, then y, then z |
| `surface/surface_NNNNNN.bin` | float32 surface height η(x, y), x fastest |
| `probes.csv` | η at the corner (0,0), the centre and the far corner; volume, KE, PE, max \|∇·u\|, Poisson iterations |
| `vtk/fields.pvd` (if `vtk_output=1`) | 3D ImageData (C, p, velocity). Use *Contour* at C = 0.5 in ParaView for the surface |

The web app turns stored frames into slices on request: pressure p = p_dyn + p_hydro, the speed, and the vorticity component normal to the slice (see `web/lib/slices.ts`).
