#ifndef FS3_POISSON_H
#define FS3_POISSON_H

#include "sim.h"

/*
 * Variable-coefficient pressure equation  -div( (1/rho) grad p ) = b  on the cell centres
 * (7-point stencil), homogeneous Neumann on all walls, p pinned to 0 in cell (ref_i, ref_j, ref_k).
 */

/* Assemble the matrix from s->rho and build the MIC(0) preconditioner. */
void poisson_setup(Sim *s);

/* Preconditioned CG on A x = b with x as the initial guess; stops when
 * max|r| <= poisson_tol * max|b|. Returns the iteration count, or -1 if not converged. */
int poisson_solve(Sim *s, const double *b, double *x);

/* y = A x on the unknowns (exposed for tests). */
void poisson_apply(const Sim *s, const double *x, double *y);

/* Projection: solve for p from (us, vs, ws) and make (u, v, w) divergence-free. */
int projection(Sim *s, double dt);

#endif
