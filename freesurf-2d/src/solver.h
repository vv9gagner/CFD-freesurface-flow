#ifndef FS_SOLVER_H
#define FS_SOLVER_H

#include "sim.h"

/* Stable time step from the advective CFL, explicit viscosity, gravity waves and dt_max. */
double solver_dt(const Sim *s);

/*
 * One time step:
 *   1. advect C with u^n (geometric VOF)       -> C^{n+1}, rho^{n+1}, mu^{n+1}
 *   2. momentum predictor with rho^{n+1}       -> u*
 *   3. pressure projection                     -> u^{n+1}, p^{n+1}
 * Position-then-velocity ordering is a symplectic-Euler splitting, so wave energy does not drift.
 * Returns the Poisson iteration count, or -1 if the pressure solve did not converge.
 */
int solver_step(Sim *s, double dt);

#endif
