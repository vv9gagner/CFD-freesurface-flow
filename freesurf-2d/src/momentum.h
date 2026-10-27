#ifndef FS_MOMENTUM_H
#define FS_MOMENTUM_H

#include "sim.h"

/* Explicit predictor: us, ws = u, w + dt * (-u.grad u + div(2 mu D)/rho + g (z - h) grad(rho)/rho).
 * Gravity uses the balanced-force "p_rgh" form: the solver's pressure is p_rgh = p + rho g (z - h),
 * so -grad p + rho g = -grad p_rgh + g (z - h) grad rho, and both terms are discretised on the
 * same faces. This removes spurious air currents at a tilted interface.
 * Uses the current rho and mu; wall-normal faces stay zero. */
void momentum_predict(Sim *s, double dt);

#endif
