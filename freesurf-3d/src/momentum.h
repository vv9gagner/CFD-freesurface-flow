#ifndef FS3_MOMENTUM_H
#define FS3_MOMENTUM_H

#include "sim.h"

/*
 * Explicit predictor: (us,vs,ws) = (u,v,w) + dt * (-u.grad u + div(2 mu D)/rho + gravity).
 * Gravity uses the balanced-force "reduced gravity" form: the solver's pressure is
 * p_rgh = p + rho g (z - h), and gravity acts as the interfacial force
 * [rho] g (z_int - h) grad(C) / rho, with z_int the PLIC plane height in each cut cell.
 * Uses the current rho and mu; wall-normal faces stay zero.
 */
void momentum_predict(Sim *s, double dt);

#endif
