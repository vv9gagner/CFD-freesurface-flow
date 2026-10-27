#ifndef FS3_BC_H
#define FS3_BC_H

#include "sim.h"

/* Closed tank: zero normal velocity on all six walls, free-slip or no-slip tangential
 * velocity, zero-gradient scalars (90 degree contact angle). */
void bc_u(const Sim *s, double *u);
void bc_v(const Sim *s, double *v);
void bc_w(const Sim *s, double *w);
void bc_scalar(const Sim *s, double *c);

#endif
