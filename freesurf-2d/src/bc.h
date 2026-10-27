#ifndef FS_BC_H
#define FS_BC_H

#include "sim.h"

/* Closed tank: zero normal velocity on all walls, free-slip or no-slip
 * tangential velocity, zero-gradient scalars (90 degree contact angle). */
void bc_u(const Sim *s, double *u);
void bc_w(const Sim *s, double *w);
void bc_scalar(const Sim *s, double *c);

#endif
