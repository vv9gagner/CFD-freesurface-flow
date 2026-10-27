#ifndef FS_INIT_H
#define FS_INIT_H

#include "sim.h"

/* Set C, rho, mu, u, w and an initial pressure guess for the configured case. Returns 0 on success. */
int init_case(Sim *s);

/* Liquid fraction of a cell for a flat surface at height h (row j). */
double flat_fraction(const Sim *s, int j, double h);

#endif
