#ifndef FS3_INIT_H
#define FS3_INIT_H

#include "sim.h"

/* Set C, rho, mu, u, v, w and the pressure guess for the configured case. Returns 0 on success. */
int init_case(Sim *s);

/* Liquid fraction of a cell in layer k for a flat surface at height h. */
double flat_fraction(const Sim *s, int k, double h);

#endif
