#ifndef FS_VOF_H
#define FS_VOF_H

#include "sim.h"

/*
 * PLIC geometry in a unit cell [0,1]^2. The interface is the line
 * nx*x + nz*z = alpha and the liquid occupies nx*x + nz*z <= alpha,
 * i.e. the normal points out of the liquid.
 */

/* Liquid fraction of the unit square for the line (nx, nz, alpha). Any normal. */
double plic_area(double nx, double nz, double alpha);

/* Line constant giving liquid fraction c; the normal must satisfy |nx|+|nz| = 1. */
double plic_alpha(double c, double nx, double nz);

/* Mixed Youngs-Centred normal of cell (i, j), normalised to |nx|+|nz| = 1. */
void plic_normal(const Sim *s, const double *c, int i, int j, double *nx, double *nz);

/* End points (x0, z0, x1, z1) of the interface segment in the unit cell. Returns 0 if degenerate. */
int plic_segment(double nx, double nz, double alpha, double pts[4]);

/* Compute normal and line constant for every interface cell of c (ghosts must be filled). */
void vof_reconstruct(Sim *s, const double *c);

/* Advance C by dt with the face velocities u, w (Weymouth & Yue 2010 split scheme). */
void vof_advect(Sim *s, double dt);

/* rho, mu from C (including ghost cells). */
void vof_update_properties(Sim *s);

#endif
