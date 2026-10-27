#ifndef FS3_VOF_H
#define FS3_VOF_H

#include "sim.h"

/*
 * PLIC geometry in a unit cube [0,1]^3. The interface is the plane
 * nx*x + ny*y + nz*z = alpha and the liquid occupies nx*x + ny*y + nz*z <= alpha,
 * i.e. the normal points out of the liquid.
 */

/* Liquid fraction of the unit square for the line n1*x + n2*y = alpha (any normal). */
double line_area(double n1, double n2, double alpha);

/* Liquid fraction of the unit cube for the plane (any normal). */
double plane_volume(double nx, double ny, double nz, double alpha);

/* Plane constant giving liquid fraction c (normal not all zero). */
double plane_alpha(double c, double nx, double ny, double nz);

/* Mixed Youngs-Centred normal of cell (i,j,k), normalised to |nx|+|ny|+|nz| = 1. */
void plic_normal(const Sim *s, const double *c, int i, int j, int k, double n[3]);

/* Height (0..1, unit-cube coordinates) of the plane above the centre of the cell's
 * horizontal cross-section, clamped to the cell; 0.5 for a near-vertical plane. */
double plane_height(double nx, double ny, double nz, double alpha);

/* Normal and plane constant for every interface cell of c (ghosts must be filled). */
void vof_reconstruct(Sim *s, const double *c);

/* Advance C by dt with the face velocities (Weymouth & Yue 2010, three split sweeps). */
void vof_advect(Sim *s, double dt);

/* rho, mu from C (including ghost cells). */
void vof_update_properties(Sim *s);

#endif
