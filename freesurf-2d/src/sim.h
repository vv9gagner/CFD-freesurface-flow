#ifndef FS_SIM_H
#define FS_SIM_H

#include <stddef.h>

#include "config.h"

/* Number of ghost layers around the domain. */
#define NG 2

/*
 * Staggered (MAC) grid on the tank [0, Lx] x [0, Lz]:
 *   cell (i, j), i = 0..nx-1, j = 0..nz-1, centre at ((i+0.5)dx, (j+0.5)dz)
 *   u(i, j) lives on the x-face at (i dx, (j+0.5)dz),   i = 0..nx
 *   w(i, j) lives on the z-face at ((i+0.5)dx, j dz),   j = 0..nz
 * Every array has the same (nx+1+2NG) x (nz+1+2NG) layout so one index
 * function serves cells, faces and corners alike.
 */
typedef struct {
    Params prm;
    int nx, nz;
    int S; /* row stride */
    size_t ntot;
    double dx, dz;
    int noslip;

    double *u, *w;    /* face velocities */
    double *us, *ws;  /* predicted (pre-projection) velocities */
    double *p;        /* solver pressure p_rgh = p + rho g (z - depth) at cell centres */
    double *C;        /* liquid volume fraction */
    double *rho, *mu; /* cell-centred density and viscosity */

    /* VOF work arrays: start-of-step C, PLIC normal + line constant, face fluxes */
    double *C0, *pn_x, *pn_z, *palpha, *flux;
    /* shear stress at cell corners; gravity jump at the interface (reduced gravity) */
    double *tau, *phi;
    /* pressure Poisson matrix, MIC(0) preconditioner and PCG vectors */
    double *Adiag, *Ax, *Az, *precon, *rhs, *res, *zv, *sv, *qv, *tmp;
    int ref_i, ref_j; /* cell where p is pinned to 0 */

    double t, dt;
    long step;
    int pois_iters;
    double pois_res;
} Sim;

static inline size_t IX(const Sim *s, int i, int j)
{
    return (size_t)(j + NG) * (size_t)s->S + (size_t)(i + NG);
}

Sim *sim_create(const Params *prm);
void sim_destroy(Sim *s);

/* Cell-centred interpolants used by diagnostics and output. */
double sim_cell_u(const Sim *s, int i, int j);
double sim_cell_w(const Sim *s, int i, int j);
double sim_cell_vorticity(const Sim *s, int i, int j);
/* Physical pressure p = p_rgh - rho g (z - depth), relative to the pinned top-left cell. */
double sim_cell_pressure(const Sim *s, int i, int j);

#endif
