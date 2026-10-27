#ifndef FS3_SIM_H
#define FS3_SIM_H

#include <stddef.h>

#include "config.h"

/* Number of ghost layers around the domain. */
#define NG 2

/*
 * Staggered (MAC) grid on the tank [0,Lx] x [0,Ly] x [0,Lz], z vertical:
 *   cell (i,j,k) centre ((i+.5)dx, (j+.5)dy, (k+.5)dz)
 *   u(i,j,k) on the x-face at x = i dx,  i = 0..nx
 *   v(i,j,k) on the y-face at y = j dy,  j = 0..ny
 *   w(i,j,k) on the z-face at z = k dz,  k = 0..nz
 * Every array uses the same (nx+1+2NG)(ny+1+2NG)(nz+1+2NG) layout, x fastest.
 */
typedef struct {
    Params prm;
    int nx, ny, nz;
    ptrdiff_t Sx, Sy; /* strides of j and k */
    size_t ntot;
    double dx, dy, dz;
    int noslip;

    double *u, *v, *w;    /* face velocities */
    double *us, *vs, *ws; /* predicted velocities */
    double *p;            /* solver pressure p_rgh = p + rho g (z - depth) */
    double *C;            /* liquid volume fraction */
    double *rho, *mu;

    /* VOF work: start-of-step C, PLIC normal + plane constant, face fluxes */
    double *C0, *pn_x, *pn_y, *pn_z, *palpha, *flux;
    /* edge shear stresses and the reduced-gravity jump at the interface */
    double *tau_xy, *tau_xz, *tau_yz, *phi;
    /* pressure matrix, MIC(0) preconditioner, PCG vectors */
    double *Adiag, *Ax, *Ay, *Az, *precon, *rhs, *res, *zv, *sv, *qv, *tmp;
    int ref_i, ref_j, ref_k; /* cell where p is pinned to 0 */
    int pc_axis, pc_blocks;  /* block-Jacobi MIC(0): slabs along this axis, one per thread */

    double t, dt;
    long step;
    int pois_iters;
    double pois_res;
} Sim;

static inline size_t IX(const Sim *s, int i, int j, int k)
{
    return (size_t)((k + NG) * s->Sy + (j + NG) * s->Sx + (i + NG));
}

Sim *sim_create(const Params *prm);
void sim_destroy(Sim *s);

/* Cell-centred interpolants (velocities need ghost-filled arrays). */
double sim_cell_u(const Sim *s, int i, int j, int k);
double sim_cell_v(const Sim *s, int i, int j, int k);
double sim_cell_w(const Sim *s, int i, int j, int k);
/* Physical pressure p = p_rgh - rho g (z - depth), relative to the pinned cell. */
double sim_cell_pressure(const Sim *s, int i, int j, int k);

#endif
