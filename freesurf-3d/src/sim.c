#include "sim.h"

#include <stdlib.h>
#include <string.h>

#define SIM_ARRAYS(s)                                                                                     \
    &(s)->u, &(s)->v, &(s)->w, &(s)->us, &(s)->vs, &(s)->ws, &(s)->p, &(s)->C, &(s)->rho, &(s)->mu,      \
        &(s)->C0, &(s)->pn_x, &(s)->pn_y, &(s)->pn_z, &(s)->palpha, &(s)->flux, &(s)->tau_xy,           \
        &(s)->tau_xz, &(s)->tau_yz, &(s)->phi, &(s)->Adiag, &(s)->Ax, &(s)->Ay, &(s)->Az, &(s)->precon, \
        &(s)->rhs, &(s)->res, &(s)->zv, &(s)->sv, &(s)->qv, &(s)->tmp

Sim *sim_create(const Params *prm)
{
    Sim *s = calloc(1, sizeof *s);
    if (!s)
        return NULL;
    s->prm = *prm;
    s->nx = prm->nx;
    s->ny = prm->ny;
    s->nz = prm->nz;
    s->Sx = prm->nx + 1 + 2 * NG;
    s->Sy = s->Sx * (prm->ny + 1 + 2 * NG);
    s->ntot = (size_t)s->Sy * (size_t)(prm->nz + 1 + 2 * NG);
    s->dx = prm->Lx / prm->nx;
    s->dy = prm->Ly / prm->ny;
    s->dz = prm->Lz / prm->nz;
    s->noslip = strcmp(prm->wall_bc, "noslip") == 0;
    s->ref_i = 0;
    s->ref_j = 0;
    s->ref_k = prm->nz - 1; /* top corner cell, in the gas */

    double **arrays[] = { SIM_ARRAYS(s) };
    for (size_t n = 0; n < sizeof arrays / sizeof arrays[0]; n++) {
        *arrays[n] = calloc(s->ntot, sizeof(double));
        if (!*arrays[n]) {
            sim_destroy(s);
            return NULL;
        }
    }
    return s;
}

void sim_destroy(Sim *s)
{
    if (!s)
        return;
    double **arrays[] = { SIM_ARRAYS(s) };
    for (size_t n = 0; n < sizeof arrays / sizeof arrays[0]; n++)
        free(*arrays[n]);
    free(s);
}

double sim_cell_u(const Sim *s, int i, int j, int k)
{
    return 0.5 * (s->u[IX(s, i, j, k)] + s->u[IX(s, i + 1, j, k)]);
}

double sim_cell_v(const Sim *s, int i, int j, int k)
{
    return 0.5 * (s->v[IX(s, i, j, k)] + s->v[IX(s, i, j + 1, k)]);
}

double sim_cell_w(const Sim *s, int i, int j, int k)
{
    return 0.5 * (s->w[IX(s, i, j, k)] + s->w[IX(s, i, j, k + 1)]);
}

static double raw_pressure(const Sim *s, int i, int j, int k)
{
    size_t n = IX(s, i, j, k);
    return s->p[n] - s->rho[n] * s->prm.gravity * ((k + 0.5) * s->dz - s->prm.depth);
}

double sim_cell_pressure(const Sim *s, int i, int j, int k)
{
    return raw_pressure(s, i, j, k) - raw_pressure(s, s->ref_i, s->ref_j, s->ref_k);
}
