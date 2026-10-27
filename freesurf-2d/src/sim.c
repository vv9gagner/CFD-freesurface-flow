#include "sim.h"

#include <stdlib.h>
#include <string.h>

Sim *sim_create(const Params *prm)
{
    Sim *s = calloc(1, sizeof *s);
    if (!s)
        return NULL;
    s->prm = *prm;
    s->nx = prm->nx;
    s->nz = prm->nz;
    s->S = prm->nx + 1 + 2 * NG;
    s->ntot = (size_t)s->S * (size_t)(prm->nz + 1 + 2 * NG);
    s->dx = prm->Lx / prm->nx;
    s->dz = prm->Lz / prm->nz;
    s->noslip = strcmp(prm->wall_bc, "noslip") == 0;
    s->ref_i = 0;
    s->ref_j = prm->nz - 1; /* top-left cell, in the gas */

    double **arrays[] = {
        &s->u, &s->w, &s->us, &s->ws, &s->p, &s->C, &s->rho, &s->mu,
        &s->C0, &s->pn_x, &s->pn_z, &s->palpha, &s->flux, &s->tau, &s->phi,
        &s->Adiag, &s->Ax, &s->Az, &s->precon, &s->rhs, &s->res, &s->zv, &s->sv, &s->qv, &s->tmp,
    };
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
    double *arrays[] = {
        s->u, s->w, s->us, s->ws, s->p, s->C, s->rho, s->mu,
        s->C0, s->pn_x, s->pn_z, s->palpha, s->flux, s->tau, s->phi,
        s->Adiag, s->Ax, s->Az, s->precon, s->rhs, s->res, s->zv, s->sv, s->qv, s->tmp,
    };
    for (size_t n = 0; n < sizeof arrays / sizeof arrays[0]; n++)
        free(arrays[n]);
    free(s);
}

double sim_cell_u(const Sim *s, int i, int j)
{
    return 0.5 * (s->u[IX(s, i, j)] + s->u[IX(s, i + 1, j)]);
}

double sim_cell_w(const Sim *s, int i, int j)
{
    return 0.5 * (s->w[IX(s, i, j)] + s->w[IX(s, i, j + 1)]);
}

/* Vorticity dw/dx - du/dz at corner (i, j); needs ghost-filled velocities. */
static double corner_vorticity(const Sim *s, int i, int j)
{
    return (s->w[IX(s, i, j)] - s->w[IX(s, i - 1, j)]) / s->dx
         - (s->u[IX(s, i, j)] - s->u[IX(s, i, j - 1)]) / s->dz;
}

double sim_cell_vorticity(const Sim *s, int i, int j)
{
    return 0.25 * (corner_vorticity(s, i, j) + corner_vorticity(s, i + 1, j)
                   + corner_vorticity(s, i, j + 1) + corner_vorticity(s, i + 1, j + 1));
}

static double raw_pressure(const Sim *s, int i, int j)
{
    size_t k = IX(s, i, j);
    return s->p[k] - s->rho[k] * s->prm.gravity * ((j + 0.5) * s->dz - s->prm.depth);
}

double sim_cell_pressure(const Sim *s, int i, int j)
{
    return raw_pressure(s, i, j) - raw_pressure(s, s->ref_i, s->ref_j);
}
