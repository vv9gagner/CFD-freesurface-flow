#include "bc.h"

void bc_u(const Sim *s, double *u)
{
    const int nx = s->nx, nz = s->nz;
    const double tang = s->noslip ? -1.0 : 1.0;

    /* walls x = 0 and x = Lx: u is the normal component, odd reflection */
    for (int j = 0; j < nz; j++) {
        u[IX(s, 0, j)] = 0.0;
        u[IX(s, nx, j)] = 0.0;
        for (int k = 1; k <= NG; k++) {
            u[IX(s, -k, j)] = -u[IX(s, k, j)];
            u[IX(s, nx + k, j)] = -u[IX(s, nx - k, j)];
        }
    }
    /* walls z = 0 and z = Lz: u is tangential */
    for (int i = -NG; i <= nx + NG; i++) {
        for (int k = 1; k <= NG; k++) {
            u[IX(s, i, -k)] = tang * u[IX(s, i, k - 1)];
            u[IX(s, i, nz - 1 + k)] = tang * u[IX(s, i, nz - k)];
        }
    }
}

void bc_w(const Sim *s, double *w)
{
    const int nx = s->nx, nz = s->nz;
    const double tang = s->noslip ? -1.0 : 1.0;

    for (int i = 0; i < nx; i++) {
        w[IX(s, i, 0)] = 0.0;
        w[IX(s, i, nz)] = 0.0;
        for (int k = 1; k <= NG; k++) {
            w[IX(s, i, -k)] = -w[IX(s, i, k)];
            w[IX(s, i, nz + k)] = -w[IX(s, i, nz - k)];
        }
    }
    for (int j = -NG; j <= nz + NG; j++) {
        for (int k = 1; k <= NG; k++) {
            w[IX(s, -k, j)] = tang * w[IX(s, k - 1, j)];
            w[IX(s, nx - 1 + k, j)] = tang * w[IX(s, nx - k, j)];
        }
    }
}

void bc_scalar(const Sim *s, double *c)
{
    const int nx = s->nx, nz = s->nz;
    for (int j = 0; j < nz; j++) {
        for (int k = 1; k <= NG; k++) {
            c[IX(s, -k, j)] = c[IX(s, k - 1, j)];
            c[IX(s, nx - 1 + k, j)] = c[IX(s, nx - k, j)];
        }
    }
    for (int i = -NG; i < nx + NG; i++) {
        for (int k = 1; k <= NG; k++) {
            c[IX(s, i, -k)] = c[IX(s, i, k - 1)];
            c[IX(s, i, nz - 1 + k)] = c[IX(s, i, nz - k)];
        }
    }
}
