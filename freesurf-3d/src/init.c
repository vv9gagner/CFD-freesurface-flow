#include "init.h"
#include "bc.h"
#include "util.h"
#include "vof.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

double flat_fraction(const Sim *s, int k, double h)
{
    return clamp01((h - k * s->dz) / s->dz);
}

/* eta(x,y) = h + A cos(m pi x / Lx) cos(n pi y / Ly): cell fractions by midpoint quadrature
 * over the cell's horizontal cross-section of the exact column fill. */
static void init_standing_wave(Sim *s)
{
    const Params *p = &s->prm;
    const double kx = p->mode_x * FS_PI / p->Lx, ky = p->mode_y * FS_PI / p->Ly;
    const int nq = 16;
#pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < s->ny; j++)
        for (int i = 0; i < s->nx; i++)
            for (int k = 0; k < s->nz; k++) {
                double sum = 0.0;
                for (int qy = 0; qy < nq; qy++)
                    for (int qx = 0; qx < nq; qx++) {
                        double x = (i + (qx + 0.5) / nq) * s->dx, y = (j + (qy + 0.5) / nq) * s->dy;
                        double eta = p->depth + p->amplitude * cos(kx * x) * cos(ky * y);
                        sum += clamp01((eta - k * s->dz) / s->dz);
                    }
                s->C[IX(s, i, j, k)] = sum / (nq * nq);
            }
}

int init_case(Sim *s)
{
    if (strcmp(s->prm.case_name, "standing_wave") == 0)
        init_standing_wave(s);
    else {
        fprintf(stderr, "unknown case '%s'\n", s->prm.case_name);
        return -1;
    }
    bc_scalar(s, s->C);
    vof_update_properties(s);
    /* fluid at rest; p_rgh is nearly uniform for a surface close to z = depth */
    for (size_t n = 0; n < s->ntot; n++)
        s->u[n] = s->v[n] = s->w[n] = s->p[n] = 0.0;
    s->t = 0.0;
    s->step = 0;
    return 0;
}
