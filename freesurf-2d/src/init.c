#include "init.h"
#include "bc.h"
#include "util.h"
#include "vof.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

double flat_fraction(const Sim *s, int j, double h)
{
    return clamp01((h - j * s->dz) / s->dz);
}

/* eta(x) = h + A cos(k x): cell fractions by midpoint quadrature in x of the exact column fill. */
static void init_standing_wave(Sim *s)
{
    const Params *p = &s->prm;
    const double k = p->mode * FS_PI / p->Lx;
    const int nq = 64;
    for (int j = 0; j < s->nz; j++)
        for (int i = 0; i < s->nx; i++) {
            double sum = 0.0;
            for (int q = 0; q < nq; q++) {
                double x = (i + (q + 0.5) / nq) * s->dx;
                double eta = p->depth + p->amplitude * cos(k * x);
                sum += clamp01((eta - j * s->dz) / s->dz);
            }
            s->C[IX(s, i, j)] = sum / nq;
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

    /* fluid at rest; p_rgh is nearly uniform for a surface close to z = depth, so 0 is a good first guess */
    for (size_t k = 0; k < s->ntot; k++)
        s->u[k] = s->w[k] = s->p[k] = 0.0;
    s->t = 0.0;
    s->step = 0;
    return 0;
}
