#include "solver.h"
#include "momentum.h"
#include "poisson.h"
#include "util.h"
#include "vof.h"

#include <math.h>

static double max_abs(const Sim *s, const double *a, int ni, int nj, int nk)
{
    double m = 0.0;
#pragma omp parallel for collapse(2) reduction(max : m) schedule(static)
    for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
            for (int i = 0; i < ni; i++)
                m = dmax(m, fabs(a[IX(s, i, j, k)]));
    return m;
}

double solver_dt(const Sim *s)
{
    const Params *p = &s->prm;
    const double umax = max_abs(s, s->u, s->nx + 1, s->ny, s->nz);
    const double vmax = max_abs(s, s->v, s->nx, s->ny + 1, s->nz);
    const double wmax = max_abs(s, s->w, s->nx, s->ny, s->nz + 1);

    const double h = dmin(s->dx, dmin(s->dy, s->dz));
    double dt = p->dt_max;
    if (umax > 0.0)
        dt = dmin(dt, p->cfl * s->dx / umax);
    if (vmax > 0.0)
        dt = dmin(dt, p->cfl * s->dy / vmax);
    if (wmax > 0.0)
        dt = dmin(dt, p->cfl * s->dz / wmax);

    /* explicit viscosity in 3D (normal stresses carry 2 mu), worst case: liquid viscosity on a gas face */
    double nu_max = dmax(p->mu_l, p->mu_g) / p->rho_g;
    if (nu_max > 0.0)
        dt = dmin(dt, 0.1 * h * h / nu_max);

    /* shortest resolved gravity wave (k = pi/h) needs omega dt < 2; keep a wide margin */
    if (p->gravity > 0.0)
        dt = dmin(dt, 0.5 * sqrt(h / p->gravity));
    return dt;
}

int solver_step(Sim *s, double dt)
{
    s->dt = dt;
    vof_advect(s, dt);
    vof_update_properties(s);
    momentum_predict(s, dt);
    int iters = projection(s, dt);
    s->t += dt;
    s->step++;
    return iters;
}
