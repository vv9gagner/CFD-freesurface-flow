#include "solver.h"
#include "momentum.h"
#include "poisson.h"
#include "util.h"
#include "vof.h"

#include <math.h>

double solver_dt(const Sim *s)
{
    const Params *p = &s->prm;
    double umax = 0.0, wmax = 0.0;
    for (int j = 0; j < s->nz; j++)
        for (int i = 0; i <= s->nx; i++)
            umax = dmax(umax, fabs(s->u[IX(s, i, j)]));
    for (int j = 0; j <= s->nz; j++)
        for (int i = 0; i < s->nx; i++)
            wmax = dmax(wmax, fabs(s->w[IX(s, i, j)]));

    const double h = dmin(s->dx, s->dz);
    double dt = p->dt_max;
    if (umax > 0.0)
        dt = dmin(dt, p->cfl * s->dx / umax);
    if (wmax > 0.0)
        dt = dmin(dt, p->cfl * s->dz / wmax);

    /* explicit viscosity, worst case: liquid viscosity on a gas-density face */
    double nu_max = dmax(p->mu_l, p->mu_g) / p->rho_g;
    if (nu_max > 0.0)
        dt = dmin(dt, 0.125 * h * h / nu_max);

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
