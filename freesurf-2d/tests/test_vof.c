/* VOF advection: one revolution of a disk in solid-body rotation.
 * The discrete velocity field is exactly divergence-free, so liquid volume must be
 * conserved to round-off and C must stay within [0, 1]. */
#include "bc.h"
#include "config.h"
#include "sim.h"
#include "util.h"
#include "vof.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(void)
{
    const int n = 100;
    const double r = 0.15, xc = 0.5, zc = 0.75, omega = 2.0 * FS_PI;
    Params p;
    params_defaults(&p);
    p.nx = p.nz = n;
    p.Lx = p.Lz = 1.0;
    Sim *s = sim_create(&p);
    if (!s)
        return 1;

    /* disk by 16x16 sub-sampling */
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            int in = 0;
            for (int a = 0; a < 16; a++)
                for (int b = 0; b < 16; b++) {
                    double x = (i + (a + 0.5) / 16) * s->dx, z = (j + (b + 0.5) / 16) * s->dz;
                    in += sq(x - xc) + sq(z - zc) < r * r;
                }
            s->C[IX(s, i, j)] = in / 256.0;
        }
    bc_scalar(s, s->C);
    double *C_init = malloc(s->ntot * sizeof(double));
    if (!C_init)
        return 1;
    for (size_t k = 0; k < s->ntot; k++)
        C_init[k] = s->C[k];

    /* u = -omega (z - 0.5) on x-faces, w = omega (x - 0.5) on z-faces */
    for (int j = -NG; j <= n + NG - 1; j++)
        for (int i = -NG; i <= n + NG; i++)
            s->u[IX(s, i, j)] = -omega * ((j + 0.5) * s->dz - 0.5);
    for (int j = -NG; j <= n + NG; j++)
        for (int i = -NG; i <= n + NG - 1; i++)
            s->w[IX(s, i, j)] = omega * ((i + 0.5) * s->dx - 0.5);

    double umax = omega * 0.5 * sqrt(2.0);
    int steps = (int)ceil(1.0 / (0.4 * s->dx / umax));
    double dt = 1.0 / steps;

    double v0 = 0.0;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            v0 += s->C[IX(s, i, j)];

    double cmin = 1.0, cmax = 0.0, worst_drift = 0.0;
    for (s->step = 0; s->step < steps; s->step++) {
        vof_advect(s, dt);
        double v = 0.0;
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                double c = s->C[IX(s, i, j)];
                v += c;
                cmin = fmin(cmin, c);
                cmax = fmax(cmax, c);
            }
        worst_drift = fmax(worst_drift, fabs(v - v0) / v0);
    }

    double l1 = 0.0;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++)
            l1 += fabs(s->C[IX(s, i, j)] - C_init[IX(s, i, j)]);
    l1 *= s->dx * s->dz / (FS_PI * r * r);

    printf("  %d steps, max volume drift %.3e, C in [%.3g, %.3g], shape error L1/area %.4f\n",
           steps, worst_drift, cmin, cmax, l1);
    int ok = worst_drift < 1e-12 && cmin >= 0.0 && cmax <= 1.0 && l1 < 0.05;
    printf(ok ? "VOF: ok\n" : "VOF: FAILED\n");
    free(C_init);
    sim_destroy(s);
    return ok ? 0 : 1;
}
