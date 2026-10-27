/* 3D VOF advection: a sphere makes one revolution about the y axis, then one about the
 * x axis, so every sweep direction is exercised. The discrete velocity fields are exactly
 * divergence-free, so liquid volume must be conserved to round-off and C must stay in [0, 1]. */
#include "bc.h"
#include "config.h"
#include "sim.h"
#include "util.h"
#include "vof.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static double total(const Sim *s)
{
    double v = 0.0;
    for (int k = 0; k < s->nz; k++)
        for (int j = 0; j < s->ny; j++)
            for (int i = 0; i < s->nx; i++)
                v += s->C[IX(s, i, j, k)];
    return v;
}

/* axis 1: rotation about y (u, w); axis 0: rotation about x (v, w) */
static void set_rotation(Sim *s, int axis, double omega)
{
    const int n = s->nx;
    for (int k = -NG; k <= n + NG; k++)
        for (int j = -NG; j <= n + NG; j++)
            for (int i = -NG; i <= n + NG; i++) {
                size_t m = IX(s, i, j, k);
                double zc = (k + 0.5) * s->dz - 0.5;
                if (axis == 1) {
                    s->u[m] = -omega * zc;
                    s->v[m] = 0.0;
                    s->w[m] = omega * ((i + 0.5) * s->dx - 0.5);
                } else {
                    s->u[m] = 0.0;
                    s->v[m] = -omega * zc;
                    s->w[m] = omega * ((j + 0.5) * s->dy - 0.5);
                }
            }
}

int main(void)
{
    const int n = 40;
    const double r = 0.15, omega = 2.0 * FS_PI;
    Params p;
    params_defaults(&p);
    p.nx = p.ny = p.nz = n;
    p.Lx = p.Ly = p.Lz = 1.0;
    Sim *s = sim_create(&p);
    if (!s)
        return 1;

    const int q = 8;
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                int in = 0;
                for (int c = 0; c < q; c++)
                    for (int b = 0; b < q; b++)
                        for (int a = 0; a < q; a++) {
                            double x = (i + (a + 0.5) / q) * s->dx - 0.5, y = (j + (b + 0.5) / q) * s->dy - 0.5;
                            double z = (k + (c + 0.5) / q) * s->dz - 0.75;
                            in += x * x + y * y + z * z < r * r;
                        }
                s->C[IX(s, i, j, k)] = in / (double)(q * q * q);
            }
    bc_scalar(s, s->C);
    double *C_init = malloc(s->ntot * sizeof(double));
    if (!C_init)
        return 1;
    for (size_t m = 0; m < s->ntot; m++)
        C_init[m] = s->C[m];

    const double umax = omega * 0.5 * sqrt(2.0);
    const int steps = (int)ceil(1.0 / (0.4 * s->dx / umax));
    const double dt = 1.0 / steps, v0 = total(s);
    double cmin = 1.0, cmax = 0.0, worst = 0.0;
    s->step = 0;
    for (int axis = 1; axis >= 0; axis--) {
        set_rotation(s, axis, omega);
        for (int m = 0; m < steps; m++, s->step++) {
            vof_advect(s, dt);
            worst = fmax(worst, fabs(total(s) - v0) / v0);
        }
    }
    double l1 = 0.0;
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                double c = s->C[IX(s, i, j, k)];
                cmin = fmin(cmin, c);
                cmax = fmax(cmax, c);
                l1 += fabs(c - C_init[IX(s, i, j, k)]);
            }
    l1 *= s->dx * s->dy * s->dz / (4.0 / 3.0 * FS_PI * r * r * r);

    printf("  2 x %d steps, max volume drift %.3e, C in [%.3g, %.3g], shape error L1/volume %.4f\n", steps, worst,
           cmin, cmax, l1);
    int ok = worst < 1e-12 && cmin >= 0.0 && cmax <= 1.0 && l1 < 0.08;
    printf(ok ? "VOF 3D: ok\n" : "VOF 3D: FAILED\n");
    free(C_init);
    sim_destroy(s);
    return ok ? 0 : 1;
}
