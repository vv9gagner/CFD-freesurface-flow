/* 3D pressure solver: second-order convergence (constant density) and an exact
 * algebraic solve across a 1000:1 density jump. */
#include "config.h"
#include "poisson.h"
#include "sim.h"
#include "util.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static Sim *make_sim(int n, double rho_low, double rho_high)
{
    Params p;
    params_defaults(&p);
    p.nx = p.ny = p.nz = n;
    p.Lx = p.Ly = p.Lz = 1.0;
    p.poisson_tol = 1e-12;
    p.poisson_maxit = 20000;
    Sim *s = sim_create(&p);
    if (!s)
        return NULL;
    for (int k = -NG; k < n + NG; k++)
        for (int j = -NG; j < n + NG; j++)
            for (int i = -NG; i < n + NG; i++)
                s->rho[IX(s, i, j, k)] = (k + 0.5) * s->dz < 0.5 ? rho_high : rho_low;
    poisson_setup(s);
    return s;
}

static int is_ref(const Sim *s, int i, int j, int k)
{
    return i == s->ref_i && j == s->ref_j && k == s->ref_k;
}

static double exact(const Sim *s, int i, int j, int k)
{
    return cos(FS_PI * (i + 0.5) * s->dx) * cos(FS_PI * (j + 0.5) * s->dy) * cos(FS_PI * (k + 0.5) * s->dz);
}

/* -lap p = 3 pi^2 p for p = cos(pi x) cos(pi y) cos(pi z), which satisfies the Neumann walls */
static double convergence_error(int n, int *iters)
{
    Sim *s = make_sim(n, 1.0, 1.0);
    if (!s)
        return NAN;
    double mean_x = 0.0, mean_e = 0.0;
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                size_t m = IX(s, i, j, k);
                s->rhs[m] = is_ref(s, i, j, k) ? 0.0 : 3.0 * FS_PI * FS_PI * exact(s, i, j, k);
                s->p[m] = 0.0;
            }
    *iters = poisson_solve(s, s->rhs, s->p);
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                mean_x += s->p[IX(s, i, j, k)];
                mean_e += exact(s, i, j, k);
            }
    double cells = (double)n * n * n;
    mean_x /= cells;
    mean_e /= cells;
    double err = 0.0;
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
                err = fmax(err, fabs((s->p[IX(s, i, j, k)] - mean_x) - (exact(s, i, j, k) - mean_e)));
    sim_destroy(s);
    return err;
}

int main(void)
{
    int ok = 1;
    printf("  constant density, p = cos(pi x) cos(pi y) cos(pi z):\n");
    double prev = NAN;
    for (int n = 16; n <= 64; n *= 2) {
        int it = -1;
        double e = convergence_error(n, &it);
        double order = isnan(prev) ? NAN : log2(prev / e);
        printf("    n=%3d  iters=%4d  max error=%.3e  order=%.2f\n", n, it, e, order);
        if (it < 0 || (!isnan(order) && order < 1.8))
            ok = 0;
        prev = e;
    }

    printf("  density jump 1000:1 at z = 0.5, smooth exact solution:\n");
    const int n = 40;
    Sim *s = make_sim(n, 1.0, 1000.0);
    if (!s)
        return 1;
    double *truth = calloc(s->ntot, sizeof(double));
    if (!truth)
        return 1;
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                double x = (i + 0.5) * s->dx, y = (j + 0.5) * s->dy, z = (k + 0.5) * s->dz;
                truth[IX(s, i, j, k)] = is_ref(s, i, j, k) ? 0.0 : sin(3.0 * x) * cos(2.0 * y) * cos(z) + z;
            }
    poisson_apply(s, truth, s->rhs);
    for (size_t m = 0; m < s->ntot; m++)
        s->p[m] = 0.0;
    int it = poisson_solve(s, s->rhs, s->p);
    double err = 0.0, scale = 0.0;
    for (int k = 0; k < n; k++)
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++) {
                size_t m = IX(s, i, j, k);
                err = fmax(err, fabs(s->p[m] - truth[m]));
                scale = fmax(scale, fabs(truth[m]));
            }
    printf("    n=%3d  iters=%4d  relative error=%.3e\n", n, it, err / scale);
    if (it < 0 || err / scale > 1e-6)
        ok = 0;
    free(truth);
    sim_destroy(s);

    printf(ok ? "Poisson 3D: ok\n" : "Poisson 3D: FAILED\n");
    return ok ? 0 : 1;
}
