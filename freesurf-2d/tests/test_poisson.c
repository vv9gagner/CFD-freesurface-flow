/* Pressure solver: second-order convergence (constant density) and an exact
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
    p.nx = p.nz = n;
    p.Lx = p.Lz = 1.0;
    p.poisson_tol = 1e-12;
    p.poisson_maxit = 20000;
    Sim *s = sim_create(&p);
    if (!s)
        return NULL;
    for (int j = -NG; j < n + NG; j++)
        for (int i = -NG; i < n + NG; i++)
            s->rho[IX(s, i, j)] = (j + 0.5) * s->dz < 0.5 ? rho_high : rho_low;
    poisson_setup(s);
    return s;
}

static int is_ref(const Sim *s, int i, int j) { return i == s->ref_i && j == s->ref_j; }

/* -lap p = 2 pi^2 p for p = cos(pi x) cos(pi z), which satisfies the Neumann walls */
static double convergence_error(int n, int *iters)
{
    Sim *s = make_sim(n, 1.0, 1.0);
    if (!s)
        return NAN;
    double mean_x = 0.0, mean_e = 0.0;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            size_t k = IX(s, i, j);
            double pe = cos(FS_PI * (i + 0.5) * s->dx) * cos(FS_PI * (j + 0.5) * s->dz);
            s->rhs[k] = is_ref(s, i, j) ? 0.0 : 2.0 * FS_PI * FS_PI * pe;
            s->p[k] = 0.0;
        }
    *iters = poisson_solve(s, s->rhs, s->p);
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            mean_x += s->p[IX(s, i, j)];
            mean_e += cos(FS_PI * (i + 0.5) * s->dx) * cos(FS_PI * (j + 0.5) * s->dz);
        }
    mean_x /= (double)n * n;
    mean_e /= (double)n * n;
    double err = 0.0;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            double pe = cos(FS_PI * (i + 0.5) * s->dx) * cos(FS_PI * (j + 0.5) * s->dz);
            err = fmax(err, fabs((s->p[IX(s, i, j)] - mean_x) - (pe - mean_e)));
        }
    sim_destroy(s);
    return err;
}

int main(void)
{
    int ok = 1;
    printf("  constant density, p = cos(pi x) cos(pi z):\n");
    double prev = NAN;
    for (int n = 32; n <= 128; n *= 2) {
        int it = -1;
        double e = convergence_error(n, &it);
        double order = isnan(prev) ? NAN : log2(prev / e);
        printf("    n=%4d  iters=%4d  max error=%.3e  order=%.2f\n", n, it, e, order);
        if (it < 0 || (!isnan(order) && order < 1.8))
            ok = 0;
        prev = e;
    }

    printf("  density jump 1000:1 at z = 0.5, smooth exact solution:\n");
    const int n = 96;
    Sim *s = make_sim(n, 1.0, 1000.0);
    if (!s)
        return 1;
    double *truth = calloc(s->ntot, sizeof(double));
    if (!truth)
        return 1;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            double x = (i + 0.5) * s->dx, z = (j + 0.5) * s->dz;
            truth[IX(s, i, j)] = is_ref(s, i, j) ? 0.0 : sin(3.0 * x) * cos(2.0 * z) + z;
        }
    poisson_apply(s, truth, s->rhs);
    for (size_t k = 0; k < s->ntot; k++)
        s->p[k] = 0.0;
    int it = poisson_solve(s, s->rhs, s->p);
    double err = 0.0, scale = 0.0;
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            err = fmax(err, fabs(s->p[IX(s, i, j)] - truth[IX(s, i, j)]));
            scale = fmax(scale, fabs(truth[IX(s, i, j)]));
        }
    printf("    n=%4d  iters=%4d  relative error=%.3e\n", n, it, err / scale);
    if (it < 0 || err / scale > 1e-6)
        ok = 0;
    free(truth);
    sim_destroy(s);

    printf(ok ? "Poisson: ok\n" : "Poisson: FAILED\n");
    return ok ? 0 : 1;
}
