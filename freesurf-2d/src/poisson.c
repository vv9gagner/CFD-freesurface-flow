#include "poisson.h"
#include "bc.h"
#include "util.h"

#include <math.h>

static inline int is_ref(const Sim *s, int i, int j)
{
    return i == s->ref_i && j == s->ref_j;
}

/* 1/(rho dx^2) on x-face i and 1/(rho dz^2) on z-face j, with arithmetic face densities */
static inline double beta_x(const Sim *s, int i, int j)
{
    return 2.0 / ((s->rho[IX(s, i - 1, j)] + s->rho[IX(s, i, j)]) * s->dx * s->dx);
}

static inline double beta_z(const Sim *s, int i, int j)
{
    return 2.0 / ((s->rho[IX(s, i, j - 1)] + s->rho[IX(s, i, j)]) * s->dz * s->dz);
}

void poisson_setup(Sim *s)
{
    const int nx = s->nx, nz = s->nz;
    for (int j = 0; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t k = IX(s, i, j);
            if (is_ref(s, i, j)) {
                s->Adiag[k] = s->Ax[k] = s->Az[k] = 0.0;
                continue;
            }
            double d = 0.0;
            if (i > 0)
                d += beta_x(s, i, j);
            if (i < nx - 1)
                d += beta_x(s, i + 1, j);
            if (j > 0)
                d += beta_z(s, i, j);
            if (j < nz - 1)
                d += beta_z(s, i, j + 1);
            s->Adiag[k] = d;
            /* couplings to the pinned cell vanish (p_ref = 0 acts as a Dirichlet value) */
            s->Ax[k] = (i < nx - 1 && !is_ref(s, i + 1, j)) ? -beta_x(s, i + 1, j) : 0.0;
            s->Az[k] = (j < nz - 1 && !is_ref(s, i, j + 1)) ? -beta_z(s, i, j + 1) : 0.0;
        }

    /* Modified incomplete Cholesky, MIC(0), as in Bridson, "Fluid Simulation for Computer Graphics" */
    const double tune = 0.97, safety = 0.25;
    double *P = s->precon;
    for (int j = 0; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t k = IX(s, i, j), kl = IX(s, i - 1, j), kb = IX(s, i, j - 1);
            if (is_ref(s, i, j)) {
                P[k] = 0.0;
                continue;
            }
            double e = s->Adiag[k]
                     - sq(s->Ax[kl] * P[kl]) - sq(s->Az[kb] * P[kb])
                     - tune * (s->Ax[kl] * s->Az[kl] * sq(P[kl]) + s->Az[kb] * s->Ax[kb] * sq(P[kb]));
            if (e < safety * s->Adiag[k])
                e = s->Adiag[k];
            P[k] = 1.0 / sqrt(e);
        }
}

void poisson_apply(const Sim *s, const double *x, double *y)
{
    for (int j = 0; j < s->nz; j++)
        for (int i = 0; i < s->nx; i++) {
            size_t k = IX(s, i, j);
            if (is_ref(s, i, j)) {
                y[k] = 0.0;
                continue;
            }
            size_t kl = IX(s, i - 1, j), kr = IX(s, i + 1, j), kb = IX(s, i, j - 1), kt = IX(s, i, j + 1);
            y[k] = s->Adiag[k] * x[k] + s->Ax[k] * x[kr] + s->Ax[kl] * x[kl] + s->Az[k] * x[kt] + s->Az[kb] * x[kb];
        }
}

/* z = M^{-1} r with M = L L^T from MIC(0); s->tmp holds the intermediate solve. */
static void apply_precon(Sim *s, const double *r, double *z)
{
    const int nx = s->nx, nz = s->nz;
    const double *P = s->precon;
    double *q = s->tmp;
    for (int j = 0; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t k = IX(s, i, j), kl = IX(s, i - 1, j), kb = IX(s, i, j - 1);
            if (is_ref(s, i, j)) {
                q[k] = 0.0;
                continue;
            }
            double t = r[k] - s->Ax[kl] * P[kl] * q[kl] - s->Az[kb] * P[kb] * q[kb];
            q[k] = t * P[k];
        }
    for (int j = nz - 1; j >= 0; j--)
        for (int i = nx - 1; i >= 0; i--) {
            size_t k = IX(s, i, j), kr = IX(s, i + 1, j), kt = IX(s, i, j + 1);
            if (is_ref(s, i, j)) {
                z[k] = 0.0;
                continue;
            }
            double t = q[k] - s->Ax[k] * P[k] * z[kr] - s->Az[k] * P[k] * z[kt];
            z[k] = t * P[k];
        }
}

static double dot(const Sim *s, const double *a, const double *b)
{
    double sum = 0.0;
    for (int j = 0; j < s->nz; j++)
        for (int i = 0; i < s->nx; i++) {
            size_t k = IX(s, i, j);
            sum += a[k] * b[k];
        }
    return sum;
}

static double norm_inf(const Sim *s, const double *a)
{
    double m = 0.0;
    for (int j = 0; j < s->nz; j++)
        for (int i = 0; i < s->nx; i++)
            m = dmax(m, fabs(a[IX(s, i, j)]));
    return m;
}

int poisson_solve(Sim *s, const double *b, double *x)
{
    const int nx = s->nx, nz = s->nz;
    double *r = s->res, *z = s->zv, *d = s->sv, *q = s->qv;

    x[IX(s, s->ref_i, s->ref_j)] = 0.0;
    double bnorm = norm_inf(s, b);
    if (bnorm == 0.0) {
        for (int j = 0; j < nz; j++)
            for (int i = 0; i < nx; i++)
                x[IX(s, i, j)] = 0.0;
        s->pois_res = 0.0;
        return 0;
    }
    const double thresh = s->prm.poisson_tol * bnorm;

    poisson_apply(s, x, q);
    for (int j = 0; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t k = IX(s, i, j);
            r[k] = is_ref(s, i, j) ? 0.0 : b[k] - q[k];
        }
    double rn = norm_inf(s, r);
    if (rn <= thresh) {
        s->pois_res = rn / bnorm;
        return 0;
    }

    apply_precon(s, r, z);
    for (size_t k = 0; k < s->ntot; k++)
        d[k] = z[k];
    double rho_old = dot(s, z, r);

    for (int it = 1; it <= s->prm.poisson_maxit; it++) {
        poisson_apply(s, d, q);
        double alpha = rho_old / dot(s, d, q);
        for (int j = 0; j < nz; j++)
            for (int i = 0; i < nx; i++) {
                size_t k = IX(s, i, j);
                x[k] += alpha * d[k];
                r[k] -= alpha * q[k];
            }
        rn = norm_inf(s, r);
        if (rn <= thresh) {
            s->pois_res = rn / bnorm;
            return it;
        }
        apply_precon(s, r, z);
        double rho_new = dot(s, z, r);
        double beta = rho_new / rho_old;
        rho_old = rho_new;
        for (int j = 0; j < nz; j++)
            for (int i = 0; i < nx; i++) {
                size_t k = IX(s, i, j);
                d[k] = z[k] + beta * d[k];
            }
    }
    s->pois_res = rn / bnorm;
    return -1;
}

int projection(Sim *s, double dt)
{
    const int nx = s->nx, nz = s->nz;
    const double dx = s->dx, dz = s->dz;

    poisson_setup(s);
    for (int j = 0; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t k = IX(s, i, j);
            double div = (s->us[IX(s, i + 1, j)] - s->us[k]) / dx + (s->ws[IX(s, i, j + 1)] - s->ws[k]) / dz;
            s->rhs[k] = is_ref(s, i, j) ? 0.0 : -div / dt;
        }

    int iters = poisson_solve(s, s->rhs, s->p);
    s->pois_iters = iters < 0 ? s->prm.poisson_maxit : iters;

    for (int j = 0; j < nz; j++) {
        s->u[IX(s, 0, j)] = 0.0;
        s->u[IX(s, nx, j)] = 0.0;
        for (int i = 1; i < nx; i++) {
            size_t k = IX(s, i, j);
            s->u[k] = s->us[k] - dt * beta_x(s, i, j) * dx * (s->p[k] - s->p[IX(s, i - 1, j)]);
        }
    }
    for (int i = 0; i < nx; i++) {
        s->w[IX(s, i, 0)] = 0.0;
        s->w[IX(s, i, nz)] = 0.0;
    }
    for (int j = 1; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t k = IX(s, i, j);
            s->w[k] = s->ws[k] - dt * beta_z(s, i, j) * dz * (s->p[k] - s->p[IX(s, i, j - 1)]);
        }
    bc_u(s, s->u);
    bc_w(s, s->w);
    return iters;
}
