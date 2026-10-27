#include "poisson.h"
#include "bc.h"
#include "util.h"

#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static inline int is_ref(const Sim *s, int i, int j, int k)
{
    return i == s->ref_i && j == s->ref_j && k == s->ref_k;
}

/* 1/(rho h^2) on the face between cell n - o and cell n, with an arithmetic face density */
static inline double beta(const Sim *s, size_t n, ptrdiff_t o, double h)
{
    return 2.0 / ((s->rho[n - o] + s->rho[n]) * h * h);
}

/*
 * Block-Jacobi MIC(0) preconditioner. The domain is cut into pc_blocks slabs along its
 * longest axis; the incomplete factorisation drops the couplings between slabs, so every
 * slab is factored and solved independently, one per thread. With one block this is the
 * plain MIC(0) of Bridson, "Fluid Simulation for Computer Graphics". The preconditioner
 * stays symmetric positive definite, so PCG converges to the same solution.
 */
typedef struct {
    int lo[3], hi[3]; /* cell range of the block */
} Block;

static Block block_of(const Sim *s, int b)
{
    Block bl = { { 0, 0, 0 }, { s->nx, s->ny, s->nz } };
    const int a = s->pc_axis, n = bl.hi[a];
    bl.lo[a] = (int)((long)n * b / s->pc_blocks);
    bl.hi[a] = (int)((long)n * (b + 1) / s->pc_blocks);
    return bl;
}

static void choose_blocks(Sim *s)
{
    const int n[3] = { s->nx, s->ny, s->nz };
    int a = 0;
    for (int d = 1; d < 3; d++)
        if (n[d] > n[a])
            a = d;
    int nb = 1;
#ifdef _OPENMP
    nb = omp_get_max_threads();
#endif
    if (nb > n[a] / 4) /* keep slabs at least 4 cells thick */
        nb = n[a] / 4;
    s->pc_axis = a;
    s->pc_blocks = nb < 1 ? 1 : nb;
}

/* Factorisation of one block: couplings that leave the block are treated as zero. */
static void factor_block(Sim *s, const Block *bl)
{
    const ptrdiff_t X = 1, Y = s->Sx, Z = s->Sy;
    const double tune = 0.97, safety = 0.25;
    const double *Ax = s->Ax, *Ay = s->Ay, *Az = s->Az;
    double *P = s->precon;
    for (int k = bl->lo[2]; k < bl->hi[2]; k++)
        for (int j = bl->lo[1]; j < bl->hi[1]; j++)
            for (int i = bl->lo[0]; i < bl->hi[0]; i++) {
                const size_t n = IX(s, i, j, k);
                if (is_ref(s, i, j, k)) {
                    P[n] = 0.0;
                    continue;
                }
                /* lower neighbours inside the block, and whether n's upper neighbours are */
                const int dn[3] = { i > bl->lo[0], j > bl->lo[1], k > bl->lo[2] };
                const int up[3] = { i < bl->hi[0] - 1, j < bl->hi[1] - 1, k < bl->hi[2] - 1 };
                const double ax = dn[0] ? Ax[n - X] : 0.0, px = dn[0] ? P[n - X] : 0.0;
                const double ay = dn[1] ? Ay[n - Y] : 0.0, py = dn[1] ? P[n - Y] : 0.0;
                const double az = dn[2] ? Az[n - Z] : 0.0, pz = dn[2] ? P[n - Z] : 0.0;
                double fill = 0.0;
                if (dn[0])
                    fill += ax * ((up[1] ? Ay[n - X] : 0.0) + (up[2] ? Az[n - X] : 0.0)) * px * px;
                if (dn[1])
                    fill += ay * ((up[0] ? Ax[n - Y] : 0.0) + (up[2] ? Az[n - Y] : 0.0)) * py * py;
                if (dn[2])
                    fill += az * ((up[0] ? Ax[n - Z] : 0.0) + (up[1] ? Ay[n - Z] : 0.0)) * pz * pz;
                double e = s->Adiag[n] - sq(ax * px) - sq(ay * py) - sq(az * pz) - tune * fill;
                if (e < safety * s->Adiag[n])
                    e = s->Adiag[n];
                P[n] = 1.0 / sqrt(e);
            }
}

void poisson_setup(Sim *s)
{
    const int nx = s->nx, ny = s->ny, nz = s->nz;
    const ptrdiff_t X = 1, Y = s->Sx, Z = s->Sy;

#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                size_t n = IX(s, i, j, k);
                if (is_ref(s, i, j, k)) {
                    s->Adiag[n] = s->Ax[n] = s->Ay[n] = s->Az[n] = 0.0;
                    continue;
                }
                double d = 0.0;
                if (i > 0)
                    d += beta(s, n, X, s->dx);
                if (i < nx - 1)
                    d += beta(s, n + X, X, s->dx);
                if (j > 0)
                    d += beta(s, n, Y, s->dy);
                if (j < ny - 1)
                    d += beta(s, n + Y, Y, s->dy);
                if (k > 0)
                    d += beta(s, n, Z, s->dz);
                if (k < nz - 1)
                    d += beta(s, n + Z, Z, s->dz);
                s->Adiag[n] = d;
                /* couplings to the pinned cell vanish (p_ref = 0 acts as a Dirichlet value) */
                s->Ax[n] = (i < nx - 1 && !is_ref(s, i + 1, j, k)) ? -beta(s, n + X, X, s->dx) : 0.0;
                s->Ay[n] = (j < ny - 1 && !is_ref(s, i, j + 1, k)) ? -beta(s, n + Y, Y, s->dy) : 0.0;
                s->Az[n] = (k < nz - 1 && !is_ref(s, i, j, k + 1)) ? -beta(s, n + Z, Z, s->dz) : 0.0;
            }

    choose_blocks(s);
#pragma omp parallel for schedule(static, 1)
    for (int b = 0; b < s->pc_blocks; b++) {
        Block bl = block_of(s, b);
        factor_block(s, &bl);
    }
}

void poisson_apply(const Sim *s, const double *x, double *y)
{
    const ptrdiff_t X = 1, Y = s->Sx, Z = s->Sy;
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < s->nz; k++)
        for (int j = 0; j < s->ny; j++)
            for (int i = 0; i < s->nx; i++) {
                size_t n = IX(s, i, j, k);
                if (is_ref(s, i, j, k)) {
                    y[n] = 0.0;
                    continue;
                }
                y[n] = s->Adiag[n] * x[n] + s->Ax[n] * x[n + X] + s->Ax[n - X] * x[n - X] + s->Ay[n] * x[n + Y]
                     + s->Ay[n - Y] * x[n - Y] + s->Az[n] * x[n + Z] + s->Az[n - Z] * x[n - Z];
            }
}

/* z = M^{-1} r for one block (forward then backward substitution); returns the block's part of z.r */
static double precon_block(Sim *s, const Block *bl, const double *r, double *z)
{
    const ptrdiff_t X = 1, Y = s->Sx, Z = s->Sy;
    const double *P = s->precon, *Ax = s->Ax, *Ay = s->Ay, *Az = s->Az;
    double *q = s->tmp;
    for (int k = bl->lo[2]; k < bl->hi[2]; k++)
        for (int j = bl->lo[1]; j < bl->hi[1]; j++)
            for (int i = bl->lo[0]; i < bl->hi[0]; i++) {
                const size_t n = IX(s, i, j, k);
                if (is_ref(s, i, j, k)) {
                    q[n] = 0.0;
                    continue;
                }
                double t = r[n];
                if (i > bl->lo[0])
                    t -= Ax[n - X] * P[n - X] * q[n - X];
                if (j > bl->lo[1])
                    t -= Ay[n - Y] * P[n - Y] * q[n - Y];
                if (k > bl->lo[2])
                    t -= Az[n - Z] * P[n - Z] * q[n - Z];
                q[n] = t * P[n];
            }
    double zr = 0.0;
    for (int k = bl->hi[2] - 1; k >= bl->lo[2]; k--)
        for (int j = bl->hi[1] - 1; j >= bl->lo[1]; j--)
            for (int i = bl->hi[0] - 1; i >= bl->lo[0]; i--) {
                const size_t n = IX(s, i, j, k);
                if (is_ref(s, i, j, k)) {
                    z[n] = 0.0;
                    continue;
                }
                double t = 0.0;
                if (i < bl->hi[0] - 1)
                    t += Ax[n] * z[n + X];
                if (j < bl->hi[1] - 1)
                    t += Ay[n] * z[n + Y];
                if (k < bl->hi[2] - 1)
                    t += Az[n] * z[n + Z];
                z[n] = (q[n] - P[n] * t) * P[n];
                zr += z[n] * r[n];
            }
    return zr;
}

/* z = M^{-1} r over all blocks in parallel; returns z.r */
static double apply_precon(Sim *s, const double *r, double *z)
{
    double zr = 0.0;
#pragma omp parallel for schedule(static, 1) reduction(+ : zr)
    for (int b = 0; b < s->pc_blocks; b++) {
        Block bl = block_of(s, b);
        zr += precon_block(s, &bl, r, z);
    }
    return zr;
}

int poisson_solve(Sim *s, const double *b, double *x)
{
    const ptrdiff_t X = 1, Y = s->Sx, Z = s->Sy;
    double *r = s->res, *z = s->zv, *d = s->sv, *q = s->qv;
    if (s->pc_blocks < 1)
        choose_blocks(s);

    x[IX(s, s->ref_i, s->ref_j, s->ref_k)] = 0.0;
    double bnorm = 0.0;
#pragma omp parallel for collapse(2) reduction(max : bnorm) schedule(static)
    for (int k = 0; k < s->nz; k++)
        for (int j = 0; j < s->ny; j++)
            for (int i = 0; i < s->nx; i++)
                bnorm = dmax(bnorm, fabs(b[IX(s, i, j, k)]));
    if (bnorm == 0.0) {
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < s->nz; k++)
            for (int j = 0; j < s->ny; j++)
                for (int i = 0; i < s->nx; i++)
                    x[IX(s, i, j, k)] = 0.0;
        s->pois_res = 0.0;
        return 0;
    }
    const double thresh = s->prm.poisson_tol * bnorm;

    poisson_apply(s, x, q);
    double rn = 0.0;
#pragma omp parallel for collapse(2) reduction(max : rn) schedule(static)
    for (int k = 0; k < s->nz; k++)
        for (int j = 0; j < s->ny; j++)
            for (int i = 0; i < s->nx; i++) {
                size_t n = IX(s, i, j, k);
                r[n] = is_ref(s, i, j, k) ? 0.0 : b[n] - q[n];
                rn = dmax(rn, fabs(r[n]));
            }
    if (rn <= thresh) {
        s->pois_res = rn / bnorm;
        return 0;
    }

    double rho_old = apply_precon(s, r, z);
#pragma omp parallel for schedule(static)
    for (size_t n = 0; n < s->ntot; n++)
        d[n] = z[n];

    for (int it = 1; it <= s->prm.poisson_maxit; it++) {
        /* q = A d together with d.q */
        double dq = 0.0;
#pragma omp parallel for collapse(2) reduction(+ : dq) schedule(static)
        for (int k = 0; k < s->nz; k++)
            for (int j = 0; j < s->ny; j++)
                for (int i = 0; i < s->nx; i++) {
                    size_t n = IX(s, i, j, k);
                    double v = is_ref(s, i, j, k) ? 0.0
                             : s->Adiag[n] * d[n] + s->Ax[n] * d[n + X] + s->Ax[n - X] * d[n - X] + s->Ay[n] * d[n + Y]
                                   + s->Ay[n - Y] * d[n - Y] + s->Az[n] * d[n + Z] + s->Az[n - Z] * d[n - Z];
                    q[n] = v;
                    dq += d[n] * v;
                }
        const double alpha = rho_old / dq;
        rn = 0.0;
#pragma omp parallel for collapse(2) reduction(max : rn) schedule(static)
        for (int k = 0; k < s->nz; k++)
            for (int j = 0; j < s->ny; j++)
                for (int i = 0; i < s->nx; i++) {
                    size_t n = IX(s, i, j, k);
                    x[n] += alpha * d[n];
                    r[n] -= alpha * q[n];
                    rn = dmax(rn, fabs(r[n]));
                }
        if (rn <= thresh) {
            s->pois_res = rn / bnorm;
            return it;
        }
        const double rho_new = apply_precon(s, r, z);
        const double beta_cg = rho_new / rho_old;
        rho_old = rho_new;
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < s->nz; k++)
            for (int j = 0; j < s->ny; j++)
                for (int i = 0; i < s->nx; i++) {
                    size_t n = IX(s, i, j, k);
                    d[n] = z[n] + beta_cg * d[n];
                }
    }
    s->pois_res = rn / bnorm;
    return -1;
}

int projection(Sim *s, double dt)
{
    const int nx = s->nx, ny = s->ny, nz = s->nz;
    const double dx = s->dx, dy = s->dy, dz = s->dz;
    const ptrdiff_t X = 1, Y = s->Sx, Z = s->Sy;

    poisson_setup(s);
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                size_t n = IX(s, i, j, k);
                double div = (s->us[n + X] - s->us[n]) / dx + (s->vs[n + Y] - s->vs[n]) / dy + (s->ws[n + Z] - s->ws[n]) / dz;
                s->rhs[n] = is_ref(s, i, j, k) ? 0.0 : -div / dt;
            }

    int iters = poisson_solve(s, s->rhs, s->p);
    s->pois_iters = iters < 0 ? s->prm.poisson_maxit : iters;

    /* u = u* - dt/rho_f grad p on interior faces; wall faces stay zero */
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k <= nz; k++)
        for (int j = 0; j <= ny; j++)
            for (int i = 0; i <= nx; i++) {
                size_t n = IX(s, i, j, k);
                const double *p = s->p;
                if (j < ny && k < nz)
                    s->u[n] = (i > 0 && i < nx) ? s->us[n] - dt * beta(s, n, X, dx) * dx * (p[n] - p[n - X]) : 0.0;
                if (i < nx && k < nz)
                    s->v[n] = (j > 0 && j < ny) ? s->vs[n] - dt * beta(s, n, Y, dy) * dy * (p[n] - p[n - Y]) : 0.0;
                if (i < nx && j < ny)
                    s->w[n] = (k > 0 && k < nz) ? s->ws[n] - dt * beta(s, n, Z, dz) * dz * (p[n] - p[n - Z]) : 0.0;
            }
    bc_u(s, s->u);
    bc_v(s, s->v);
    bc_w(s, s->w);
    return iters;
}
