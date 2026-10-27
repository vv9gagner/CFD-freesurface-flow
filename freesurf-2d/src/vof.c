#include "vof.h"
#include "bc.h"
#include "util.h"

#include <math.h>
#include <string.h>

#define VOF_EPS 1e-12

double plic_area(double nx, double nz, double alpha)
{
    /* reflect so both normal components are non-negative */
    if (nx < 0.0) {
        alpha -= nx;
        nx = -nx;
    }
    if (nz < 0.0) {
        alpha -= nz;
        nz = -nz;
    }
    if (alpha <= 0.0)
        return 0.0;
    if (alpha >= nx + nz)
        return 1.0;

    /* piecewise form (triangle / trapezoid / complement of a triangle) avoids cancellation
     * for nearly axis-aligned lines; n1 > 0 whenever its branch is reached */
    double n1 = dmin(nx, nz), n2 = dmax(nx, nz), area;
    if (alpha <= n1)
        area = alpha * alpha / (2.0 * n1 * n2);
    else if (alpha <= n2)
        area = (alpha - 0.5 * n1) / n2;
    else {
        double b = n1 + n2 - alpha;
        area = 1.0 - b * b / (2.0 * n1 * n2);
    }
    return clamp01(area);
}

double plic_alpha(double c, double nx, double nz)
{
    double n1 = fabs(nx), n2 = fabs(nz);
    if (n1 > n2) {
        double t = n1;
        n1 = n2;
        n2 = t;
    }
    c = clamp01(c);
    double v1 = 0.5 * n1 / n2; /* fraction at which the line reaches the far corner of the short side */
    double alpha;
    if (c <= v1)
        alpha = sqrt(2.0 * c * n1 * n2);
    else if (c <= 1.0 - v1)
        alpha = c * n2 + 0.5 * n1;
    else
        alpha = n1 + n2 - sqrt(2.0 * n1 * n2 * (1.0 - c));

    /* undo the reflection used for negative components */
    if (nx < 0.0)
        alpha += nx;
    if (nz < 0.0)
        alpha += nz;
    return alpha;
}

void plic_normal(const Sim *s, const double *c, int i, int j, double *pnx, double *pnz)
{
#define CC(a, b) c[IX(s, i + (a), j + (b))]
    /* centred-columns estimate: heights from 3x3 column/row sums */
    double c_t = CC(-1, 1) + CC(0, 1) + CC(1, 1);
    double c_b = CC(-1, -1) + CC(0, -1) + CC(1, -1);
    double c_r = CC(1, -1) + CC(1, 0) + CC(1, 1);
    double c_l = CC(-1, -1) + CC(-1, 0) + CC(-1, 1);
    double mx0 = 0.5 * (c_l - c_r);
    double mz0 = 0.5 * (c_b - c_t);
    int horizontal; /* interface closer to horizontal: use heights along z */
    if (fabs(mx0) <= fabs(mz0)) {
        mz0 = mz0 > 0.0 ? 1.0 : -1.0;
        horizontal = 1;
    } else {
        mx0 = mx0 > 0.0 ? 1.0 : -1.0;
        horizontal = 0;
    }

    /* Youngs' estimate: weighted 3x3 gradient */
    double mx1 = (CC(-1, -1) + 2.0 * CC(-1, 0) + CC(-1, 1)) - (CC(1, -1) + 2.0 * CC(1, 0) + CC(1, 1)) + 1e-30;
    double mz1 = (CC(-1, -1) + 2.0 * CC(0, -1) + CC(1, -1)) - (CC(-1, 1) + 2.0 * CC(0, 1) + CC(1, 1)) + 1e-30;
#undef CC

    /* MYC selection rule (Aulisa et al. 2007): take Youngs' normal if its slope is steeper */
    if (horizontal) {
        if (fabs(mx1) / fabs(mz1) > fabs(mx0)) {
            mx0 = mx1;
            mz0 = mz1;
        }
    } else {
        if (fabs(mz1) / fabs(mx1) > fabs(mz0)) {
            mx0 = mx1;
            mz0 = mz1;
        }
    }
    double m = fabs(mx0) + fabs(mz0);
    *pnx = mx0 / m;
    *pnz = mz0 / m;
}

int plic_segment(double nx, double nz, double alpha, double pts[4])
{
    const double tol = 1e-12;
    double cx[4], cz[4];
    int n = 0;
    if (fabs(nz) > tol) {
        double z0 = alpha / nz, z1 = (alpha - nx) / nz; /* at x = 0 and x = 1 */
        if (z0 >= -tol && z0 <= 1.0 + tol) {
            cx[n] = 0.0;
            cz[n++] = clamp01(z0);
        }
        if (z1 >= -tol && z1 <= 1.0 + tol) {
            cx[n] = 1.0;
            cz[n++] = clamp01(z1);
        }
    }
    if (fabs(nx) > tol) {
        double x0 = alpha / nx, x1 = (alpha - nz) / nx; /* at z = 0 and z = 1 */
        if (x0 >= -tol && x0 <= 1.0 + tol) {
            cx[n] = clamp01(x0);
            cz[n++] = 0.0;
        }
        if (x1 >= -tol && x1 <= 1.0 + tol) {
            cx[n] = clamp01(x1);
            cz[n++] = 1.0;
        }
    }
    /* corners can be found twice: keep the most distant pair */
    double best = -1.0;
    for (int a = 0; a < n; a++)
        for (int b = a + 1; b < n; b++) {
            double d = sq(cx[a] - cx[b]) + sq(cz[a] - cz[b]);
            if (d > best) {
                best = d;
                pts[0] = cx[a];
                pts[1] = cz[a];
                pts[2] = cx[b];
                pts[3] = cz[b];
            }
        }
    return best > 1e-20;
}

void vof_reconstruct(Sim *s, const double *c)
{
    for (int j = 0; j < s->nz; j++)
        for (int i = 0; i < s->nx; i++) {
            size_t k = IX(s, i, j);
            if (c[k] > VOF_EPS && c[k] < 1.0 - VOF_EPS) {
                double nx, nz;
                plic_normal(s, c, i, j, &nx, &nz);
                s->pn_x[k] = nx;
                s->pn_z[k] = nz;
                s->palpha[k] = plic_alpha(c[k], nx, nz);
            }
        }
}

/* Liquid volume, as a fraction of the donor cell, inside the slab x in [x0, x1] (unit coordinates). */
static double slab_x(const Sim *s, size_t k, double c, double x0, double x1)
{
    if (c <= VOF_EPS)
        return 0.0;
    if (c >= 1.0 - VOF_EPS)
        return x1 - x0;
    double nx = s->pn_x[k], nz = s->pn_z[k];
    return (x1 - x0) * plic_area(nx * (x1 - x0), nz, s->palpha[k] - nx * x0);
}

static double slab_z(const Sim *s, size_t k, double c, double z0, double z1)
{
    if (c <= VOF_EPS)
        return 0.0;
    if (c >= 1.0 - VOF_EPS)
        return z1 - z0;
    double nx = s->pn_x[k], nz = s->pn_z[k];
    return (z1 - z0) * plic_area(nx, nz * (z1 - z0), s->palpha[k] - nz * z0);
}

/*
 * One directional sweep of the Weymouth & Yue (2010) scheme:
 *   C* = C - (F_{+} - F_{-}) + c_c (s_{+} - s_{-}),
 * with F the geometric PLIC fluxes, s the face Courant numbers and c_c = [C0 > 1/2]
 * frozen over the whole step. Exactly conservative for div-free u, and bounded for CFL <= 1/2.
 */
static void sweep_x(Sim *s, double dt)
{
    const int nx = s->nx, nz = s->nz;
    double *C = s->C;
    vof_reconstruct(s, C);
    for (int j = 0; j < nz; j++)
        for (int i = 0; i <= nx; i++) {
            double cr = s->u[IX(s, i, j)] * dt / s->dx;
            double f = 0.0;
            if (i > 0 && i < nx) {
                if (cr > 0.0) {
                    size_t kd = IX(s, i - 1, j);
                    f = slab_x(s, kd, C[kd], 1.0 - cr, 1.0);
                } else if (cr < 0.0) {
                    size_t kd = IX(s, i, j);
                    f = -slab_x(s, kd, C[kd], 0.0, -cr);
                }
            }
            s->flux[IX(s, i, j)] = f;
        }
    for (int j = 0; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t k = IX(s, i, j), kr = IX(s, i + 1, j);
            double dil = (s->u[kr] - s->u[k]) * dt / s->dx;
            double cc = s->C0[k] > 0.5 ? 1.0 : 0.0;
            C[k] = clamp01(C[k] - (s->flux[kr] - s->flux[k]) + cc * dil);
        }
    bc_scalar(s, C);
}

static void sweep_z(Sim *s, double dt)
{
    const int nx = s->nx, nz = s->nz;
    double *C = s->C;
    vof_reconstruct(s, C);
    for (int j = 0; j <= nz; j++)
        for (int i = 0; i < nx; i++) {
            double cr = s->w[IX(s, i, j)] * dt / s->dz;
            double f = 0.0;
            if (j > 0 && j < nz) {
                if (cr > 0.0) {
                    size_t kd = IX(s, i, j - 1);
                    f = slab_z(s, kd, C[kd], 1.0 - cr, 1.0);
                } else if (cr < 0.0) {
                    size_t kd = IX(s, i, j);
                    f = -slab_z(s, kd, C[kd], 0.0, -cr);
                }
            }
            s->flux[IX(s, i, j)] = f;
        }
    for (int j = 0; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t k = IX(s, i, j), kt = IX(s, i, j + 1);
            double dil = (s->w[kt] - s->w[k]) * dt / s->dz;
            double cc = s->C0[k] > 0.5 ? 1.0 : 0.0;
            C[k] = clamp01(C[k] - (s->flux[kt] - s->flux[k]) + cc * dil);
        }
    bc_scalar(s, C);
}

void vof_advect(Sim *s, double dt)
{
    memcpy(s->C0, s->C, s->ntot * sizeof(double));
    /* alternate the sweep order (Strang-like) to avoid a directional bias */
    if (s->step % 2 == 0) {
        sweep_x(s, dt);
        sweep_z(s, dt);
    } else {
        sweep_z(s, dt);
        sweep_x(s, dt);
    }
}

void vof_update_properties(Sim *s)
{
    const Params *p = &s->prm;
    for (size_t k = 0; k < s->ntot; k++) {
        double c = clamp01(s->C[k]);
        s->rho[k] = c * p->rho_l + (1.0 - c) * p->rho_g;
        s->mu[k] = c * p->mu_l + (1.0 - c) * p->mu_g;
    }
}
