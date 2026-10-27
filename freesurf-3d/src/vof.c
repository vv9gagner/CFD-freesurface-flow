#include "vof.h"
#include "bc.h"
#include "util.h"

#include <math.h>
#include <string.h>

#define VOF_EPS 1e-12

double line_area(double n1, double n2, double alpha)
{
    if (n1 < 0.0) {
        alpha -= n1;
        n1 = -n1;
    }
    if (n2 < 0.0) {
        alpha -= n2;
        n2 = -n2;
    }
    if (alpha <= 0.0)
        return 0.0;
    if (alpha >= n1 + n2)
        return 1.0;
    /* triangle / trapezoid / complement of a triangle; the small component is > 0 whenever its branch is reached */
    double a = dmin(n1, n2), b = dmax(n1, n2), area;
    if (alpha <= a)
        area = alpha * alpha / (2.0 * a * b);
    else if (alpha <= b)
        area = (alpha - 0.5 * a) / b;
    else {
        double r = a + b - alpha;
        area = 1.0 - r * r / (2.0 * a * b);
    }
    return clamp01(area);
}

/*
 * Volume for non-negative normals, with na the smallest component:
 *   V = (1/na) * integral_{alpha-na}^{alpha} line_area(nb, nc, s) ds.
 * line_area is piecewise quadratic in s with kinks at 0, min, max and nb+nc, so Simpson's
 * rule on each piece is exact. Unlike the closed-form cubic, this has no cancellation when
 * one normal component is tiny (nearly 2D interfaces). The integral runs over t = alpha - s
 * in [0, na], so the piece widths add up to na exactly even when na << alpha.
 */
static double section_average(double na, double nb, double nc, double alpha)
{
    const double kinks[4] = { nb + nc, dmax(nb, nc), dmin(nb, nc), 0.0 }; /* ascending in t */
    double t[6];
    int n = 0;
    t[n++] = 0.0;
    for (int m = 0; m < 4; m++) {
        double tk = alpha - kinks[m];
        if (tk > t[n - 1] && tk < na)
            t[n++] = tk;
    }
    t[n++] = na;
    double sum = 0.0;
    for (int m = 0; m + 1 < n; m++) {
        double t0 = t[m], t1 = t[m + 1];
        sum += (t1 - t0) * (line_area(nb, nc, alpha - t0) + 4.0 * line_area(nb, nc, alpha - 0.5 * (t0 + t1))
                            + line_area(nb, nc, alpha - t1));
    }
    return clamp01(sum / (6.0 * na));
}

/* Order the (non-negative) components so that n[0] is the smallest. */
static void smallest_first(double n[3])
{
    int a = 0;
    if (n[1] < n[a])
        a = 1;
    if (n[2] < n[a])
        a = 2;
    double t = n[0];
    n[0] = n[a];
    n[a] = t;
}

double plane_volume(double nx, double ny, double nz, double alpha)
{
    double n[3] = { nx, ny, nz };
    for (int d = 0; d < 3; d++)
        if (n[d] < 0.0) {
            alpha -= n[d];
            n[d] = -n[d];
        }
    double sum = n[0] + n[1] + n[2];
    if (alpha <= 0.0)
        return 0.0;
    if (alpha >= sum)
        return 1.0;
    smallest_first(n);
    if (n[0] <= 1e-14 * sum)
        return line_area(n[1], n[2], alpha);
    return section_average(n[0], n[1], n[2], alpha);
}

/* Inverse of line_area for non-negative normals (not both zero). */
static double line_alpha(double c, double n1, double n2)
{
    double s = n1 + n2, a = dmin(n1, n2) / s, b = dmax(n1, n2) / s;
    double v1 = 0.5 * a / b, alpha;
    if (c <= v1)
        alpha = sqrt(2.0 * c * a * b);
    else if (c <= 1.0 - v1)
        alpha = c * b + 0.5 * a;
    else
        alpha = 1.0 - sqrt(2.0 * a * b * (1.0 - c));
    return alpha * s;
}

double plane_alpha(double c, double nx, double ny, double nz)
{
    double n[3] = { nx, ny, nz }, shift = 0.0;
    for (int d = 0; d < 3; d++)
        if (n[d] < 0.0) {
            shift += n[d];
            n[d] = -n[d];
        }
    const double sum = n[0] + n[1] + n[2];
    c = clamp01(c);
    if (c <= 0.0)
        return shift;
    if (c >= 1.0)
        return shift + sum;
    smallest_first(n);
    if (n[0] <= 1e-14 * sum)
        return shift + line_alpha(c, n[1], n[2]);

    /* safeguarded Newton on the monotone V(alpha); dV/dalpha = (A(alpha) - A(alpha - na)) / na */
    double lo = 0.0, hi = sum, x = c * sum;
    for (int it = 0; it < 100; it++) {
        double f = section_average(n[0], n[1], n[2], x) - c;
        if (fabs(f) < 1e-15)
            break;
        if (f > 0.0)
            hi = x;
        else
            lo = x;
        double df = (line_area(n[1], n[2], x) - line_area(n[1], n[2], x - n[0])) / n[0];
        double xn = df > 0.0 ? x - f / df : 0.5 * (lo + hi);
        if (!(xn > lo && xn < hi))
            xn = 0.5 * (lo + hi);
        if (fabs(xn - x) < 1e-16 * sum) {
            x = xn;
            break;
        }
        x = xn;
    }
    return shift + x;
}

double plane_height(double nx, double ny, double nz, double alpha)
{
    double m = fabs(nx) + fabs(ny) + fabs(nz);
    if (fabs(nz) < 1e-6 * m)
        return 0.5;
    return clamp01((alpha - 0.5 * nx - 0.5 * ny) / nz);
}

void plic_normal(const Sim *s, const double *c, int i, int j, int k, double out[3])
{
    const ptrdiff_t st[3] = { 1, s->Sx, s->Sy };
    const double *p = c + IX(s, i, j, k);
    double cand[4][3];

    /* centred columns along each axis d: 3-cell column sums give heights H(e, f) */
    for (int d = 0; d < 3; d++) {
        const int e = (d + 1) % 3, f = (d + 2) % 3;
        double H[3][3], lower = 0.0, upper = 0.0;
        for (int oe = -1; oe <= 1; oe++)
            for (int of = -1; of <= 1; of++) {
                double h = 0.0;
                for (int m = -1; m <= 1; m++)
                    h += p[m * st[d] + oe * st[e] + of * st[f]];
                H[oe + 1][of + 1] = h;
                lower += p[-st[d] + oe * st[e] + of * st[f]];
                upper += p[st[d] + oe * st[e] + of * st[f]];
            }
        cand[d][d] = lower > upper ? 1.0 : -1.0; /* more liquid on the low side: normal points +d */
        cand[d][e] = -0.5 * (H[2][1] - H[0][1]);
        cand[d][f] = -0.5 * (H[1][2] - H[1][0]);
        double m = fabs(cand[d][0]) + fabs(cand[d][1]) + fabs(cand[d][2]);
        for (int a = 0; a < 3; a++)
            cand[d][a] /= m;
    }
    int cn = 0;
    for (int d = 1; d < 3; d++)
        if (fabs(cand[d][d]) > fabs(cand[cn][cn]))
            cn = d;

    /* Youngs: weighted 3x3x3 gradient, n = -grad C */
    static const double w[3] = { 1.0, 2.0, 1.0 };
    double m = 0.0;
    for (int d = 0; d < 3; d++) {
        const int e = (d + 1) % 3, f = (d + 2) % 3;
        double g = 0.0;
        for (int oe = -1; oe <= 1; oe++)
            for (int of = -1; of <= 1; of++)
                g += w[oe + 1] * w[of + 1] * (p[-st[d] + oe * st[e] + of * st[f]] - p[st[d] + oe * st[e] + of * st[f]]);
        cand[3][d] = g;
        m += fabs(g);
    }
    if (m > 1e-30) {
        double t0 = 0.0;
        for (int d = 0; d < 3; d++) {
            cand[3][d] /= m;
            t0 = dmax(t0, fabs(cand[3][d]));
        }
        /* MYC rule (Aulisa et al. 2007): Youngs if it gives the steeper interface */
        if (fabs(cand[cn][cn]) > t0)
            cn = 3;
    }
    for (int d = 0; d < 3; d++)
        out[d] = cand[cn][d];
}

void vof_reconstruct(Sim *s, const double *c)
{
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < s->nz; k++)
        for (int j = 0; j < s->ny; j++)
            for (int i = 0; i < s->nx; i++) {
                size_t n = IX(s, i, j, k);
                if (c[n] > VOF_EPS && c[n] < 1.0 - VOF_EPS) {
                    double nv[3];
                    plic_normal(s, c, i, j, k, nv);
                    s->pn_x[n] = nv[0];
                    s->pn_y[n] = nv[1];
                    s->pn_z[n] = nv[2];
                    s->palpha[n] = plane_alpha(c[n], nv[0], nv[1], nv[2]);
                }
            }
}

/* Liquid volume, as a fraction of the donor cell, inside the slab [x0, x1] along axis d (unit coordinates). */
static double slab(const Sim *s, size_t n, int d, double x0, double x1)
{
    double c = s->C[n];
    if (c <= VOF_EPS)
        return 0.0;
    if (c >= 1.0 - VOF_EPS)
        return x1 - x0;
    double nv[3] = { s->pn_x[n], s->pn_y[n], s->pn_z[n] };
    double alpha = s->palpha[n] - nv[d] * x0;
    nv[d] *= x1 - x0;
    return (x1 - x0) * plane_volume(nv[0], nv[1], nv[2], alpha);
}

/*
 * One directional sweep of the Weymouth & Yue (2010) scheme:
 *   C* = C - (F+ - F-) + c_c (s+ - s-),
 * F the geometric PLIC fluxes, s the face Courant numbers and c_c = [C0 > 1/2] frozen over
 * the step. Exactly conservative for div-free u, bounded for CFL <= 1/2.
 */
static void sweep(Sim *s, int d, double dt)
{
    const int n[3] = { s->nx, s->ny, s->nz };
    const double h[3] = { s->dx, s->dy, s->dz };
    const ptrdiff_t st[3] = { 1, s->Sx, s->Sy };
    const double *vel = d == 0 ? s->u : d == 1 ? s->v : s->w;
    double *C = s->C, *flux = s->flux;
    const int ni = n[0] + (d == 0), nj = n[1] + (d == 1), nk = n[2] + (d == 2);

    vof_reconstruct(s, C);
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < nk; k++)
        for (int j = 0; j < nj; j++)
            for (int i = 0; i < ni; i++) {
                const size_t f = IX(s, i, j, k);
                const int idx = d == 0 ? i : d == 1 ? j : k;
                const double cr = vel[f] * dt / h[d];
                double fl = 0.0;
                if (idx > 0 && idx < n[d]) {
                    if (cr > 0.0)
                        fl = slab(s, f - (size_t)st[d], d, 1.0 - cr, 1.0);
                    else if (cr < 0.0)
                        fl = -slab(s, f, d, 0.0, -cr);
                }
                flux[f] = fl;
            }
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < n[2]; k++)
        for (int j = 0; j < n[1]; j++)
            for (int i = 0; i < n[0]; i++) {
                const size_t c = IX(s, i, j, k), cp = c + (size_t)st[d];
                const double dil = (vel[cp] - vel[c]) * dt / h[d];
                const double cc = s->C0[c] > 0.5 ? 1.0 : 0.0;
                C[c] = clamp01(C[c] - (flux[cp] - flux[c]) + cc * dil);
            }
    bc_scalar(s, C);
}

void vof_advect(Sim *s, double dt)
{
    memcpy(s->C0, s->C, s->ntot * sizeof(double));
    /* cyclic sweep order xyz, yzx, zxy to avoid a directional bias */
    const int first = (int)(s->step % 3);
    for (int m = 0; m < 3; m++)
        sweep(s, (first + m) % 3, dt);
}

void vof_update_properties(Sim *s)
{
    const Params *p = &s->prm;
#pragma omp parallel for schedule(static)
    for (size_t n = 0; n < s->ntot; n++) {
        double c = clamp01(s->C[n]);
        s->rho[n] = c * p->rho_l + (1.0 - c) * p->rho_g;
        s->mu[n] = c * p->mu_l + (1.0 - c) * p->mu_g;
    }
}
