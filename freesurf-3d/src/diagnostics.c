#include "diagnostics.h"
#include "init.h"
#include "util.h"

#include <math.h>
#include <stdlib.h>

double analytic_period(const Params *p)
{
    double kx = p->mode_x * FS_PI / p->Lx, ky = p->mode_y * FS_PI / p->Ly;
    double k = sqrt(kx * kx + ky * ky);
    double num = p->gravity * k * (p->rho_l - p->rho_g);
    double den = p->rho_l / tanh(k * p->depth) + p->rho_g / tanh(k * (p->Lz - p->depth));
    if (!(num > 0.0))
        return NAN;
    return 2.0 * FS_PI / sqrt(num / den);
}

double diag_column_height(const Sim *s, int i, int j)
{
    double h = 0.0;
    for (int k = 0; k < s->nz; k++)
        h += s->C[IX(s, i, j, k)];
    return h * s->dz;
}

double diag_volume(const Sim *s)
{
    double v = 0.0;
#pragma omp parallel for collapse(2) reduction(+ : v) schedule(static)
    for (int k = 0; k < s->nz; k++)
        for (int j = 0; j < s->ny; j++)
            for (int i = 0; i < s->nx; i++)
                v += s->C[IX(s, i, j, k)];
    return v * s->dx * s->dy * s->dz;
}

double diag_kinetic(const Sim *s)
{
    const ptrdiff_t X = 1, Y = s->Sx, Z = s->Sy;
    double e = 0.0;
#pragma omp parallel for collapse(2) reduction(+ : e) schedule(static)
    for (int k = 0; k < s->nz; k++)
        for (int j = 0; j < s->ny; j++)
            for (int i = 0; i < s->nx; i++) {
                size_t n = IX(s, i, j, k);
                /* each cell owns its low faces; wall faces carry no normal velocity */
                if (i > 0)
                    e += 0.25 * (s->rho[n - X] + s->rho[n]) * sq(s->u[n]);
                if (j > 0)
                    e += 0.25 * (s->rho[n - Y] + s->rho[n]) * sq(s->v[n]);
                if (k > 0)
                    e += 0.25 * (s->rho[n - Z] + s->rho[n]) * sq(s->w[n]);
            }
    return e * s->dx * s->dy * s->dz;
}

/* Potential energy of one cell per g dV, with the liquid fraction c at the bottom of the
 * cell and the gas on top (exact for a horizontal interface). */
static double cell_pe(const Params *p, int k, double dz, double c)
{
    double zb = k * dz;
    return p->rho_l * c * (zb + 0.5 * c * dz) + p->rho_g * (1.0 - c) * (zb + 0.5 * (1.0 + c) * dz);
}

double diag_potential(const Sim *s, double pe_flat)
{
    double e = 0.0;
#pragma omp parallel for collapse(2) reduction(+ : e) schedule(static)
    for (int k = 0; k < s->nz; k++)
        for (int j = 0; j < s->ny; j++)
            for (int i = 0; i < s->nx; i++)
                e += cell_pe(&s->prm, k, s->dz, s->C[IX(s, i, j, k)]);
    return e * s->prm.gravity * s->dx * s->dy * s->dz - pe_flat;
}

double diag_max_div(const Sim *s)
{
    const ptrdiff_t X = 1, Y = s->Sx, Z = s->Sy;
    double m = 0.0;
#pragma omp parallel for collapse(2) reduction(max : m) schedule(static)
    for (int k = 0; k < s->nz; k++)
        for (int j = 0; j < s->ny; j++)
            for (int i = 0; i < s->nx; i++) {
                size_t n = IX(s, i, j, k);
                double div = (s->u[n + X] - s->u[n]) / s->dx + (s->v[n + Y] - s->v[n]) / s->dy
                           + (s->w[n + Z] - s->w[n]) / s->dz;
                m = dmax(m, fabs(div));
            }
    return m;
}

int diag_init(Diagnostics *d, const Sim *s, const char *csv_path)
{
    const Params *p = &s->prm;
    d->V0 = diag_volume(s);
    d->pe_flat = 0.0;
    for (int k = 0; k < s->nz; k++)
        d->pe_flat += cell_pe(p, k, s->dz, flat_fraction(s, k, p->depth));
    d->pe_flat *= p->gravity * p->Lx * p->Ly * s->dz;
    d->volume_drift = 0.0;
    d->max_div = 0.0;
    d->n = 0;
    d->cap = 1024;
    d->t = malloc(d->cap * sizeof(double));
    d->eta = malloc(d->cap * sizeof(double));
    d->energy = malloc(d->cap * sizeof(double));
    d->csv = fopen(csv_path, "w");
    if (!d->t || !d->eta || !d->energy || !d->csv) {
        diag_free(d);
        return -1;
    }
    fprintf(d->csv, "t,eta_corner,eta_center,eta_far_corner,volume,kinetic,potential,max_div,poisson_iters\n");
    return 0;
}

void diag_sample(Diagnostics *d, const Sim *s)
{
    if (d->n == d->cap) {
        size_t cap = d->cap * 2;
        double *arrays[3] = { d->t, d->eta, d->energy };
        for (int a = 0; a < 3; a++) {
            double *grown = realloc(arrays[a], cap * sizeof(double));
            if (!grown)
                return;
            arrays[a] = grown;
        }
        d->t = arrays[0];
        d->eta = arrays[1];
        d->energy = arrays[2];
        d->cap = cap;
    }
    d->t[d->n] = s->t;
    d->eta[d->n] = diag_column_height(s, 0, 0);
    d->energy[d->n] = diag_kinetic(s) + diag_potential(s, d->pe_flat);
    d->n++;
    d->volume_drift = dmax(d->volume_drift, fabs(diag_volume(s) - d->V0) / d->V0);
    if (s->step > 0)
        d->max_div = dmax(d->max_div, diag_max_div(s));
}

void diag_write_row(Diagnostics *d, const Sim *s)
{
    fprintf(d->csv, "%.8g,%.10g,%.10g,%.10g,%.15g,%.8g,%.8g,%.4g,%d\n", s->t, diag_column_height(s, 0, 0),
            diag_column_height(s, s->nx / 2, s->ny / 2), diag_column_height(s, s->nx - 1, s->ny - 1), diag_volume(s),
            diag_kinetic(s), diag_potential(s, d->pe_flat), diag_max_div(s), s->pois_iters);
    fflush(d->csv);
}

void diag_summary(const Diagnostics *d, const Sim *s, Summary *out)
{
    const Params *p = &s->prm;
    const double h = p->depth;
    out->period_analytic = analytic_period(p);
    out->volume_drift = d->volume_drift;
    out->max_div = d->max_div;

    /* zero crossings of eta - h, half a period apart: least-squares slope of t_k vs k is T/2.
     * Samples in a small band around zero are skipped so round-off noise is not counted. */
    double emax = 0.0;
    for (size_t m = 0; m < d->n; m++)
        emax = dmax(emax, fabs(d->eta[m] - h));
    const double band = 1e-3 * emax;
    int n = 0;
    double sk = 0, st = 0, skk = 0, skt = 0;
    long prev = -1;
    for (size_t m = 0; m < d->n && emax > 1e-6 * s->dz; m++) {
        double b = d->eta[m] - h;
        if (fabs(b) <= band)
            continue;
        if (prev >= 0) {
            double a = d->eta[prev] - h;
            if ((a > 0.0) != (b > 0.0)) {
                double tc = d->t[prev] + (d->t[m] - d->t[prev]) * a / (a - b);
                sk += n;
                st += tc;
                skk += (double)n * n;
                skt += n * tc;
                n++;
            }
        }
        prev = (long)m;
    }
    out->crossings = n;
    if (n >= 3) {
        out->period_measured = 2.0 * (n * skt - sk * st) / (n * skk - sk * sk);
        out->period_error = (out->period_measured - out->period_analytic) / out->period_analytic;
    } else {
        out->period_measured = NAN;
        out->period_error = NAN;
    }

    /* wave energy (KE + PE above the flat state) averaged over the first and the last period */
    double T = out->period_analytic, t_last = d->n ? d->t[d->n - 1] : 0.0;
    out->energy_ratio = NAN;
    if (isfinite(T) && t_last >= 2.0 * T) {
        double e0 = 0.0, e1 = 0.0;
        long n0 = 0, n1 = 0;
        for (size_t m = 0; m < d->n; m++) {
            if (d->t[m] <= T) {
                e0 += d->energy[m];
                n0++;
            }
            if (d->t[m] >= t_last - T) {
                e1 += d->energy[m];
                n1++;
            }
        }
        if (n0 && n1 && e0 > 0.0)
            out->energy_ratio = (e1 / n1) / (e0 / n0);
    }
}

void diag_free(Diagnostics *d)
{
    free(d->t);
    free(d->eta);
    free(d->energy);
    if (d->csv)
        fclose(d->csv);
    d->t = d->eta = d->energy = NULL;
    d->csv = NULL;
}
