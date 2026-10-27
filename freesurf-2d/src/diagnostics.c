#include "diagnostics.h"
#include "init.h"
#include "util.h"

#include <math.h>
#include <stdlib.h>

double analytic_period(const Params *p)
{
    double k = p->mode * FS_PI / p->Lx;
    double num = p->gravity * k * (p->rho_l - p->rho_g);
    double den = p->rho_l / tanh(k * p->depth) + p->rho_g / tanh(k * (p->Lz - p->depth));
    if (!(num > 0.0))
        return NAN;
    return 2.0 * FS_PI / sqrt(num / den);
}

double diag_column_height(const Sim *s, int i)
{
    double h = 0.0;
    for (int j = 0; j < s->nz; j++)
        h += s->C[IX(s, i, j)];
    return h * s->dz;
}

double diag_volume(const Sim *s)
{
    double v = 0.0;
    for (int j = 0; j < s->nz; j++)
        for (int i = 0; i < s->nx; i++)
            v += s->C[IX(s, i, j)];
    return v * s->dx * s->dz;
}

double diag_kinetic(const Sim *s)
{
    double e = 0.0;
    for (int j = 0; j < s->nz; j++)
        for (int i = 1; i < s->nx; i++) {
            double rho_f = 0.5 * (s->rho[IX(s, i - 1, j)] + s->rho[IX(s, i, j)]);
            e += 0.5 * rho_f * sq(s->u[IX(s, i, j)]);
        }
    for (int j = 1; j < s->nz; j++)
        for (int i = 0; i < s->nx; i++) {
            double rho_f = 0.5 * (s->rho[IX(s, i, j - 1)] + s->rho[IX(s, i, j)]);
            e += 0.5 * rho_f * sq(s->w[IX(s, i, j)]);
        }
    return e * s->dx * s->dz;
}

/* Potential energy of one cell per g dx dz, with the liquid fraction C sitting at the
 * bottom of the cell and the gas on top (exact for a horizontal interface). */
static double cell_pe(const Params *p, int j, double dz, double c)
{
    double zb = j * dz;
    return p->rho_l * c * (zb + 0.5 * c * dz) + p->rho_g * (1.0 - c) * (zb + 0.5 * (1.0 + c) * dz);
}

/* Gravitational potential energy (per unit width) relative to the flat state. */
double diag_potential(const Sim *s, double pe_flat)
{
    double e = 0.0;
    for (int j = 0; j < s->nz; j++)
        for (int i = 0; i < s->nx; i++)
            e += cell_pe(&s->prm, j, s->dz, s->C[IX(s, i, j)]);
    return e * s->prm.gravity * s->dx * s->dz - pe_flat;
}

double diag_max_div(const Sim *s)
{
    double m = 0.0;
    for (int j = 0; j < s->nz; j++)
        for (int i = 0; i < s->nx; i++) {
            double div = (s->u[IX(s, i + 1, j)] - s->u[IX(s, i, j)]) / s->dx
                       + (s->w[IX(s, i, j + 1)] - s->w[IX(s, i, j)]) / s->dz;
            m = dmax(m, fabs(div));
        }
    return m;
}

int diag_init(Diagnostics *d, const Sim *s, const char *csv_path)
{
    const Params *p = &s->prm;
    d->V0 = diag_volume(s);
    d->pe_flat = 0.0;
    for (int j = 0; j < s->nz; j++)
        d->pe_flat += cell_pe(p, j, s->dz, flat_fraction(s, j, p->depth));
    d->pe_flat *= p->gravity * p->Lx * s->dz;
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
    fprintf(d->csv, "t,eta_left,eta_center,eta_right,volume,kinetic,potential,max_div,poisson_iters\n");
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
    d->eta[d->n] = diag_column_height(s, 0);
    d->energy[d->n] = diag_kinetic(s) + diag_potential(s, d->pe_flat);
    d->n++;
    d->volume_drift = dmax(d->volume_drift, fabs(diag_volume(s) - d->V0) / d->V0);
    if (s->step > 0)
        d->max_div = dmax(d->max_div, diag_max_div(s));
}

void diag_write_row(Diagnostics *d, const Sim *s)
{
    fprintf(d->csv, "%.8g,%.10g,%.10g,%.10g,%.15g,%.8g,%.8g,%.4g,%d\n",
            s->t, diag_column_height(s, 0), diag_column_height(s, s->nx / 2), diag_column_height(s, s->nx - 1),
            diag_volume(s), diag_kinetic(s), diag_potential(s, d->pe_flat), diag_max_div(s), s->pois_iters);
    fflush(d->csv);
}

void diag_summary(const Diagnostics *d, const Sim *s, Summary *out)
{
    const Params *p = &s->prm;
    const double h = p->depth;
    out->period_analytic = analytic_period(p);
    out->volume_drift = d->volume_drift;
    out->max_div = d->max_div;

    /* zero crossings of eta - h with linear interpolation; they are half a period apart,
     * so a least-squares line through (k, t_k) has slope T/2. Samples within a small
     * band around zero are skipped (hysteresis) so round-off noise is not counted. */
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
        double slope = (n * skt - sk * st) / (n * skk - sk * sk);
        out->period_measured = 2.0 * slope;
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
