#ifndef FS3_DIAGNOSTICS_H
#define FS3_DIAGNOSTICS_H

#include <stdio.h>

#include "sim.h"

typedef struct {
    double period_measured; /* from zero crossings of eta(corner) - h (NaN if too few) */
    double period_analytic; /* linear two-layer theory */
    double period_error;    /* relative */
    int crossings;
    double volume_drift;    /* max |V(t) - V0| / V0 */
    double energy_ratio;    /* mean wave energy over the last period / over the first period */
    double max_div;         /* max |div u| over the run */
} Summary;

typedef struct {
    double V0;
    double pe_flat;
    double volume_drift;
    double max_div;
    size_t n, cap;
    double *t, *eta, *energy; /* every step: time, corner surface height, wave energy */
    FILE *csv;
} Diagnostics;

/* Linear standing-wave period for two layers in a closed tank, k = pi sqrt((m/Lx)^2 + (n/Ly)^2):
 * omega^2 = g k (rho_l - rho_g) / (rho_l coth(k h) + rho_g coth(k (Lz - h))). */
double analytic_period(const Params *p);

double diag_column_height(const Sim *s, int i, int j);
double diag_volume(const Sim *s);
double diag_kinetic(const Sim *s);
double diag_potential(const Sim *s, double pe_flat);
double diag_max_div(const Sim *s);

int diag_init(Diagnostics *d, const Sim *s, const char *csv_path);
void diag_sample(Diagnostics *d, const Sim *s);
void diag_write_row(Diagnostics *d, const Sim *s);
void diag_summary(const Diagnostics *d, const Sim *s, Summary *out);
void diag_free(Diagnostics *d);

#endif
