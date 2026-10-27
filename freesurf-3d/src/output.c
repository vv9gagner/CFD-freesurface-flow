#include "output.h"
#include "init.h"
#include "util.h"
#include "vtk_io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

_Static_assert(sizeof(fs_half) == 2, "fs_half must be IEEE binary16");

static const char *FIELD_NAME[OUT_NDISPLAY] = { "C", "p", "p_dyn", "speed", "vorticity" };
static const char *FIELD_LABEL[OUT_NDISPLAY] = {
    "Liquid fraction", "Pressure", "Dynamic pressure", "Speed", "Vorticity (normal to slice)",
};
static const char *FIELD_UNIT[OUT_NDISPLAY] = { "", "Pa", "Pa", "m/s", "1/s" };
static const char *STORED_NAME[OUT_NSTORED] = { "C", "p_dyn", "u", "v", "w" };

int output_init(Output *o, const Sim *s)
{
    const Params *p = &s->prm;
    const size_t N = (size_t)s->nx * s->ny * s->nz;
    o->n = 0;
    o->cap = 64;
    o->times = malloc(o->cap * sizeof(double));
    o->p_hydro = malloc((size_t)s->nz * sizeof(double));
    o->cu = malloc(N * sizeof(float));
    o->cv = malloc(N * sizeof(float));
    o->cw = malloc(N * sizeof(float));
    o->half = malloc(OUT_NSTORED * N * sizeof(fs_half));
    o->surface = malloc((size_t)s->nx * s->ny * sizeof(float));
    o->wall_start = wall_seconds();
    if (!o->times || !o->p_hydro || !o->cu || !o->cv || !o->cw || !o->half || !o->surface) {
        output_free(o);
        return -1;
    }
    for (int f = 0; f < OUT_NDISPLAY; f++) {
        o->fmin[f] = INFINITY;
        o->fmax[f] = -INFINITY;
    }
    o->eta_min = INFINITY;
    o->eta_max = -INFINITY;

    /* discrete hydrostatic profile of the flat surface, p = 0 in the top layer like the solver */
    o->p_hydro[s->nz - 1] = 0.0;
    for (int k = s->nz - 2; k >= 0; k--) {
        double c0 = flat_fraction(s, k, p->depth), c1 = flat_fraction(s, k + 1, p->depth);
        double rho_f = 0.5 * ((c0 + c1) * p->rho_l + (2.0 - c0 - c1) * p->rho_g);
        o->p_hydro[k] = o->p_hydro[k + 1] + rho_f * p->gravity * s->dz;
    }

    if (!host_is_little_endian()) {
        fprintf(stderr, "web output requires a little-endian host\n");
        output_free(o);
        return -1;
    }
    char dir[FS_PATH_MAX + 16];
    snprintf(dir, sizeof dir, "%s/frames", p->output_dir);
    int rc = mkdir_p(dir);
    snprintf(dir, sizeof dir, "%s/surface", p->output_dir);
    if (rc != 0 || mkdir_p(dir) != 0) {
        fprintf(stderr, "cannot create output directories in '%s'\n", p->output_dir);
        output_free(o);
        return -1;
    }
    return 0;
}

static fs_half to_half(double v)
{
    if (v > 65504.0)
        v = 65504.0;
    else if (v < -65504.0)
        v = -65504.0;
    return (fs_half)v;
}

/* central difference of a cell-centred float field, one-sided at the walls */
static double cdiff(const float *a, size_t m, ptrdiff_t o, int idx, int n, double h)
{
    if (idx == 0)
        return (a[m + o] - a[m]) / h;
    if (idx == n - 1)
        return (a[m] - a[m - o]) / h;
    return (a[m + o] - a[m - o]) / (2.0 * h);
}

static int write_atomic(const char *path, const void *data, size_t bytes)
{
    char tmp[FS_PATH_MAX + 80];
    FILE *f = atomic_open(path, tmp, sizeof tmp);
    if (!f)
        return -1;
    if (fwrite(data, 1, bytes, f) != bytes) {
        fclose(f);
        remove(tmp);
        return -1;
    }
    return atomic_commit(f, tmp, path);
}

static int write_frame_files(Output *o, const Sim *s, size_t index)
{
    const int nx = s->nx, ny = s->ny, nz = s->nz;
    const size_t N = (size_t)nx * ny * nz;
    const ptrdiff_t SX = 1, SY = nx, SZ = (ptrdiff_t)nx * ny;
    const double depth = s->prm.depth;
    fs_half *H = o->half;

    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                size_t m = (size_t)k * SZ + (size_t)j * SY + i;
                o->cu[m] = (float)sim_cell_u(s, i, j, k);
                o->cv[m] = (float)sim_cell_v(s, i, j, k);
                o->cw[m] = (float)sim_cell_w(s, i, j, k);
            }

    double lo[OUT_NDISPLAY], hi[OUT_NDISPLAY];
    for (int f = 0; f < OUT_NDISPLAY; f++) {
        lo[f] = INFINITY;
        hi[f] = -INFINITY;
    }
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                size_t m = (size_t)k * SZ + (size_t)j * SY + i;
                double c = s->C[IX(s, i, j, k)];
                double p = sim_cell_pressure(s, i, j, k);
                double pd = p - o->p_hydro[k];
                double speed = sqrt(sq(o->cu[m]) + sq(o->cv[m]) + sq(o->cw[m]));
                /* vorticity components from cell-centred velocities, as the viewer computes them */
                double wx = cdiff(o->cw, m, SY, j, ny, s->dy) - cdiff(o->cv, m, SZ, k, nz, s->dz);
                double wy = cdiff(o->cu, m, SZ, k, nz, s->dz) - cdiff(o->cw, m, SX, i, nx, s->dx);
                double wz = cdiff(o->cv, m, SX, i, nx, s->dx) - cdiff(o->cu, m, SY, j, ny, s->dy);
                double vort = dmax(fabs(wx), dmax(fabs(wy), fabs(wz)));
                double vals[OUT_NDISPLAY] = { c, p, pd, speed, vort };
                for (int f = 0; f < OUT_NDISPLAY; f++) {
                    lo[f] = dmin(lo[f], vals[f]);
                    hi[f] = dmax(hi[f], vals[f]);
                }
                H[m] = to_half(c);
                H[N + m] = to_half(pd);
                H[2 * N + m] = to_half(o->cu[m]);
                H[3 * N + m] = to_half(o->cv[m]);
                H[4 * N + m] = to_half(o->cw[m]);
            }
    for (int f = 0; f < OUT_NDISPLAY; f++) {
        o->fmin[f] = dmin(o->fmin[f], lo[f]);
        o->fmax[f] = dmax(o->fmax[f], hi[f]);
    }

    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            double eta = 0.0;
            for (int k = 0; k < nz; k++)
                eta += s->C[IX(s, i, j, k)];
            eta *= s->dz;
            o->surface[(size_t)j * nx + i] = (float)eta;
            o->eta_min = dmin(o->eta_min, eta - depth);
            o->eta_max = dmax(o->eta_max, eta - depth);
        }

    char path[FS_PATH_MAX + 64];
    snprintf(path, sizeof path, "%s/frames/frame_%06zu.bin", s->prm.output_dir, index);
    if (write_atomic(path, H, OUT_NSTORED * N * sizeof(fs_half)) != 0)
        return -1;
    snprintf(path, sizeof path, "%s/surface/surface_%06zu.bin", s->prm.output_dir, index);
    return write_atomic(path, o->surface, (size_t)nx * ny * sizeof(float));
}

int output_frame(Output *o, Sim *s)
{
    if (o->n == o->cap) {
        double *t = realloc(o->times, 2 * o->cap * sizeof(double));
        if (!t)
            return -1;
        o->times = t;
        o->cap *= 2;
    }
    size_t index = o->n;
    if (write_frame_files(o, s, index) != 0) {
        fprintf(stderr, "failed to write frame %zu\n", index);
        return -1;
    }
    o->times[o->n++] = s->t;
    if (s->prm.vtk_output && vtk_write_frame(s, (int)index, o->times, o->n) != 0) {
        fprintf(stderr, "failed to write VTK frame %zu\n", index);
        return -1;
    }
    return output_manifest(o, s, "running", NULL, NULL);
}

static void write_range(FILE *f, double lo, double hi)
{
    fputs("\"min\": ", f);
    json_write_number(f, lo);
    fputs(", \"max\": ", f);
    json_write_number(f, hi);
}

int output_manifest(const Output *o, const Sim *s, const char *status, const char *message, const Summary *sum)
{
    const Params *p = &s->prm;
    const size_t N = (size_t)s->nx * s->ny * s->nz;
    char path[FS_PATH_MAX + 32], tmp[FS_PATH_MAX + 40];
    snprintf(path, sizeof path, "%s/manifest.json", p->output_dir);
    FILE *f = atomic_open(path, tmp, sizeof tmp);
    if (!f)
        return -1;

    fputs("{\n  \"format\": \"freesurf3d-web-1\",\n  \"dim\": 3,\n  \"status\": ", f);
    json_write_string(f, status);
    fputs(",\n  \"message\": ", f);
    if (message)
        json_write_string(f, message);
    else
        fputs("null", f);
    fputs(",\n  \"case\": ", f);
    json_write_string(f, p->case_name);
    fprintf(f, ",\n  \"nx\": %d, \"ny\": %d, \"nz\": %d,\n  \"Lx\": ", s->nx, s->ny, s->nz);
    json_write_number(f, p->Lx);
    fputs(", \"Ly\": ", f);
    json_write_number(f, p->Ly);
    fputs(", \"Lz\": ", f);
    json_write_number(f, p->Lz);
    fputs(", \"depth\": ", f);
    json_write_number(f, p->depth);
    fputs(",\n  \"t\": ", f);
    json_write_number(f, s->t);
    fputs(", \"t_end\": ", f);
    json_write_number(f, p->t_end);
    fprintf(f, ", \"step\": %ld, \"progress\": ", s->step);
    json_write_number(f, dmin(1.0, s->t / p->t_end));
    fputs(", \"wall_time\": ", f);
    json_write_number(f, wall_seconds() - o->wall_start);
    fputs(",\n  \"analytic_period\": ", f);
    json_write_number(f, analytic_period(p));
    fputs(",\n  \"params\": ", f);
    params_write_json(p, f);

    fputs(",\n  \"storage\": {\"type\": \"float16\", \"layout\": \"field-major, x fastest, then y, then z\", \"fields\": [", f);
    for (int k = 0; k < OUT_NSTORED; k++) {
        if (k)
            fputs(", ", f);
        json_write_string(f, STORED_NAME[k]);
    }
    fputs("]},\n  \"p_hydro\": [", f);
    for (int k = 0; k < s->nz; k++) {
        if (k)
            fputs(", ", f);
        fprintf(f, "%.10g", o->p_hydro[k]);
    }
    fputs("],\n  \"fields\": [\n", f);
    for (int k = 0; k < OUT_NDISPLAY; k++) {
        double lo = o->fmin[k], hi = o->fmax[k];
        if (k == 0) { /* volume fraction: fixed range */
            lo = 0.0;
            hi = 1.0;
        } else if (k == 4) { /* vorticity components: symmetric about zero */
            lo = -hi;
        }
        fputs("    {\"name\": ", f);
        json_write_string(f, FIELD_NAME[k]);
        fputs(", \"label\": ", f);
        json_write_string(f, FIELD_LABEL[k]);
        fputs(", \"unit\": ", f);
        json_write_string(f, FIELD_UNIT[k]);
        fputs(", ", f);
        write_range(f, lo, hi);
        fputs(k + 1 < OUT_NDISPLAY ? "},\n" : "}\n", f);
    }
    fputs("  ],\n  \"surface_range\": {", f);
    write_range(f, o->eta_min, o->eta_max);
    fprintf(f, "},\n  \"frame_bytes\": %zu, \"surface_bytes\": %zu,\n  \"frame_times\": [",
            OUT_NSTORED * N * sizeof(fs_half), (size_t)s->nx * s->ny * sizeof(float));
    for (size_t m = 0; m < o->n; m++) {
        if (m)
            fputs(m % 10 ? ", " : ",\n    ", f);
        fprintf(f, "%.10g", o->times[m]);
    }
    fputs("],\n  \"summary\": ", f);
    if (sum) {
        fputs("{\"period_measured\": ", f);
        json_write_number(f, sum->period_measured);
        fputs(", \"period_analytic\": ", f);
        json_write_number(f, sum->period_analytic);
        fputs(", \"period_error\": ", f);
        json_write_number(f, sum->period_error);
        fprintf(f, ", \"crossings\": %d, \"volume_drift\": ", sum->crossings);
        json_write_number(f, sum->volume_drift);
        fputs(", \"energy_ratio\": ", f);
        json_write_number(f, sum->energy_ratio);
        fputs(", \"max_div\": ", f);
        json_write_number(f, sum->max_div);
        fputs("}", f);
    } else
        fputs("null", f);
    fputs("\n}\n", f);
    return atomic_commit(f, tmp, path);
}

void output_free(Output *o)
{
    free(o->times);
    free(o->p_hydro);
    free(o->cu);
    free(o->cv);
    free(o->cw);
    free(o->half);
    free(o->surface);
    o->times = o->p_hydro = NULL;
    o->cu = o->cv = o->cw = o->surface = NULL;
    o->half = NULL;
}
