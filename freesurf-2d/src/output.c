#include "output.h"
#include "init.h"
#include "util.h"
#include "vtk_io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static const char *FIELD_NAME[OUT_NFIELDS] = { "C", "p", "p_dyn", "speed", "vorticity" };
static const char *FIELD_LABEL[OUT_NFIELDS] = {
    "Liquid fraction", "Pressure", "Dynamic pressure", "Speed", "Vorticity",
};
static const char *FIELD_UNIT[OUT_NFIELDS] = { "", "Pa", "Pa", "m/s", "1/s" };

int output_init(Output *o, const Sim *s)
{
    const Params *p = &s->prm;
    o->n = 0;
    o->cap = 64;
    o->times = malloc(o->cap * sizeof(double));
    o->p_hydro = malloc((size_t)s->nz * sizeof(double));
    o->buf = malloc((size_t)OUT_NFIELDS * s->nx * s->nz * sizeof(float));
    o->wall_start = wall_seconds();
    if (!o->times || !o->p_hydro || !o->buf) {
        output_free(o);
        return -1;
    }
    for (int f = 0; f < OUT_NFIELDS; f++) {
        o->fmin[f] = INFINITY;
        o->fmax[f] = -INFINITY;
    }

    /* discrete hydrostatic profile of the flat surface, p = 0 in the top row like the solver */
    o->p_hydro[s->nz - 1] = 0.0;
    for (int j = s->nz - 2; j >= 0; j--) {
        double c0 = flat_fraction(s, j, p->depth), c1 = flat_fraction(s, j + 1, p->depth);
        double rho_f = 0.5 * ((c0 + c1) * p->rho_l + (2.0 - c0 - c1) * p->rho_g);
        o->p_hydro[j] = o->p_hydro[j + 1] + rho_f * p->gravity * s->dz;
    }

    if (!host_is_little_endian()) {
        fprintf(stderr, "web output requires a little-endian host\n");
        output_free(o);
        return -1;
    }
    char dir[FS_PATH_MAX + 16];
    snprintf(dir, sizeof dir, "%s/frames", p->output_dir);
    if (mkdir_p(dir) != 0) {
        fprintf(stderr, "cannot create output directory '%s'\n", dir);
        output_free(o);
        return -1;
    }
    return 0;
}

static int write_frame_bin(Output *o, const Sim *s, size_t index)
{
    const int nx = s->nx, nz = s->nz;
    const size_t n = (size_t)nx * nz;
    float *F[OUT_NFIELDS];
    for (int f = 0; f < OUT_NFIELDS; f++)
        F[f] = o->buf + (size_t)f * n;

    for (int j = 0; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t m = (size_t)j * nx + i, k = IX(s, i, j);
            double v[OUT_NFIELDS];
            v[0] = s->C[k];
            v[1] = sim_cell_pressure(s, i, j);
            v[2] = v[1] - o->p_hydro[j];
            v[3] = hypot(sim_cell_u(s, i, j), sim_cell_w(s, i, j));
            v[4] = sim_cell_vorticity(s, i, j);
            for (int f = 0; f < OUT_NFIELDS; f++) {
                F[f][m] = (float)v[f];
                o->fmin[f] = dmin(o->fmin[f], v[f]);
                o->fmax[f] = dmax(o->fmax[f], v[f]);
            }
        }

    char path[FS_PATH_MAX + 64], tmp[FS_PATH_MAX + 72];
    snprintf(path, sizeof path, "%s/frames/frame_%06zu.bin", s->prm.output_dir, index);
    FILE *f = atomic_open(path, tmp, sizeof tmp);
    if (!f)
        return -1;
    if (fwrite(o->buf, sizeof(float), OUT_NFIELDS * n, f) != OUT_NFIELDS * n) {
        fclose(f);
        remove(tmp);
        return -1;
    }
    return atomic_commit(f, tmp, path);
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
    if (write_frame_bin(o, s, index) != 0) {
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

static void write_range(FILE *f, const Output *o, int field)
{
    double lo = o->fmin[field], hi = o->fmax[field];
    if (field == 0) { /* volume fraction: fixed range */
        lo = 0.0;
        hi = 1.0;
    } else if (field == 4) { /* vorticity: symmetric about zero */
        hi = dmax(fabs(lo), fabs(hi));
        lo = -hi;
    }
    fputs("\"min\": ", f);
    json_write_number(f, lo);
    fputs(", \"max\": ", f);
    json_write_number(f, hi);
}

int output_manifest(const Output *o, const Sim *s, const char *status, const char *message, const Summary *sum)
{
    const Params *p = &s->prm;
    char path[FS_PATH_MAX + 32], tmp[FS_PATH_MAX + 40];
    snprintf(path, sizeof path, "%s/manifest.json", p->output_dir);
    FILE *f = atomic_open(path, tmp, sizeof tmp);
    if (!f)
        return -1;

    fputs("{\n  \"format\": \"freesurf-web-1\",\n  \"dim\": 2,\n  \"status\": ", f);
    json_write_string(f, status);
    fputs(",\n  \"message\": ", f);
    if (message)
        json_write_string(f, message);
    else
        fputs("null", f);
    fputs(",\n  \"case\": ", f);
    json_write_string(f, p->case_name);
    fprintf(f, ",\n  \"nx\": %d, \"nz\": %d,\n  \"Lx\": ", s->nx, s->nz);
    json_write_number(f, p->Lx);
    fputs(", \"Lz\": ", f);
    json_write_number(f, p->Lz);
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

    fputs(",\n  \"fields\": [\n", f);
    for (int k = 0; k < OUT_NFIELDS; k++) {
        fputs("    {\"name\": ", f);
        json_write_string(f, FIELD_NAME[k]);
        fputs(", \"label\": ", f);
        json_write_string(f, FIELD_LABEL[k]);
        fputs(", \"unit\": ", f);
        json_write_string(f, FIELD_UNIT[k]);
        fputs(", ", f);
        write_range(f, o, k);
        fputs(k + 1 < OUT_NFIELDS ? "},\n" : "}\n", f);
    }
    fprintf(f, "  ],\n  \"frame_bytes\": %zu,\n  \"frame_times\": [",
            (size_t)OUT_NFIELDS * (size_t)s->nx * (size_t)s->nz * sizeof(float));
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
    free(o->buf);
    o->times = o->p_hydro = NULL;
    o->buf = NULL;
}
