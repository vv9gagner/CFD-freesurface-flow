#define _POSIX_C_SOURCE 200809L
#include "config.h"
#include "util.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum { P_INT, P_DOUBLE, P_BOOL, P_STRING } PType;

typedef struct {
    const char *name;
    PType type;
    size_t offset;
    size_t size;
} PDesc;

#define FIELD(key, type, member) { key, type, offsetof(Params, member), sizeof(((Params *)0)->member) }

static const PDesc PARAMS[] = {
    FIELD("case", P_STRING, case_name),
    FIELD("depth", P_DOUBLE, depth),
    FIELD("amplitude", P_DOUBLE, amplitude),
    FIELD("mode_x", P_INT, mode_x),
    FIELD("mode_y", P_INT, mode_y),
    FIELD("Lx", P_DOUBLE, Lx),
    FIELD("Ly", P_DOUBLE, Ly),
    FIELD("Lz", P_DOUBLE, Lz),
    FIELD("nx", P_INT, nx),
    FIELD("ny", P_INT, ny),
    FIELD("nz", P_INT, nz),
    FIELD("rho_liquid", P_DOUBLE, rho_l),
    FIELD("rho_gas", P_DOUBLE, rho_g),
    FIELD("mu_liquid", P_DOUBLE, mu_l),
    FIELD("mu_gas", P_DOUBLE, mu_g),
    FIELD("gravity", P_DOUBLE, gravity),
    FIELD("wall_bc", P_STRING, wall_bc),
    FIELD("t_end", P_DOUBLE, t_end),
    FIELD("cfl", P_DOUBLE, cfl),
    FIELD("dt_max", P_DOUBLE, dt_max),
    FIELD("output_dt", P_DOUBLE, output_dt),
    FIELD("poisson_tol", P_DOUBLE, poisson_tol),
    FIELD("poisson_maxit", P_INT, poisson_maxit),
    FIELD("output_dir", P_STRING, output_dir),
    FIELD("vtk_output", P_BOOL, vtk_output),
    FIELD("print_every", P_INT, print_every),
    FIELD("threads", P_INT, threads),
};

#define NPARAMS (sizeof PARAMS / sizeof PARAMS[0])

void params_defaults(Params *p)
{
    memset(p, 0, sizeof *p);
    snprintf(p->case_name, sizeof p->case_name, "standing_wave");
    p->depth = 0.25;
    p->amplitude = 0.01;
    p->mode_x = 1;
    p->mode_y = 1;
    p->Lx = 1.0;
    p->Ly = 0.5;
    p->Lz = 0.5;
    p->nx = 64;
    p->ny = 32;
    p->nz = 32;
    p->rho_l = 1000.0;
    p->rho_g = 1.2;
    p->mu_l = 1.0e-3;
    p->mu_g = 1.8e-5;
    p->gravity = 9.81;
    snprintf(p->wall_bc, sizeof p->wall_bc, "slip");
    p->t_end = 10.0;
    p->cfl = 0.3;
    p->dt_max = 5.0e-3;
    p->output_dt = 0.05;
    p->poisson_tol = 1.0e-9;
    p->poisson_maxit = 5000;
    snprintf(p->output_dir, sizeof p->output_dir, "output/standing_wave_3d");
    p->vtk_output = 0;
    p->print_every = 100;
    p->threads = 0;
}

static const PDesc *find_param(const char *key)
{
    for (size_t n = 0; n < NPARAMS; n++)
        if (strcmp(PARAMS[n].name, key) == 0)
            return &PARAMS[n];
    return NULL;
}

int params_set(Params *p, const char *key, const char *value)
{
    const PDesc *d = find_param(key);
    if (!d) {
        fprintf(stderr, "unknown parameter '%s'\n", key);
        return -1;
    }
    char *base = (char *)p + d->offset;
    char *end = NULL;
    errno = 0;
    switch (d->type) {
    case P_INT: {
        long v = strtol(value, &end, 10);
        if (errno || end == value || *end != '\0' || v < -1000000000L || v > 1000000000L) {
            fprintf(stderr, "parameter '%s': '%s' is not an integer\n", key, value);
            return -1;
        }
        *(int *)base = (int)v;
        break;
    }
    case P_DOUBLE: {
        double v = strtod(value, &end);
        if (errno || end == value || *end != '\0' || !isfinite(v)) {
            fprintf(stderr, "parameter '%s': '%s' is not a number\n", key, value);
            return -1;
        }
        *(double *)base = v;
        break;
    }
    case P_BOOL:
        if (!strcmp(value, "1") || !strcmp(value, "true") || !strcmp(value, "yes") || !strcmp(value, "on"))
            *(int *)base = 1;
        else if (!strcmp(value, "0") || !strcmp(value, "false") || !strcmp(value, "no") || !strcmp(value, "off"))
            *(int *)base = 0;
        else {
            fprintf(stderr, "parameter '%s': '%s' is not a boolean (0/1)\n", key, value);
            return -1;
        }
        break;
    case P_STRING:
        if (strlen(value) >= d->size) {
            fprintf(stderr, "parameter '%s': value too long\n", key);
            return -1;
        }
        memcpy(base, value, strlen(value) + 1);
        break;
    }
    return 0;
}

static char *trim(char *s)
{
    while (isspace((unsigned char)*s))
        s++;
    char *e = s + strlen(s);
    while (e > s && isspace((unsigned char)e[-1]))
        *--e = '\0';
    return s;
}

int params_load_file(Params *p, const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "cannot open case file '%s': %s\n", path, strerror(errno));
        return -1;
    }
    char line[2048];
    int lineno = 0, rc = 0;
    while (fgets(line, sizeof line, f)) {
        lineno++;
        char *hash = strchr(line, '#');
        if (hash)
            *hash = '\0';
        char *s = trim(line);
        if (*s == '\0')
            continue;
        char *eq = strchr(s, '=');
        if (!eq) {
            fprintf(stderr, "%s:%d: expected 'key = value'\n", path, lineno);
            rc = -1;
            break;
        }
        *eq = '\0';
        if (params_set(p, trim(s), trim(eq + 1)) != 0) {
            fprintf(stderr, "%s:%d: invalid line\n", path, lineno);
            rc = -1;
            break;
        }
    }
    fclose(f);
    return rc;
}

int params_apply_arg(Params *p, const char *kv)
{
    char buf[2048];
    if (strlen(kv) >= sizeof buf) {
        fprintf(stderr, "argument too long\n");
        return -1;
    }
    strcpy(buf, kv);
    char *eq = strchr(buf, '=');
    if (!eq) {
        fprintf(stderr, "expected key=value, got '%s'\n", kv);
        return -1;
    }
    *eq = '\0';
    return params_set(p, trim(buf), trim(eq + 1));
}

#define FAIL(...)                           \
    do {                                    \
        snprintf(err, errlen, __VA_ARGS__); \
        return -1;                          \
    } while (0)

int params_validate(const Params *p, char *err, size_t errlen)
{
    if (strcmp(p->case_name, "standing_wave") != 0)
        FAIL("case must be 'standing_wave'");
    if (p->nx < 4 || p->ny < 4 || p->nz < 8)
        FAIL("need nx, ny >= 4 and nz >= 8");
    if ((long)p->nx * p->ny * p->nz > 8L * 1024 * 1024)
        FAIL("grid too large (nx*ny*nz > 8M)");
    if (!(p->Lx > 0) || !(p->Ly > 0) || !(p->Lz > 0))
        FAIL("Lx, Ly and Lz must be positive");
    if (p->amplitude < 0)
        FAIL("amplitude must be >= 0");
    if (!(p->depth - p->amplitude > 0) || !(p->depth + p->amplitude < p->Lz))
        FAIL("need 0 < depth - amplitude and depth + amplitude < Lz");
    if (p->mode_x < 0 || p->mode_y < 0 || p->mode_x + p->mode_y == 0)
        FAIL("modes must be >= 0 and not both zero");
    if (!(p->rho_l > 0) || !(p->rho_g > 0) || p->rho_g >= p->rho_l)
        FAIL("need 0 < rho_gas < rho_liquid");
    if (p->mu_l < 0 || p->mu_g < 0)
        FAIL("viscosities must be >= 0");
    if (p->gravity < 0)
        FAIL("gravity must be >= 0");
    if (strcmp(p->wall_bc, "slip") != 0 && strcmp(p->wall_bc, "noslip") != 0)
        FAIL("wall_bc must be 'slip' or 'noslip'");
    if (!(p->t_end > 0))
        FAIL("t_end must be positive");
    if (!(p->cfl > 0) || p->cfl > 0.5)
        FAIL("cfl must be in (0, 0.5] (VOF boundedness)");
    if (!(p->dt_max > 0) || !(p->output_dt > 0))
        FAIL("dt_max and output_dt must be positive");
    if (p->t_end / p->output_dt > 100000)
        FAIL("too many output frames (t_end/output_dt > 100000)");
    if (!(p->poisson_tol > 0) || p->poisson_maxit < 1)
        FAIL("poisson_tol must be > 0 and poisson_maxit >= 1");
    if (p->output_dir[0] == '\0')
        FAIL("output_dir must not be empty");
    if (p->print_every < 1 || p->threads < 0)
        FAIL("print_every must be >= 1 and threads >= 0");
    return 0;
}

void params_print(const Params *p, FILE *f)
{
    for (size_t n = 0; n < NPARAMS; n++) {
        const char *base = (const char *)p + PARAMS[n].offset;
        fprintf(f, "  %-14s = ", PARAMS[n].name);
        switch (PARAMS[n].type) {
        case P_INT:
        case P_BOOL: fprintf(f, "%d\n", *(const int *)base); break;
        case P_DOUBLE: fprintf(f, "%g\n", *(const double *)base); break;
        case P_STRING: fprintf(f, "%s\n", base); break;
        }
    }
}

void params_write_json(const Params *p, FILE *f)
{
    fputc('{', f);
    for (size_t n = 0; n < NPARAMS; n++) {
        const char *base = (const char *)p + PARAMS[n].offset;
        if (n)
            fputs(", ", f);
        json_write_string(f, PARAMS[n].name);
        fputs(": ", f);
        switch (PARAMS[n].type) {
        case P_INT: fprintf(f, "%d", *(const int *)base); break;
        case P_BOOL: fputs(*(const int *)base ? "true" : "false", f); break;
        case P_DOUBLE: json_write_number(f, *(const double *)base); break;
        case P_STRING: json_write_string(f, base); break;
        }
    }
    fputc('}', f);
}
