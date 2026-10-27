#ifndef FS3_CONFIG_H
#define FS3_CONFIG_H

#include <stddef.h>
#include <stdio.h>

#define FS_PATH_MAX 1024

typedef struct {
    /* scenario */
    char case_name[64];
    double depth;       /* mean liquid depth h [m] */
    double amplitude;   /* initial surface amplitude A [m] */
    int mode_x, mode_y; /* standing-wave mode (m, n): k = pi sqrt((m/Lx)^2 + (n/Ly)^2) */

    /* domain and grid; z is vertical */
    double Lx, Ly, Lz;
    int nx, ny, nz;

    /* fluids */
    double rho_l, rho_g;
    double mu_l, mu_g;
    double gravity;   /* acts in -z */
    char wall_bc[16]; /* "slip" or "noslip" */

    /* time stepping */
    double t_end, cfl, dt_max, output_dt;

    /* pressure solver */
    double poisson_tol;
    int poisson_maxit;

    /* output and execution */
    char output_dir[FS_PATH_MAX];
    int vtk_output;
    int print_every;
    int threads; /* OpenMP threads, 0 = runtime default */
} Params;

void params_defaults(Params *p);
int params_load_file(Params *p, const char *path);
int params_set(Params *p, const char *key, const char *value);
int params_apply_arg(Params *p, const char *key_eq_value);
int params_validate(const Params *p, char *err, size_t errlen);
void params_print(const Params *p, FILE *f);
void params_write_json(const Params *p, FILE *f);

#endif
