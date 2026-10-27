#ifndef FS_CONFIG_H
#define FS_CONFIG_H

#include <stddef.h>
#include <stdio.h>

#define FS_PATH_MAX 1024

typedef struct {
    /* scenario */
    char case_name[64];
    double depth;       /* mean liquid depth h [m] */
    double amplitude;   /* initial surface amplitude A [m] */
    int mode;           /* standing-wave mode n, k = n*pi/Lx */

    /* domain and grid */
    double Lx, Lz;      /* tank length and height [m] */
    int nx, nz;

    /* fluids */
    double rho_l, rho_g;  /* densities [kg/m^3] */
    double mu_l, mu_g;    /* dynamic viscosities [Pa s] */
    double gravity;       /* [m/s^2], acts in -z */
    char wall_bc[16];     /* "slip" or "noslip" */

    /* time stepping */
    double t_end, cfl, dt_max, output_dt;

    /* pressure solver */
    double poisson_tol;  /* relative residual (max-norm) */
    int poisson_maxit;

    /* output */
    char output_dir[FS_PATH_MAX];
    int vtk_output;
    int print_every;
} Params;

void params_defaults(Params *p);

/* All return 0 on success and print a message to stderr on failure. */
int params_load_file(Params *p, const char *path);
int params_set(Params *p, const char *key, const char *value);
int params_apply_arg(Params *p, const char *key_eq_value);

/* Returns 0 if the parameters are usable, otherwise writes a reason into err. */
int params_validate(const Params *p, char *err, size_t errlen);

void params_print(const Params *p, FILE *f);
void params_write_json(const Params *p, FILE *f);

#endif
