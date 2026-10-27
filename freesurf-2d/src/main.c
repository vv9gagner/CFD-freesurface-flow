#define _POSIX_C_SOURCE 200809L
#include "config.h"
#include "diagnostics.h"
#include "init.h"
#include "output.h"
#include "sim.h"
#include "solver.h"
#include "util.h"

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static volatile sig_atomic_t g_stop = 0;

static void on_signal(int sig)
{
    (void)sig;
    g_stop = 1;
}

static void usage(const char *prog)
{
    fprintf(stderr,
            "usage: %s <case.cfg> [key=value ...]\n"
            "  2D two-phase (VOF) Navier-Stokes solver for free-surface waves.\n"
            "  Parameters in the case file can be overridden on the command line,\n"
            "  e.g.  %s cases/standing_wave.cfg nx=256 nz=128 output_dir=output/fine\n",
            prog, prog);
}

static void print_summary(const Summary *sum)
{
    printf("\n=== summary ===\n");
    printf("  analytic period (linear, two-layer) : %.6f s\n", sum->period_analytic);
    if (isfinite(sum->period_measured))
        printf("  measured period (%2d zero crossings) : %.6f s   (error %+.3f %%)\n",
               sum->crossings, sum->period_measured, 100.0 * sum->period_error);
    else
        printf("  measured period                     : n/a (too few zero crossings; run longer)\n");
    if (isfinite(sum->energy_ratio))
        printf("  wave energy ratio last/first period : %.4f\n", sum->energy_ratio);
    printf("  max liquid volume drift             : %.3e (relative)\n", sum->volume_drift);
    printf("  max |div u|                         : %.3e 1/s\n", sum->max_div);
}

int main(int argc, char **argv)
{
    if (argc < 2 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        usage(argv[0]);
        return argc < 2 ? 1 : 0;
    }

    Params prm;
    params_defaults(&prm);
    if (params_load_file(&prm, argv[1]) != 0)
        return 1;
    for (int a = 2; a < argc; a++)
        if (params_apply_arg(&prm, argv[a]) != 0)
            return 1;
    char err[256];
    if (params_validate(&prm, err, sizeof err) != 0) {
        fprintf(stderr, "invalid parameters: %s\n", err);
        return 1;
    }
    printf("freesurf: 2D two-phase VOF Navier-Stokes\n");
    params_print(&prm, stdout);

    Sim *s = sim_create(&prm);
    if (!s) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }
    if (init_case(s) != 0) {
        sim_destroy(s);
        return 1;
    }

    Output out;
    if (output_init(&out, s) != 0) {
        sim_destroy(s);
        return 1;
    }
    Diagnostics diag;
    char csv_path[FS_PATH_MAX + 16];
    snprintf(csv_path, sizeof csv_path, "%s/probes.csv", prm.output_dir);
    if (diag_init(&diag, s, csv_path) != 0) {
        fprintf(stderr, "cannot open '%s'\n", csv_path);
        output_free(&out);
        sim_destroy(s);
        return 1;
    }

    struct sigaction sa;
    memset(&sa, 0, sizeof sa);
    sa.sa_handler = on_signal;
    sigaction(SIGINT, &sa, NULL);
    sigaction(SIGTERM, &sa, NULL);

    printf("\ngrid %dx%d, dx=%.4g m, dz=%.4g m, analytic period %.5f s\n",
           s->nx, s->nz, s->dx, s->dz, analytic_period(&prm));
    diag_sample(&diag, s);
    diag_write_row(&diag, s);
    const char *status = "done", *message = NULL;
    int rc = 0;
    if (output_frame(&out, s) != 0) {
        status = "failed";
        message = "could not write output";
        rc = 2;
    }

    long frames_written = 1, poisson_failures = 0;
    double t0 = wall_seconds();
    while (rc == 0 && s->t < prm.t_end - 1e-12) {
        if (g_stop) {
            status = "cancelled";
            message = "stopped by signal";
            rc = 3;
            break;
        }
        double next_out = dmin(frames_written * prm.output_dt, prm.t_end);
        double dt = solver_dt(s);
        if (s->t + dt > next_out - 1e-9 * prm.output_dt)
            dt = next_out - s->t; /* land exactly on the output time */

        int iters = solver_step(s, dt);
        if (iters < 0 && poisson_failures++ < 10)
            fprintf(stderr, "warning: pressure solve not converged at step %ld (residual %.2e)\n",
                    s->step, s->pois_res);
        diag_sample(&diag, s);

        double probe = diag_column_height(s, 0) - prm.depth;
        if (!isfinite(probe) || !isfinite(s->p[IX(s, 0, 0)])) {
            status = "failed";
            message = "solution diverged (non-finite values)";
            fprintf(stderr, "error: %s at step %ld, t=%g\n", message, s->step, s->t);
            rc = 2;
            break;
        }
        if (s->step % prm.print_every == 0)
            printf("step %7ld  t=%8.4f  dt=%.3e  poisson=%4d  div=%.2e  dV/V0=%+.2e  eta_L-h=%+.4e\n",
                   s->step, s->t, dt, s->pois_iters, diag_max_div(s),
                   (diag_volume(s) - diag.V0) / diag.V0, probe);

        if (s->t >= next_out - 1e-9 * prm.output_dt) {
            diag_write_row(&diag, s);
            if (output_frame(&out, s) != 0) {
                status = "failed";
                message = "could not write output";
                rc = 2;
                break;
            }
            frames_written++;
        }
    }

    Summary sum;
    diag_summary(&diag, s, &sum);
    output_manifest(&out, s, status, message, &sum);
    printf("\n%s after %ld steps, t=%.4f s, %ld frames, %.1f s wall time\n",
           status, s->step, s->t, frames_written, wall_seconds() - t0);
    print_summary(&sum);
    printf("\noutput: %s\n", prm.output_dir);

    diag_free(&diag);
    output_free(&out);
    sim_destroy(s);
    return rc;
}
