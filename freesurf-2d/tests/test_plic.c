/* PLIC geometry: volume <-> line-constant round trip and slab splitting. */
#include "util.h"
#include "vof.h"

#include <math.h>
#include <stdio.h>

static double urand(unsigned *state)
{
    *state = *state * 1664525u + 1013904223u;
    return (double)(*state >> 8) / 16777216.0;
}

static int check(const char *what, double got, double want, double tol)
{
    int ok = fabs(got - want) <= tol;
    printf("  %-42s %.15g (expected %.15g) %s\n", what, got, want, ok ? "ok" : "FAIL");
    return ok;
}

int main(void)
{
    int ok = 1;
    ok &= check("horizontal line z = 0.3", plic_area(0.0, 1.0, 0.3), 0.3, 1e-15);
    ok &= check("vertical line x = 0.25", plic_area(1.0, 0.0, 0.25), 0.25, 1e-15);
    ok &= check("diagonal corner triangle", plic_area(0.5, 0.5, 0.25), 0.125, 1e-15);
    ok &= check("liquid above (normal -z)", plic_area(0.0, -1.0, -0.3), 0.7, 1e-15);

    unsigned seed = 12345;
    double err_inv = 0.0, err_split_x = 0.0, err_split_z = 0.0;
    for (int n = 0; n < 200000; n++) {
        double th = 2.0 * FS_PI * urand(&seed);
        double nx = cos(th), nz = sin(th), m = fabs(nx) + fabs(nz);
        nx /= m;
        nz /= m;
        double c = n % 50 == 0 ? (double)(n % 100 == 0) : urand(&seed);
        double a = plic_alpha(c, nx, nz);
        err_inv = fmax(err_inv, fabs(plic_area(nx, nz, a) - c));

        /* the two slabs [0, s] and [s, 1] must hold all the liquid */
        double s = urand(&seed);
        double lx = s * plic_area(nx * s, nz, a) + (1.0 - s) * plic_area(nx * (1.0 - s), nz, a - nx * s);
        double lz = s * plic_area(nx, nz * s, a) + (1.0 - s) * plic_area(nx, nz * (1.0 - s), a - nz * s);
        err_split_x = fmax(err_split_x, fabs(lx - c));
        err_split_z = fmax(err_split_z, fabs(lz - c));
    }
    ok &= check("max |area(alpha(c)) - c|, 2e5 random", err_inv, 0.0, 1e-12);
    ok &= check("max x-slab split error", err_split_x, 0.0, 1e-12);
    ok &= check("max z-slab split error", err_split_z, 0.0, 1e-12);

    double pts[4];
    int seg = plic_segment(0.0, 1.0, 0.3, pts);
    ok &= check("segment of z = 0.3 exists", seg, 1, 0);
    ok &= check("segment length", hypot(pts[2] - pts[0], pts[3] - pts[1]), 1.0, 1e-12);

    printf(ok ? "PLIC: ok\n" : "PLIC: FAILED\n");
    return ok ? 0 : 1;
}
