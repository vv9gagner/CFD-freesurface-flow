/* 3D PLIC geometry: plane/cube volume against numerical integration, alpha <-> volume
 * round trip, slab splitting, and exact normals for planar interfaces. */
#include "config.h"
#include "sim.h"
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
    printf("  %-46s %.3e / %.6g %s\n", what, fabs(got - want), want, ok ? "ok" : "FAIL");
    return ok;
}

/* reference volume: composite Simpson in x over the exact 2D section area */
static double volume_reference(double nx, double ny, double nz, double alpha)
{
    const int M = 20000;
    double sum = 0.0;
    for (int m = 0; m <= M; m++) {
        double x = (double)m / M, wgt = (m == 0 || m == M) ? 1.0 : (m % 2 ? 4.0 : 2.0);
        sum += wgt * line_area(ny, nz, alpha - nx * x);
    }
    return sum / (3.0 * M);
}

static void random_normal(unsigned *seed, double n[3])
{
    double m = 0.0;
    for (int d = 0; d < 3; d++) {
        n[d] = 2.0 * urand(seed) - 1.0;
        m += fabs(n[d]);
    }
    for (int d = 0; d < 3; d++)
        n[d] /= m;
    /* now and then make it nearly or exactly 2D */
    unsigned pick = *seed % 7;
    if (pick == 0)
        n[*seed % 3] = 0.0;
    else if (pick == 1)
        n[*seed % 3] *= 1e-9;
}

/* C of cell (i,j,k) for the plane nv . (x - x0) = 0 in index coordinates, liquid on the negative side */
static double plane_fraction(const double nv[3], const double x0[3], int i, int j, int k)
{
    double alpha = nv[0] * (x0[0] - i) + nv[1] * (x0[1] - j) + nv[2] * (x0[2] - k);
    return plane_volume(nv[0], nv[1], nv[2], alpha);
}

static double normal_error(Sim *s, const double nv[3])
{
    const double x0[3] = { 4.5, 4.5, 4.5 };
    for (int k = 0; k < 9; k++)
        for (int j = 0; j < 9; j++)
            for (int i = 0; i < 9; i++)
                s->C[IX(s, i, j, k)] = plane_fraction(nv, x0, i, j, k);
    double got[3], m = fabs(nv[0]) + fabs(nv[1]) + fabs(nv[2]), err = 0.0;
    plic_normal(s, s->C, 4, 4, 4, got);
    for (int d = 0; d < 3; d++)
        err = fmax(err, fabs(got[d] - nv[d] / m));
    return err;
}

int main(void)
{
    int ok = 1;
    ok &= check("axis plane x = 0.3", plane_volume(1, 0, 0, 0.3), 0.3, 1e-15);
    ok &= check("axis plane, liquid above z = 0.25", plane_volume(0, 0, -1, -0.25), 0.75, 1e-15);
    ok &= check("corner tetrahedron x+y+z <= 0.5", plane_volume(1, 1, 1, 0.5), 0.125 / 6.0, 1e-15);
    ok &= check("half cube x+y+z <= 1.5", plane_volume(1, 1, 1, 1.5), 0.5, 1e-15);
    ok &= check("2D diagonal x+y <= 0.5", plane_volume(1, 1, 0, 0.5), 0.125, 1e-15);

    /* Simpson's reference is only accurate when no normal component is tiny (a tiny
     * component makes the section area nearly non-differentiable); degenerate normals are
     * covered by the round-trip and split checks below. */
    unsigned seed = 2024;
    double err_ref = 0.0;
    for (int n = 0; n < 300;) {
        double nv[3];
        random_normal(&seed, nv);
        if (fabs(nv[0]) < 0.05 || fabs(nv[1]) < 0.05 || fabs(nv[2]) < 0.05)
            continue;
        n++;
        double lo = fmin(nv[0], 0) + fmin(nv[1], 0) + fmin(nv[2], 0);
        double hi = fmax(nv[0], 0) + fmax(nv[1], 0) + fmax(nv[2], 0);
        double alpha = lo + urand(&seed) * (hi - lo);
        err_ref = fmax(err_ref, fabs(plane_volume(nv[0], nv[1], nv[2], alpha) - volume_reference(nv[0], nv[1], nv[2], alpha)));
    }
    ok &= check("max |V - numerical integral|, 300 random", err_ref, 0.0, 1e-9);

    double err_inv = 0.0, err_split = 0.0;
    for (int n = 0; n < 100000; n++) {
        double nv[3];
        random_normal(&seed, nv);
        double c = n % 50 == 0 ? (double)(n % 100 == 0) : urand(&seed);
        double a = plane_alpha(c, nv[0], nv[1], nv[2]);
        err_inv = fmax(err_inv, fabs(plane_volume(nv[0], nv[1], nv[2], a) - c));
        /* the two slabs [0, s] and [s, 1] along axis d hold all the liquid */
        int d = n % 3;
        double s = urand(&seed), n1[3] = { nv[0], nv[1], nv[2] }, n2[3] = { nv[0], nv[1], nv[2] };
        n1[d] *= s;
        n2[d] *= 1.0 - s;
        double split = s * plane_volume(n1[0], n1[1], n1[2], a)
                     + (1.0 - s) * plane_volume(n2[0], n2[1], n2[2], a - nv[d] * s);
        err_split = fmax(err_split, fabs(split - c));
    }
    ok &= check("max |V(alpha(c)) - c|, 1e5 random", err_inv, 0.0, 1e-12);
    ok &= check("max slab split error", err_split, 0.0, 1e-12);

    Params p;
    params_defaults(&p);
    p.nx = p.ny = p.nz = 8;
    Sim *s = sim_create(&p);
    if (!s)
        return 1;
    const double gentle[3] = { -0.2, 0.15, 1.0 };  /* nearly horizontal: centred columns are exact */
    const double sideways[3] = { 1.0, -0.25, 0.1 }; /* nearly vertical, facing +x */
    const double liquid_up[3] = { 0.1, 0.2, -1.0 }; /* liquid above the plane */
    ok &= check("normal of a gentle plane", normal_error(s, gentle), 0.0, 1e-10);
    ok &= check("normal of a near-vertical plane", normal_error(s, sideways), 0.0, 1e-10);
    ok &= check("normal with liquid above", normal_error(s, liquid_up), 0.0, 1e-10);
    const double steep[3] = { 0.6, 0.5, 0.7 };
    ok &= check("normal of a steep plane (Youngs, approx.)", normal_error(s, steep), 0.0, 0.05);
    sim_destroy(s);

    printf(ok ? "PLIC 3D: ok\n" : "PLIC 3D: FAILED\n");
    return ok ? 0 : 1;
}
