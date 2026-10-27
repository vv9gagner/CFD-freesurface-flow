#include "bc.h"

/*
 * Fill the ghost layers of `a` along direction d (0 = x, 1 = y, 2 = z), over the full
 * extent (ghosts included) of the two other directions, so successive passes also fill
 * edges and corners consistently.
 *   normal = 1: `a` is the velocity component normal to these walls; it lives on the
 *               wall points 0 and n, which are set to zero, with odd reflection beyond.
 *   normal = 0: `a` is cell-centred along d; ghosts mirror the interior times `sign`.
 */
static void ghosts(const Sim *s, double *a, int d, int normal, double sign)
{
    const int n[3] = { s->nx, s->ny, s->nz };
    const ptrdiff_t stride[3] = { 1, s->Sx, s->Sy };
    const int d1 = (d + 1) % 3, d2 = (d + 2) % 3;
    const int nd = n[d];
    const ptrdiff_t st = stride[d];

    for (int b = -NG; b <= n[d2] + NG; b++)
        for (int a1 = -NG; a1 <= n[d1] + NG; a1++) {
            int c[3];
            c[d] = 0;
            c[d1] = a1;
            c[d2] = b;
            double *p = a + IX(s, c[0], c[1], c[2]);
            if (normal) {
                p[0] = 0.0;
                p[nd * st] = 0.0;
                for (int m = 1; m <= NG; m++) {
                    p[-m * st] = -p[m * st];
                    p[(nd + m) * st] = -p[(nd - m) * st];
                }
            } else {
                for (int m = 1; m <= NG; m++) {
                    p[-m * st] = sign * p[(m - 1) * st];
                    p[(nd - 1 + m) * st] = sign * p[(nd - m) * st];
                }
            }
        }
}

void bc_u(const Sim *s, double *u)
{
    const double tang = s->noslip ? -1.0 : 1.0;
    ghosts(s, u, 0, 1, 0.0);
    ghosts(s, u, 1, 0, tang);
    ghosts(s, u, 2, 0, tang);
}

void bc_v(const Sim *s, double *v)
{
    const double tang = s->noslip ? -1.0 : 1.0;
    ghosts(s, v, 1, 1, 0.0);
    ghosts(s, v, 2, 0, tang);
    ghosts(s, v, 0, 0, tang);
}

void bc_w(const Sim *s, double *w)
{
    const double tang = s->noslip ? -1.0 : 1.0;
    ghosts(s, w, 2, 1, 0.0);
    ghosts(s, w, 0, 0, tang);
    ghosts(s, w, 1, 0, tang);
}

void bc_scalar(const Sim *s, double *c)
{
    ghosts(s, c, 0, 0, 1.0);
    ghosts(s, c, 1, 0, 1.0);
    ghosts(s, c, 2, 0, 1.0);
}
