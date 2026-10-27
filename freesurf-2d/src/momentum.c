#include "momentum.h"
#include "bc.h"
#include "vof.h"

#include <math.h>

/* ENO2 upwind derivative at f0 from samples f[-2..2] spaced h, upwinded by vel. */
static inline double eno2(double vel, double fm2, double fm1, double f0, double fp1, double fp2, double h)
{
    if (vel > 0.0) {
        double a = fm2 - 2.0 * fm1 + f0;
        double b = fm1 - 2.0 * f0 + fp1;
        double c = fabs(a) < fabs(b) ? a : b;
        return (f0 - fm1) / h + 0.5 * c / h;
    } else {
        double a = fm1 - 2.0 * f0 + fp1;
        double b = f0 - 2.0 * fp1 + fp2;
        double c = fabs(a) < fabs(b) ? a : b;
        return (fp1 - f0) / h - 0.5 * c / h;
    }
}

/* Gravity jump on a face: the interface values of the two cells, or the face height if neither is cut. */
static inline double face_phi(double a, double b, double phi_face)
{
    if (!isnan(a) && !isnan(b))
        return 0.5 * (a + b);
    if (!isnan(a))
        return a;
    if (!isnan(b))
        return b;
    return phi_face;
}

void momentum_predict(Sim *s, double dt)
{
    const int nx = s->nx, nz = s->nz;
    const double dx = s->dx, dz = s->dz, g = s->prm.gravity, zref = s->prm.depth;
    const double drho_g = (s->prm.rho_l - s->prm.rho_g) * g;
    const double *u = s->u, *w = s->w, *mu = s->mu, *C = s->C;
    double *tau = s->tau, *phi = s->phi;

    /* reduced gravity: phi = [rho] g (z_interface - depth) at the PLIC segment midpoint of each cut cell */
    vof_reconstruct(s, C);
    for (int j = -NG; j < nz + NG; j++)
        for (int i = -NG; i < nx + NG; i++)
            phi[IX(s, i, j)] = NAN;
    for (int j = 0; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t k = IX(s, i, j);
            double pts[4];
            if (C[k] > 1e-6 && C[k] < 1.0 - 1e-6 && plic_segment(s->pn_x[k], s->pn_z[k], s->palpha[k], pts))
                phi[k] = drho_g * ((j + 0.5 * (pts[1] + pts[3])) * dz - zref);
        }
    const double *rho = s->rho;

    /* shear stress mu (du/dz + dw/dx) at corners (i dx, j dz) */
    for (int j = 0; j <= nz; j++)
        for (int i = 0; i <= nx; i++) {
            double muc = 0.25 * (mu[IX(s, i - 1, j - 1)] + mu[IX(s, i, j - 1)] + mu[IX(s, i - 1, j)] + mu[IX(s, i, j)]);
            tau[IX(s, i, j)] = muc * ((u[IX(s, i, j)] - u[IX(s, i, j - 1)]) / dz
                                      + (w[IX(s, i, j)] - w[IX(s, i - 1, j)]) / dx);
        }

    /* x-momentum on interior x-faces */
    for (int j = 0; j < nz; j++) {
        s->us[IX(s, 0, j)] = 0.0;
        s->us[IX(s, nx, j)] = 0.0;
        for (int i = 1; i < nx; i++) {
            size_t k = IX(s, i, j);
            double uc = u[k];
            double wb = 0.25 * (w[IX(s, i - 1, j)] + w[IX(s, i, j)] + w[IX(s, i - 1, j + 1)] + w[IX(s, i, j + 1)]);
            double dudx = eno2(uc, u[IX(s, i - 2, j)], u[IX(s, i - 1, j)], uc, u[IX(s, i + 1, j)], u[IX(s, i + 2, j)], dx);
            double dudz = eno2(wb, u[IX(s, i, j - 2)], u[IX(s, i, j - 1)], uc, u[IX(s, i, j + 1)], u[IX(s, i, j + 2)], dz);

            double txx_r = 2.0 * mu[IX(s, i, j)] * (u[IX(s, i + 1, j)] - uc) / dx;
            double txx_l = 2.0 * mu[IX(s, i - 1, j)] * (uc - u[IX(s, i - 1, j)]) / dx;
            double rho_l = rho[IX(s, i - 1, j)], rho_r = rho[IX(s, i, j)];
            double rho_f = 0.5 * (rho_l + rho_r);
            double visc = ((txx_r - txx_l) / dx + (tau[IX(s, i, j + 1)] - tau[IX(s, i, j)]) / dz) / rho_f;
            double dC = C[IX(s, i, j)] - C[IX(s, i - 1, j)];
            double grav = dC == 0.0 ? 0.0
                        : face_phi(phi[IX(s, i - 1, j)], phi[IX(s, i, j)], drho_g * ((j + 0.5) * dz - zref)) * dC / (dx * rho_f);

            s->us[k] = uc + dt * (-(uc * dudx + wb * dudz) + visc + grav);
        }
    }

    /* z-momentum on interior z-faces */
    for (int i = 0; i < nx; i++) {
        s->ws[IX(s, i, 0)] = 0.0;
        s->ws[IX(s, i, nz)] = 0.0;
    }
    for (int j = 1; j < nz; j++)
        for (int i = 0; i < nx; i++) {
            size_t k = IX(s, i, j);
            double wc = w[k];
            double ub = 0.25 * (u[IX(s, i, j - 1)] + u[IX(s, i + 1, j - 1)] + u[IX(s, i, j)] + u[IX(s, i + 1, j)]);
            double dwdx = eno2(ub, w[IX(s, i - 2, j)], w[IX(s, i - 1, j)], wc, w[IX(s, i + 1, j)], w[IX(s, i + 2, j)], dx);
            double dwdz = eno2(wc, w[IX(s, i, j - 2)], w[IX(s, i, j - 1)], wc, w[IX(s, i, j + 1)], w[IX(s, i, j + 2)], dz);

            double tzz_t = 2.0 * mu[IX(s, i, j)] * (w[IX(s, i, j + 1)] - wc) / dz;
            double tzz_b = 2.0 * mu[IX(s, i, j - 1)] * (wc - w[IX(s, i, j - 1)]) / dz;
            double rho_b = rho[IX(s, i, j - 1)], rho_t = rho[IX(s, i, j)];
            double rho_f = 0.5 * (rho_b + rho_t);
            double visc = ((tau[IX(s, i + 1, j)] - tau[IX(s, i, j)]) / dx + (tzz_t - tzz_b) / dz) / rho_f;
            double dC = C[IX(s, i, j)] - C[IX(s, i, j - 1)];
            double grav = dC == 0.0 ? 0.0
                        : face_phi(phi[IX(s, i, j - 1)], phi[IX(s, i, j)], drho_g * (j * dz - zref)) * dC / (dz * rho_f);

            s->ws[k] = wc + dt * (-(ub * dwdx + wc * dwdz) + visc + grav);
        }

    bc_u(s, s->us);
    bc_w(s, s->ws);
}
