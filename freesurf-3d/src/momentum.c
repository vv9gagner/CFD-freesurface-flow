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

/* derivative of field a at index n along stride o (spacing h), upwinded by vel */
static inline double deriv(const double *a, size_t n, ptrdiff_t o, double vel, double h)
{
    return eno2(vel, a[n - 2 * o], a[n - o], a[n], a[n + o], a[n + 2 * o], h);
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
    const int nx = s->nx, ny = s->ny, nz = s->nz;
    const double dx = s->dx, dy = s->dy, dz = s->dz;
    const double zref = s->prm.depth, drho_g = (s->prm.rho_l - s->prm.rho_g) * s->prm.gravity;
    const ptrdiff_t X = 1, Y = s->Sx, Z = s->Sy;
    const double *u = s->u, *v = s->v, *w = s->w, *mu = s->mu, *rho = s->rho, *C = s->C;
    double *txy = s->tau_xy, *txz = s->tau_xz, *tyz = s->tau_yz, *phi = s->phi;

    /* reduced gravity: phi = [rho] g (z_interface - depth) in every cut cell */
    vof_reconstruct(s, C);
#pragma omp parallel for schedule(static)
    for (size_t n = 0; n < s->ntot; n++)
        phi[n] = NAN;
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                size_t n = IX(s, i, j, k);
                if (C[n] > 1e-6 && C[n] < 1.0 - 1e-6)
                    phi[n] = drho_g * ((k + plane_height(s->pn_x[n], s->pn_y[n], s->pn_z[n], s->palpha[n])) * dz - zref);
            }

    /* shear stresses on cell edges, viscosity averaged over the 4 surrounding cells */
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k <= nz; k++)
        for (int j = 0; j <= ny; j++)
            for (int i = 0; i <= nx; i++) {
                size_t n = IX(s, i, j, k);
                if (j < ny) /* tau_xz at (i dx, (j+.5) dy, k dz) */
                    txz[n] = 0.25 * (mu[n - X - Z] + mu[n - Z] + mu[n - X] + mu[n])
                           * ((u[n] - u[n - Z]) / dz + (w[n] - w[n - X]) / dx);
                if (k < nz) /* tau_xy at (i dx, j dy, (k+.5) dz) */
                    txy[n] = 0.25 * (mu[n - X - Y] + mu[n - Y] + mu[n - X] + mu[n])
                           * ((u[n] - u[n - Y]) / dy + (v[n] - v[n - X]) / dx);
                if (i < nx) /* tau_yz at ((i+.5) dx, j dy, k dz) */
                    tyz[n] = 0.25 * (mu[n - Y - Z] + mu[n - Z] + mu[n - Y] + mu[n])
                           * ((v[n] - v[n - Z]) / dz + (w[n] - w[n - Y]) / dy);
            }

    /* x-momentum on interior x-faces */
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++) {
            s->us[IX(s, 0, j, k)] = 0.0;
            s->us[IX(s, nx, j, k)] = 0.0;
            for (int i = 1; i < nx; i++) {
                size_t n = IX(s, i, j, k);
                double uc = u[n];
                double vb = 0.25 * (v[n - X] + v[n] + v[n - X + Y] + v[n + Y]);
                double wb = 0.25 * (w[n - X] + w[n] + w[n - X + Z] + w[n + Z]);
                double adv = uc * deriv(u, n, X, uc, dx) + vb * deriv(u, n, Y, vb, dy) + wb * deriv(u, n, Z, wb, dz);
                double txx_r = 2.0 * mu[n] * (u[n + X] - uc) / dx;
                double txx_l = 2.0 * mu[n - X] * (uc - u[n - X]) / dx;
                double rho_f = 0.5 * (rho[n - X] + rho[n]);
                double visc = ((txx_r - txx_l) / dx + (txy[n + Y] - txy[n]) / dy + (txz[n + Z] - txz[n]) / dz) / rho_f;
                double dC = C[n] - C[n - X];
                double grav = dC == 0.0 ? 0.0
                            : face_phi(phi[n - X], phi[n], drho_g * ((k + 0.5) * dz - zref)) * dC / (dx * rho_f);
                s->us[n] = uc + dt * (-adv + visc + grav);
            }
        }

    /* y-momentum on interior y-faces */
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < nz; k++)
        for (int j = 0; j <= ny; j++)
            for (int i = 0; i < nx; i++) {
                size_t n = IX(s, i, j, k);
                if (j == 0 || j == ny) {
                    s->vs[n] = 0.0;
                    continue;
                }
                double vc = v[n];
                double ub = 0.25 * (u[n - Y] + u[n + X - Y] + u[n] + u[n + X]);
                double wb = 0.25 * (w[n - Y] + w[n] + w[n - Y + Z] + w[n + Z]);
                double adv = ub * deriv(v, n, X, ub, dx) + vc * deriv(v, n, Y, vc, dy) + wb * deriv(v, n, Z, wb, dz);
                double tyy_t = 2.0 * mu[n] * (v[n + Y] - vc) / dy;
                double tyy_b = 2.0 * mu[n - Y] * (vc - v[n - Y]) / dy;
                double rho_f = 0.5 * (rho[n - Y] + rho[n]);
                double visc = ((txy[n + X] - txy[n]) / dx + (tyy_t - tyy_b) / dy + (tyz[n + Z] - tyz[n]) / dz) / rho_f;
                double dC = C[n] - C[n - Y];
                double grav = dC == 0.0 ? 0.0
                            : face_phi(phi[n - Y], phi[n], drho_g * ((k + 0.5) * dz - zref)) * dC / (dy * rho_f);
                s->vs[n] = vc + dt * (-adv + visc + grav);
            }

    /* z-momentum on interior z-faces */
#pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k <= nz; k++)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                size_t n = IX(s, i, j, k);
                if (k == 0 || k == nz) {
                    s->ws[n] = 0.0;
                    continue;
                }
                double wc = w[n];
                double ub = 0.25 * (u[n - Z] + u[n + X - Z] + u[n] + u[n + X]);
                double vb = 0.25 * (v[n - Z] + v[n + Y - Z] + v[n] + v[n + Y]);
                double adv = ub * deriv(w, n, X, ub, dx) + vb * deriv(w, n, Y, vb, dy) + wc * deriv(w, n, Z, wc, dz);
                double tzz_t = 2.0 * mu[n] * (w[n + Z] - wc) / dz;
                double tzz_b = 2.0 * mu[n - Z] * (wc - w[n - Z]) / dz;
                double rho_f = 0.5 * (rho[n - Z] + rho[n]);
                double visc = ((txz[n + X] - txz[n]) / dx + (tyz[n + Y] - tyz[n]) / dy + (tzz_t - tzz_b) / dz) / rho_f;
                double dC = C[n] - C[n - Z];
                double grav = dC == 0.0 ? 0.0
                            : face_phi(phi[n - Z], phi[n], drho_g * (k * dz - zref)) * dC / (dz * rho_f);
                s->ws[n] = wc + dt * (-adv + visc + grav);
            }

    bc_u(s, s->us);
    bc_v(s, s->vs);
    bc_w(s, s->ws);
}
