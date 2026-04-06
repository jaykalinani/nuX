#ifndef NUX_M0_KERNEL_HXX
#define NUX_M0_KERNEL_HXX

#include <cmath>
#include <memory>

#include "cctk.h"
#include <loop_device.hxx>
#include "nuX_metric.hxx"
#include "nuX_sph_grid.hxx"
#include "nuX_valencia.hxx"

namespace nuX_M0 {

using nuX_Utils::sph_grid::SphericalGrid;

extern std::unique_ptr<SphericalGrid> M0Grid;

struct OwnedBox {
  CCTK_REAL xmin[3];
  CCTK_REAL xmax[3];
  int upper_inclusive[3];
};

CCTK_HOST CCTK_DEVICE inline int index(int const nrad, int const irad,
                                       int const iray) {
  return irad + nrad * iray;
}

CCTK_HOST CCTK_DEVICE inline int index(SphericalGrid const &grid, int const irad,
                                       int const itheta, int const iphi) {
  return irad + grid.nrad * (itheta + grid.ntheta * iphi);
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL signum(CCTK_REAL const x) {
  return x < 0.0 ? -1.0 : 1.0;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL dot3(CCTK_REAL const a[3],
                                            CCTK_REAL const b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL
dot4_metric(CCTK_REAL const g[16], CCTK_REAL const a[4], CCTK_REAL const b[4]) {
  CCTK_REAL sum = 0.0;
  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu)
      sum += g[4 * mu + nu] * a[mu] * b[nu];
  return sum;
}

CCTK_HOST CCTK_DEVICE inline void matvec3(CCTK_REAL const A[9],
                                          CCTK_REAL const x[3],
                                          CCTK_REAL y[3]) {
  for (int i = 0; i < 3; ++i) {
    y[i] = 0.0;
    for (int j = 0; j < 3; ++j)
      y[i] += A[3 * i + j] * x[j];
  }
}

CCTK_HOST CCTK_DEVICE inline bool inside_owned_box(OwnedBox const &box,
                                                   CCTK_REAL const x,
                                                   CCTK_REAL const y,
                                                   CCTK_REAL const z) {
  CCTK_REAL const epsx = 1.0e-12 * (box.xmax[0] - box.xmin[0] + 1.0);
  CCTK_REAL const epsy = 1.0e-12 * (box.xmax[1] - box.xmin[1] + 1.0);
  CCTK_REAL const epsz = 1.0e-12 * (box.xmax[2] - box.xmin[2] + 1.0);

  bool const in_x =
      x >= box.xmin[0] - epsx &&
      (box.upper_inclusive[0] ? x <= box.xmax[0] + epsx
                              : x < box.xmax[0] - epsx);
  bool const in_y =
      y >= box.xmin[1] - epsy &&
      (box.upper_inclusive[1] ? y <= box.xmax[1] + epsy
                              : y < box.xmax[1] - epsy);
  bool const in_z =
      z >= box.xmin[2] - epsz &&
      (box.upper_inclusive[2] ? z <= box.xmax[2] + epsz
                              : z < box.xmax[2] - epsz);
  return in_x && in_y && in_z;
}

template <typename T>
CCTK_HOST CCTK_DEVICE inline T sample_cart_linear(
    Loop::GF3D2<const T> const &gf, Loop::GridDescBase const &grid,
    CCTK_REAL const x, CCTK_REAL const y, CCTK_REAL const z,
    bool const cell_centered) {
  const int nx = grid.lsh[0] - (cell_centered ? 1 : 0);
  const int ny = grid.lsh[1] - (cell_centered ? 1 : 0);
  const int nz = grid.lsh[2] - (cell_centered ? 1 : 0);

  int i0, j0, k0;
  CCTK_REAL tx, ty, tz;
  {
    const CCTK_REAL shift = cell_centered ? 0.0 : 0.5;
    const CCTK_REAL q = (x - grid.x0[0]) / grid.dx[0] - grid.lbnd[0] + shift;
    if (q <= 0.0) {
      i0 = 0;
      tx = 0.0;
    } else if (q >= nx - 1) {
      i0 = nx - 2;
      tx = 1.0;
    } else {
      i0 = static_cast<int>(::floor(q));
      tx = q - i0;
    }
  }
  {
    const CCTK_REAL shift = cell_centered ? 0.0 : 0.5;
    const CCTK_REAL q = (y - grid.x0[1]) / grid.dx[1] - grid.lbnd[1] + shift;
    if (q <= 0.0) {
      j0 = 0;
      ty = 0.0;
    } else if (q >= ny - 1) {
      j0 = ny - 2;
      ty = 1.0;
    } else {
      j0 = static_cast<int>(::floor(q));
      ty = q - j0;
    }
  }
  {
    const CCTK_REAL shift = cell_centered ? 0.0 : 0.5;
    const CCTK_REAL q = (z - grid.x0[2]) / grid.dx[2] - grid.lbnd[2] + shift;
    if (q <= 0.0) {
      k0 = 0;
      tz = 0.0;
    } else if (q >= nz - 1) {
      k0 = nz - 2;
      tz = 1.0;
    } else {
      k0 = static_cast<int>(::floor(q));
      tz = q - k0;
    }
  }

  const T c000 = gf(i0, j0, k0);
  const T c100 = gf(i0 + 1, j0, k0);
  const T c010 = gf(i0, j0 + 1, k0);
  const T c110 = gf(i0 + 1, j0 + 1, k0);
  const T c001 = gf(i0, j0, k0 + 1);
  const T c101 = gf(i0 + 1, j0, k0 + 1);
  const T c011 = gf(i0, j0 + 1, k0 + 1);
  const T c111 = gf(i0 + 1, j0 + 1, k0 + 1);

  const T c00 = (T(1) - tx) * c000 + tx * c100;
  const T c10 = (T(1) - tx) * c010 + tx * c110;
  const T c01 = (T(1) - tx) * c001 + tx * c101;
  const T c11 = (T(1) - tx) * c011 + tx * c111;
  const T c0 = (T(1) - ty) * c00 + ty * c10;
  const T c1 = (T(1) - ty) * c01 + ty * c11;
  return (T(1) - tz) * c0 + tz * c1;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL
sample_spherical_linear(SphericalGrid const &grid, CCTK_REAL const *data,
                        CCTK_REAL r, CCTK_REAL theta, CCTK_REAL phi) {
  if (r < 0.0)
    r = 0.0;
  if (r > grid.rmax)
    r = grid.rmax;
  if (theta < 0.0)
    theta = 0.0;
  if (theta > nuX_Utils::sph_grid::kPi)
    theta = nuX_Utils::sph_grid::kPi;
  if (phi < 0.0)
    phi = 0.0;
  if (phi > nuX_Utils::sph_grid::kTwoPi)
    phi = nuX_Utils::sph_grid::kTwoPi;

  int ir0, it0, ip0;
  CCTK_REAL tr, tt, tp;
  {
    const CCTK_REAL q = r / grid.dr;
    if (q <= 0.0) {
      ir0 = 0;
      tr = 0.0;
    } else if (q >= grid.nrad - 1) {
      ir0 = grid.nrad - 2;
      tr = 1.0;
    } else {
      ir0 = static_cast<int>(::floor(q));
      tr = q - ir0;
    }
  }
  {
    const CCTK_REAL q = theta / grid.dtheta;
    if (q <= 0.0) {
      it0 = 0;
      tt = 0.0;
    } else if (q >= grid.ntheta - 1) {
      it0 = grid.ntheta - 2;
      tt = 1.0;
    } else {
      it0 = static_cast<int>(::floor(q));
      tt = q - it0;
    }
  }
  {
    const CCTK_REAL q = phi / grid.dphi;
    if (q <= 0.0) {
      ip0 = 0;
      tp = 0.0;
    } else if (q >= grid.nphi - 1) {
      ip0 = grid.nphi - 2;
      tp = 1.0;
    } else {
      ip0 = static_cast<int>(::floor(q));
      tp = q - ip0;
    }
  }

  const int ir1 = ir0 + 1;
  const int it1 = it0 + 1;
  const int ip1 = ip0 + 1;

  const CCTK_REAL c000 = data[index(grid, ir0, it0, ip0)];
  const CCTK_REAL c100 = data[index(grid, ir1, it0, ip0)];
  const CCTK_REAL c010 = data[index(grid, ir0, it1, ip0)];
  const CCTK_REAL c110 = data[index(grid, ir1, it1, ip0)];
  const CCTK_REAL c001 = data[index(grid, ir0, it0, ip1)];
  const CCTK_REAL c101 = data[index(grid, ir1, it0, ip1)];
  const CCTK_REAL c011 = data[index(grid, ir0, it1, ip1)];
  const CCTK_REAL c111 = data[index(grid, ir1, it1, ip1)];

  const CCTK_REAL c00 = (1.0 - tr) * c000 + tr * c100;
  const CCTK_REAL c10 = (1.0 - tr) * c010 + tr * c110;
  const CCTK_REAL c01 = (1.0 - tr) * c001 + tr * c101;
  const CCTK_REAL c11 = (1.0 - tr) * c011 + tr * c111;
  const CCTK_REAL c0 = (1.0 - tt) * c00 + tt * c10;
  const CCTK_REAL c1 = (1.0 - tt) * c01 + tt * c11;
  return (1.0 - tp) * c0 + tp * c1;
}

CCTK_HOST CCTK_DEVICE inline void
rad_null_flat(CCTK_REAL const r, CCTK_REAL const theta, CCTK_REAL &kt,
              CCTK_REAL &kr, CCTK_REAL &chi, CCTK_REAL &sqrt_det_g) {
  kt = 1.0;
  kr = 1.0;
  chi = 1.0;
  sqrt_det_g = r * r * ::sin(theta);
}

CCTK_HOST CCTK_DEVICE inline void
rad_null(SphericalGrid const &grid, int const irad, int const iray,
         CCTK_REAL const alp, CCTK_REAL const betax, CCTK_REAL const betay,
         CCTK_REAL const betaz, CCTK_REAL const gxx, CCTK_REAL const gxy,
         CCTK_REAL const gxz, CCTK_REAL const gyy, CCTK_REAL const gyz,
         CCTK_REAL const gzz, CCTK_REAL const zvecx, CCTK_REAL const zvecy,
         CCTK_REAL const zvecz, CCTK_INT &mask, CCTK_REAL &kt, CCTK_REAL &kr,
         CCTK_REAL &chi, CCTK_REAL &sqrt_det_g) {
  constexpr CCTK_REAL grid_eps = 0.05;

  CCTK_REAL r, theta, phi;
  nuX_Utils::sph_grid::get_r_theta_phi(grid, irad, iray, r, theta, phi);
  r = (r > grid_eps * grid.dr ? r : grid_eps * grid.dr);
  theta = (theta > grid_eps * grid.dtheta ? theta : grid_eps * grid.dtheta);
  theta = (theta < nuX_Utils::sph_grid::kPi - grid_eps * grid.dtheta
               ? theta
               : nuX_Utils::sph_grid::kPi - grid_eps * grid.dtheta);

  if (mask) {
    rad_null_flat(r, theta, kt, kr, chi, sqrt_det_g);
    return;
  }

  CCTK_REAL Jac[9], InvJac[9];
  nuX_Utils::sph_grid::coord_sph_to_cart_jacobian(r, theta, phi, Jac);
  nuX_Utils::sph_grid::coord_cart_to_sph_jacobian(r, theta, phi, InvJac);

  CCTK_REAL gamma_c[9], gamma_s[9], g_s[16];
  nuX_Utils::metric::space(gxx, gxy, gxz, gyy, gyz, gzz, gamma_c);
  for (int a1 = 0; a1 < 3; ++a1) {
    for (int b1 = 0; b1 < 3; ++b1) {
      gamma_s[3 * a1 + b1] = 0.0;
      for (int a2 = 0; a2 < 3; ++a2)
        for (int b2 = 0; b2 < 3; ++b2)
          gamma_s[3 * a1 + b1] +=
              Jac[3 * a2 + a1] * Jac[3 * b2 + b1] * gamma_c[3 * a2 + b2];
    }
  }

  CCTK_REAL shift_c[3] = {betax, betay, betaz};
  CCTK_REAL shift_s[3];
  matvec3(InvJac, shift_c, shift_s);

  nuX_Utils::metric::spacetime(alp, shift_s[0], shift_s[1], shift_s[2],
                               gamma_s[0], gamma_s[1], gamma_s[2], gamma_s[4],
                               gamma_s[5], gamma_s[8], g_s);

  CCTK_REAL zvec_c[3] = {zvecx, zvecy, zvecz};
  CCTK_REAL W2 = 1.0;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      W2 += gamma_c[3 * a + b] * zvec_c[a] * zvec_c[b];
  if (!(W2 > 0.0)) {
    mask = 1;
    rad_null_flat(r, theta, kt, kr, chi, sqrt_det_g);
    return;
  }
  CCTK_REAL const w_lorentz = ::sqrt(W2);

  CCTK_REAL vel_c[3] = {zvecx / w_lorentz, zvecy / w_lorentz,
                        zvecz / w_lorentz};
  CCTK_REAL vel_s[3];
  matvec3(InvJac, vel_c, vel_s);

  CCTK_REAL u_s[4];
  nuX_Utils::valencia::uvel(alp, shift_s[0], shift_s[1], shift_s[2], w_lorentz,
                            vel_s[0], vel_s[1], vel_s[2], u_s);

  CCTK_REAL const partial_r[4] = {0.0, 1.0, 0.0, 0.0};
  CCTK_REAL const partial_t[4] = {1.0, 0.0, 0.0, 0.0};
  CCTK_REAL e_r[4];

  CCTK_REAL xi = dot4_metric(g_s, partial_r, u_s);
  for (int a = 0; a < 4; ++a)
    e_r[a] = partial_r[a] + xi * u_s[a];

  CCTK_REAL const dir = dot4_metric(g_s, partial_r, e_r);
  CCTK_REAL const e_norm = dot4_metric(g_s, e_r, e_r);
  if (!(e_norm > 0.0)) {
    mask = 1;
    rad_null_flat(r, theta, kt, kr, chi, sqrt_det_g);
    return;
  }
  xi = ::sqrt(e_norm);
  for (int a = 0; a < 4; ++a)
    e_r[a] = e_r[a] / xi * signum(dir);

  CCTK_REAL kappa[4];
  for (int a = 0; a < 4; ++a)
    kappa[a] = u_s[a] + e_r[a];

  kt = kappa[0];
  kr = kappa[1];
  chi = -dot4_metric(g_s, kappa, partial_t);
  if (!(kt > 0.0) || !(kr > 0.0) || !(chi > 0.0)) {
    mask = 1;
    rad_null_flat(r, theta, kt, kr, chi, sqrt_det_g);
    return;
  }

  CCTK_REAL const det_gamma = nuX_Utils::metric::spatial_det(gamma_s);
  if (!(det_gamma > 0.0)) {
    mask = 1;
    rad_null_flat(r, theta, kt, kr, chi, sqrt_det_g);
    return;
  }
  sqrt_det_g = alp * ::sqrt(det_gamma);
}

CCTK_HOST CCTK_DEVICE inline void
evol_density_slice(SphericalGrid const &grid, CCTK_REAL const dt,
                   CCTK_INT const *mask, CCTK_REAL const *theta,
                   CCTK_REAL const *sqrt_det_g, CCTK_REAL const *kt,
                   CCTK_REAL const *R_eff, CCTK_REAL const *abs_coeff,
                   CCTK_REAL const *N_old, CCTK_REAL *N, CCTK_REAL *ndens) {
  CCTK_REAL const lambda = dt / grid.dr;

  N[0] = 0.0;
  ndens[0] = 0.0;

  for (int i = 1; i < grid.nrad; ++i) {
    if (mask[i]) {
      N[i] = 0.0;
      ndens[i] = 0.0;
      continue;
    }

    CCTK_REAL const eta = R_eff[i] * sqrt_det_g[i];
    CCTK_REAL const mu =
        (abs_coeff && kt[i] > 0.0) ? abs_coeff[i] / kt[i] : 0.0;
    CCTK_REAL const a = 1.0 + lambda * theta[i] + dt * mu;
    CCTK_REAL const b = N_old[i] + lambda * theta[i - 1] * N[i - 1] + dt * eta;
    N[i] = a > 0.0 ? b / a : 0.0;

    CCTK_REAL const denom = sqrt_det_g[i] * kt[i];
    ndens[i] = denom > 0.0 ? N[i] / denom : 0.0;
  }
}

CCTK_HOST CCTK_DEVICE inline void
evol_energy_slice(SphericalGrid const &grid, CCTK_REAL const dt,
                  CCTK_INT const *mask, CCTK_REAL const *theta,
                  CCTK_REAL const *kt, CCTK_REAL const *chi,
                  CCTK_REAL const *R_eff, CCTK_REAL const *Q_eff,
                  CCTK_REAL const *ndens, CCTK_REAL const *E_old, CCTK_REAL *E,
                  CCTK_REAL *eave) {
  CCTK_REAL const lambda = dt / grid.dr;
  CCTK_REAL const eps = 1.0e-300;

  if (!mask[0] && R_eff[0] > eps) {
    E[0] = Q_eff[0] * chi[0] / R_eff[0];
    eave[0] = chi[0] > eps ? E[0] / chi[0] : 0.0;
  } else {
    E[0] = 0.0;
    eave[0] = 0.0;
  }

  for (int i = 1; i < grid.nrad; ++i) {
    if (mask[i]) {
      E[i] = 0.0;
      eave[i] = 0.0;
      continue;
    }

    CCTK_REAL const eta = kt[i] > eps ? Q_eff[i] * chi[i] / kt[i] : 0.0;
    CCTK_REAL const mu = kt[i] > eps ? R_eff[i] / kt[i] : 0.0;
    CCTK_REAL const a = ndens[i] * (1.0 + lambda * theta[i]) + dt * mu;
    CCTK_REAL const b =
        ndens[i] * (E_old[i] + lambda * theta[i] * E[i - 1]) + dt * eta;
    E[i] = ::fabs(a) > eps ? b / a : E[i - 1];
    eave[i] = chi[i] > eps ? E[i] / chi[i] : 0.0;
  }
}

} // namespace nuX_M0

#endif
