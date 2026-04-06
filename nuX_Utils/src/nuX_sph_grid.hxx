#ifndef NUX_SPH_GRID_HXX
#define NUX_SPH_GRID_HXX

#include <cassert>
#include <cmath>

#include <cctk.h>
#include <loop_device.hxx>

namespace nuX_Utils {
namespace sph_grid {

constexpr CCTK_REAL kPi = CCTK_REAL(3.14159265358979323846264338327950288L);
constexpr CCTK_REAL kTwoPi = CCTK_REAL(2.0) * kPi;
constexpr CCTK_REAL kHalfPi = CCTK_REAL(0.5) * kPi;

struct SphericalGrid {
  CCTK_REAL origin[3] = {0.0, 0.0, 0.0};
  CCTK_REAL rmax = 0.0;
  int nrad = 0;
  int nray = 0;
  int ntheta = 0;
  int nphi = 0;
  bool endpoint = false;
  CCTK_REAL dr = 0.0;
  CCTK_REAL dtheta = 0.0;
  CCTK_REAL dphi = 0.0;
};

CCTK_HOST CCTK_DEVICE inline int index(SphericalGrid const &grid,
                                       int const irad, int const iray) {
  return iray * grid.nrad + irad;
}

CCTK_HOST CCTK_DEVICE inline void
init(SphericalGrid &grid, CCTK_REAL const origin[3], CCTK_REAL const rmax,
     int const nrad, int const ntheta, int const nphi, bool const endpoint) {
  assert(rmax > 0.0);
  assert(nrad > 1);
  assert(ntheta > 2);
  assert(nphi > 5);

  grid.origin[0] = origin[0];
  grid.origin[1] = origin[1];
  grid.origin[2] = origin[2];
  grid.rmax = rmax;
  grid.nrad = nrad;
  grid.nray = ntheta * nphi;
  grid.ntheta = ntheta;
  grid.nphi = nphi;
  grid.endpoint = endpoint;
  grid.dr = rmax / (nrad - 1);
  grid.dtheta = kPi / (ntheta - (endpoint ? 1 : 0));
  grid.dphi = kTwoPi / (nphi - (endpoint ? 1 : 0));
}

CCTK_HOST CCTK_DEVICE inline SphericalGrid
make_grid(CCTK_REAL const origin[3], CCTK_REAL const rmax, int const nrad,
          int const ntheta, int const nphi, bool const endpoint) {
  SphericalGrid grid;
  init(grid, origin, rmax, nrad, ntheta, nphi, endpoint);
  return grid;
}

CCTK_HOST CCTK_DEVICE inline int get_iray(SphericalGrid const &grid,
                                          int const itheta, int const iphi) {
  return itheta + iphi * grid.ntheta;
}

CCTK_HOST CCTK_DEVICE inline int get_itheta(SphericalGrid const &grid,
                                            int const iray) {
  assert(iray >= 0 && iray < grid.nray);
  return iray % grid.ntheta;
}

CCTK_HOST CCTK_DEVICE inline int get_iphi(SphericalGrid const &grid,
                                          int const iray) {
  assert(iray >= 0 && iray < grid.nray);
  return iray / grid.ntheta;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL get_r(SphericalGrid const &grid,
                                             int const irad) {
  assert(irad >= 0 && irad < grid.nrad);
  return grid.dr * irad;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL get_theta(SphericalGrid const &grid,
                                                 int const iray) {
  return grid.dtheta * (1.0 + 0.5 * (grid.endpoint ? 0 : 1)) *
         get_itheta(grid, iray);
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL get_phi(SphericalGrid const &grid,
                                               int const iray) {
  return grid.dphi * get_iphi(grid, iray);
}

CCTK_HOST CCTK_DEVICE inline void
coord_sph_to_cart(CCTK_REAL const r, CCTK_REAL const theta, CCTK_REAL const phi,
                  CCTK_REAL &x, CCTK_REAL &y, CCTK_REAL &z) {
  assert(r >= 0.0);
  assert(theta >= 0.0 && theta <= kPi);
  assert(phi >= 0.0 && phi <= kTwoPi);
  x = r * ::sin(theta) * ::cos(phi);
  y = r * ::sin(theta) * ::sin(phi);
  z = r * ::cos(theta);
}

CCTK_HOST CCTK_DEVICE inline void
coord_cart_to_sph(CCTK_REAL const x, CCTK_REAL const y, CCTK_REAL const z,
                  CCTK_REAL &r, CCTK_REAL &theta, CCTK_REAL &phi) {
  r = ::sqrt(x * x + y * y + z * z);
  if (r > CCTK_REAL(0.0)) {
    theta = ::acos(z / r);
    phi = ::atan2(y, x);
    if (phi < 0.0) {
      phi += kTwoPi;
    }
  } else {
    theta = kHalfPi;
    phi = 0.0;
  }
}

CCTK_HOST CCTK_DEVICE inline void
coord_sph_to_cart_jacobian(CCTK_REAL const r, CCTK_REAL const theta,
                           CCTK_REAL const phi, CCTK_REAL J[9]) {
  CCTK_REAL const sin_theta = ::sin(theta);
  CCTK_REAL const cos_theta = ::cos(theta);
  CCTK_REAL const sin_phi = ::sin(phi);
  CCTK_REAL const cos_phi = ::cos(phi);

  J[0] = sin_theta * cos_phi;
  J[1] = r * cos_theta * cos_phi;
  J[2] = -r * sin_theta * sin_phi;
  J[3] = sin_theta * sin_phi;
  J[4] = r * cos_theta * sin_phi;
  J[5] = r * sin_theta * cos_phi;
  J[6] = cos_theta;
  J[7] = -r * sin_theta;
  J[8] = 0.0;
}

CCTK_HOST CCTK_DEVICE inline void
coord_cart_to_sph_jacobian(CCTK_REAL const r, CCTK_REAL const theta,
                           CCTK_REAL const phi, CCTK_REAL J[9]) {
  CCTK_REAL const safe_r = (r > CCTK_REAL(1.0e-300) ? r : CCTK_REAL(1.0e-300));
  CCTK_REAL const ir = 1.0 / safe_r;
  CCTK_REAL const sin_theta = ::sin(theta);
  CCTK_REAL const cos_theta = ::cos(theta);
  CCTK_REAL const abs_sin_theta = ::fabs(sin_theta);
  CCTK_REAL const safe_sin_theta =
      (abs_sin_theta > CCTK_REAL(1.0e-300) ? abs_sin_theta
                                           : CCTK_REAL(1.0e-300));
  CCTK_REAL const csc_theta = 1.0 / safe_sin_theta;
  CCTK_REAL const sin_phi = ::sin(phi);
  CCTK_REAL const cos_phi = ::cos(phi);

  J[0] = sin_theta * cos_phi;
  J[1] = sin_theta * sin_phi;
  J[2] = cos_theta;
  J[3] = ir * cos_theta * cos_phi;
  J[4] = ir * cos_theta * sin_phi;
  J[5] = -ir * sin_theta;
  J[6] = -ir * csc_theta * sin_phi;
  J[7] = ir * csc_theta * cos_phi;
  J[8] = 0.0;
}

CCTK_HOST CCTK_DEVICE inline void
get_r_theta_phi(SphericalGrid const &grid, int const irad, int const iray,
                CCTK_REAL &r, CCTK_REAL &theta, CCTK_REAL &phi) {
  r = get_r(grid, irad);
  theta = get_theta(grid, iray);
  phi = get_phi(grid, iray);
}

CCTK_HOST CCTK_DEVICE inline void get_x_y_z(SphericalGrid const &grid,
                                            int const irad, int const iray,
                                            CCTK_REAL &x, CCTK_REAL &y,
                                            CCTK_REAL &z) {
  CCTK_REAL r, theta, phi;
  get_r_theta_phi(grid, irad, iray, r, theta, phi);
  coord_sph_to_cart(r, theta, phi, x, y, z);
  x += grid.origin[0];
  y += grid.origin[1];
  z += grid.origin[2];
}

} // namespace sph_grid
} // namespace nuX_Utils

#endif
