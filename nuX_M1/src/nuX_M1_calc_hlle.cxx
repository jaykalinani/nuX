#include <cassert>
#include <cmath>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define NUX_M1_CLOSURE_IMPLEMENTATION
#include "nuX_M1_hlle.hxx"

namespace nuX_M1 {
namespace {

// Use nuX's vertex-to-cell interpolation for every stencil metric. In
// particular, do not reinterpret vertex-centred ADM arrays as cell-centred.
template <int dir>
CCTK_HOST CCTK_DEVICE inline tensor::metric<3>
hlle_cell_metric(const tensor::slicing_geometry_const &geom,
                 const Loop::GF3D2layout &layout_vc, const Loop::PointDesc &p,
                 const int offset) {
  tensor::metric<3> gamma;
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b) {
      const Loop::GF3D2<const CCTK_REAL> gf(
          layout_vc, geom.get_space_metric_comp(a, b));
      gamma(a, b) = tensor::interp_v2c(gf, p, offset * (dir == 0),
                                      offset * (dir == 1), offset * (dir == 2));
    }
  return gamma;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL
hlle_primitive(const int component, const int index, const CCTK_REAL volume,
               const CCTK_REAL *rN, const CCTK_REAL *rFx,
               const CCTK_REAL *rFy, const CCTK_REAL *rFz,
               const CCTK_REAL *rE, const CCTK_REAL E_floor,
               const CCTK_REAL N_floor) {
  const CCTK_REAL raw_E = rE[index];
  // fmax and the slope limiter can hide NaNs; check raw stencil values.
  if (!isfinite(raw_E))
    device_abort();
  const CCTK_REAL E = fmax(raw_E, E_floor);
  if (!(isfinite(volume) && volume > 0 && isfinite(E) && E >= 0))
    device_abort();
  if (E == 0) {
    // A true vacuum has defined zero ratios by its regular limit. A positive
    // N or F at zero energy has no finite N/E, Fi/E representation.
    if (rN[index] != 0 || rFx[index] != 0 || rFy[index] != 0 || rFz[index] != 0 ||
        N_floor != 0)
      device_abort();
    return 0;
  }
  CCTK_REAL primitive;
  if (component == HLLE_E) {
    primitive = E / volume;
  } else if (component == HLLE_N) {
    const CCTK_REAL N = rN[index];
    if (!isfinite(N))
      device_abort();
    primitive = fmax(N, N_floor) / E;
  } else {
    const CCTK_REAL *F = component == HLLE_FX ? rFx
                         : component == HLLE_FY ? rFy : rFz;
    const CCTK_REAL flux = F[index];
    if (!isfinite(flux))
      device_abort();
    primitive = flux / E;
  }
  if (!isfinite(primitive))
    device_abort();
  return primitive;
}

template <int dir> void UpdateRHSFromHLLE(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_UpdateRHSFromHLLE;
  DECLARE_CCTK_PARAMETERS;
  using namespace Loop;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  const tensor::slicing_geometry_const geom(
      layout_vc, layout_cc, alp, betax, betay, betaz, gxx, gxy, gxz, gyy,
      gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz);
  const CCTK_REAL dx = CCTK_DELTA_SPACE(dir);
  const CCTK_REAL theta = CCTK_Equals(hlle_limiter, "mc") ? 2 : 1;
  const int groupspec = ngroups * nspecies;
  CCTK_REAL *r_rhs[5] = {rN_rhs, rFx_rhs, rFy_rhs, rFz_rhs, rE_rhs};
  closure_t closure_fun = CLOSURE_MINERBO;
  if (CCTK_Equals(closure, "Eddington"))
    closure_fun = CLOSURE_EDDINGTON;
  else if (CCTK_Equals(closure, "Kershaw"))
    closure_fun = CLOSURE_KERSHAW;
  else if (CCTK_Equals(closure, "thin"))
    closure_fun = CLOSURE_THIN;

  assert(cctk_nghostzones[dir] >= 2);
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int center = layout_cc.linear(p.i, p.j, p.k);
        if (nuX_m1_mask[center])
          return;
        // Each cell independently evaluates the same two shared-face fluxes.
        // Only one face and one group are live; no additional grid functions
        // or full-grid temporary face arrays are needed.
        for (int face = 0; face < 2; ++face) {
          const int left = face - 1;
          int indices[4];
          CCTK_REAL volumes[4];
          tensor::metric<3> gamma;
          for (int s = 0; s < 4; ++s) {
            const int offset = left + s - 1;
            indices[s] = layout_cc.linear(p.i + offset * (dir == 0),
                                           p.j + offset * (dir == 1),
                                           p.k + offset * (dir == 2));
            const auto gc = hlle_cell_metric<dir>(geom, layout_vc, p, offset);
            volumes[s] = sqrt(gc.det());
            if (!(isfinite(volumes[s]) && volumes[s] > 0))
              device_abort();
            if (s == 1 || s == 2)
              for (int a = 0; a < 3; ++a)
                for (int b = a; b < 3; ++b)
                  gamma(a, b) = s == 1 ? 0.5 * gc(a, b)
                                         : gamma(a, b) + 0.5 * gc(a, b);
          }
          tensor::generic<CCTK_REAL, 4, 1> beta_l, beta_r;
          geom.get_shift_vec(p, left * (dir == 0), left * (dir == 1),
                              left * (dir == 2), &beta_l);
          geom.get_shift_vec(p, (left + 1) * (dir == 0),
                              (left + 1) * (dir == 1),
                              (left + 1) * (dir == 2), &beta_r);
          for (int a = 0; a < 4; ++a)
            beta_l(a) = 0.5 * beta_l(a) + 0.5 * beta_r(a);
          const CCTK_REAL alpha =
              0.5 * geom.get_lapse(p, left * (dir == 0), left * (dir == 1),
                                     left * (dir == 2)) +
              0.5 * geom.get_lapse(p, (left + 1) * (dir == 0),
                                     (left + 1) * (dir == 1),
                                     (left + 1) * (dir == 2));
          hlle_geometry_t g;
          g.set(alpha, beta_l, gamma,
                0.5 * fidu_velx[indices[1]] + 0.5 * fidu_velx[indices[2]],
                0.5 * fidu_vely[indices[1]] + 0.5 * fidu_vely[indices[2]],
                0.5 * fidu_velz[indices[1]] + 0.5 * fidu_velz[indices[2]]);
          // Distance between coordinate faces measured along their unit
          // normal: dx/sqrt(gamma^{dd}), including non-diagonal metrics.
          const CCTK_REAL proper_dx = dx / sqrt(g.gamma_uu(dir, dir));

          for (int ig = 0; ig < groupspec; ++ig) {
            int ig_indices[4];
            for (int s = 0; s < 4; ++s) {
              const int offset = left + s - 1;
              ig_indices[s] = layout_cc.linear(p.i + offset * (dir == 0),
                                               p.j + offset * (dir == 1),
                                               p.k + offset * (dir == 2), ig);
            }
            CCTK_REAL ql[5], qr[5];
            for (int a = 0; a < 5; ++a) {
              CCTK_REAL q[4];
              for (int s = 0; s < 4; ++s)
                q[s] = hlle_primitive(a, ig_indices[s], volumes[s], rN, rFx,
                                       rFy, rFz, rE, rad_E_floor, rad_N_floor);
              ql[a] = q[1] + 0.5 * hlle_slope(q[0], q[1], q[2], theta);
              qr[a] = q[2] - 0.5 * hlle_slope(q[1], q[2], q[3], theta);
            }
            const auto l = hlle_state(cctkGH, p.i, p.j, p.k, ig, dir, g, ql,
                                       closure_fun, rad_E_floor, rad_N_floor,
                                       rad_eps, closure_epsilon,
                                       closure_maxiter, use_fallback != 0);
            const auto r = hlle_state(cctkGH, p.i, p.j, p.k, ig, dir, g, qr,
                                       closure_fun, rad_E_floor, rad_N_floor,
                                       rad_eps, closure_epsilon,
                                       closure_maxiter, use_fallback != 0);
            const CCTK_REAL sm = fmin(CCTK_REAL(0), fmin(l.lm, r.lm));
            const CCTK_REAL sp = fmax(CCTK_REAL(0), fmax(l.lp, r.lp));
            // N's transport current is tied to H/J; its spatial relaxation
            // is controlled by the energy transport opacity, not abs_0 alone.
            const CCTK_REAL A = hlle_diffusion_factor(
                abs_1[ig_indices[1]] + scat_1[ig_indices[1]],
                abs_1[ig_indices[2]] + scat_1[ig_indices[2]], proper_dx);
            const int out = layout_cc.linear(p.i, p.j, p.k, ig);
            const CCTK_REAL sign_over_dx = face == 0 ? 1 / dx : -1 / dx;
            for (int a = 0; a < 5; ++a)
              r_rhs[a][out] += sign_over_dx *
                  hlle_flux(l.u[a], r.u[a], l.flux[a], r.flux[a], sm, sp, A,
                             a >= HLLE_FX && a <= HLLE_FZ);
          }
        }
      });
}

} // namespace

extern "C" void nuX_M1_UpdateRHSFromHLLE(CCTK_ARGUMENTS) {
  UpdateRHSFromHLLE<0>(cctkGH);
  UpdateRHSFromHLLE<1>(cctkGH);
  UpdateRHSFromHLLE<2>(cctkGH);
}

} // namespace nuX_M1
