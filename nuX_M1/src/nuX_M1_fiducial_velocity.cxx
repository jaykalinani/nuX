#include <algorithm>
#include <cmath>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "aster_utils.hxx"
#include "nuX_utils.hxx"

#define CGS_GCC (1.619100425158886e-18)

namespace nuX_M1 {

using namespace Arith;
using namespace Loop;
using namespace AsterUtils;

namespace {

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
repair_velocity_and_compute_W(const smat<CCTK_REAL, 3> &g,
                              vec<CCTK_REAL, 3> *const v) {
  for (int a = 0; a < 3; ++a) {
    if (!isfinite((*v)(a))) {
      *v = vec<CCTK_REAL, 3>::pure(CCTK_REAL(0));
      return CCTK_REAL(1);
    }
  }

  const vec<CCTK_REAL, 3> v_low = calc_contraction(g, *v);
  CCTK_REAL v2 = calc_contraction(*v, v_low);
  constexpr CCTK_REAL v2_limit = CCTK_REAL(1) - CCTK_REAL(1.0e-12);
  if (!isfinite(v2) || v2 < CCTK_REAL(0)) {
    *v = vec<CCTK_REAL, 3>::pure(CCTK_REAL(0));
    return CCTK_REAL(1);
  }
  if (v2 >= v2_limit) {
    *v *= sqrt(v2_limit / v2);
    // Recompute from the rounded repaired components so W is consistent with
    // the velocity actually stored, as in the vertex Tmunu path.
    const vec<CCTK_REAL, 3> repaired_v_low = calc_contraction(g, *v);
    v2 = calc_contraction(*v, repaired_v_low);
    if (!isfinite(v2) || v2 < CCTK_REAL(0) || v2 >= CCTK_REAL(1)) {
      *v = vec<CCTK_REAL, 3>::pure(CCTK_REAL(0));
      v2 = CCTK_REAL(0);
    }
  }
  return CCTK_REAL(1) / sqrt(CCTK_REAL(1) - v2);
}

} // namespace

extern "C" void nuX_M1_FiducialVelocity(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_FiducialVelocity;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_FiducialVelocity");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  const smat<GF3D2<const CCTK_REAL8>, 3> gf_g{
      GF3D2<const CCTK_REAL8>(layout_vc, gxx),
      GF3D2<const CCTK_REAL8>(layout_vc, gxy),
      GF3D2<const CCTK_REAL8>(layout_vc, gxz),
      GF3D2<const CCTK_REAL8>(layout_vc, gyy),
      GF3D2<const CCTK_REAL8>(layout_vc, gyz),
      GF3D2<const CCTK_REAL8>(layout_vc, gzz)};

  if (CCTK_Equals(fiducial_velocity, "fluid")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const int ijk = layout_cc.linear(p.i, p.j, p.k);

          const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
            return calc_avg_v2c(gf_g(i, j), p);
          });

          vec<CCTK_REAL, 3> v_up{velx[ijk], vely[ijk], velz[ijk]};
          fidu_w_lorentz[ijk] = repair_velocity_and_compute_W(g_avg, &v_up);
          fidu_velx[ijk] = v_up(0);
          fidu_vely[ijk] = v_up(1);
          fidu_velz[ijk] = v_up(2);
        });

  } else if (CCTK_Equals(fiducial_velocity, "mixed")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const int ijk = layout_cc.linear(p.i, p.j, p.k);

          const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
            return calc_avg_v2c(gf_g(i, j), p);
          });

          // Weight continuously between the fluid velocity and zero in the
          // atmosphere, then enforce the same timelike-domain invariant used
          // by the vertex stress-energy construction.
          const CCTK_REAL local_dens = isfinite(dens[ijk])
                                           ? fmax(dens[ijk], CCTK_REAL(0))
                                           : CCTK_REAL(0);
          const CCTK_REAL transition_dens =
              fmax(fiducial_velocity_rho_fluid * CGS_GCC, CCTK_REAL(0));
          const CCTK_REAL denom = fmax(local_dens, transition_dens);
          const CCTK_REAL weight =
              denom > CCTK_REAL(0) ? local_dens / denom : CCTK_REAL(0);
          vec<CCTK_REAL, 3> v_up{weight * velx[ijk], weight * vely[ijk],
                                 weight * velz[ijk]};
          fidu_w_lorentz[ijk] = repair_velocity_and_compute_W(g_avg, &v_up);
          fidu_velx[ijk] = v_up(0);
          fidu_vely[ijk] = v_up(1);
          fidu_velz[ijk] = v_up(2);
        });

  } else {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const int ijk = layout_cc.linear(p.i, p.j, p.k);

          fidu_velx[ijk] = CCTK_REAL(0);
          fidu_vely[ijk] = CCTK_REAL(0);
          fidu_velz[ijk] = CCTK_REAL(0);
          fidu_w_lorentz[ijk] = CCTK_REAL(1);
        });
  }
}

} // namespace nuX_M1
