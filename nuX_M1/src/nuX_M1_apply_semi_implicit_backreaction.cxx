#include <cassert>
#include <cmath>

#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_baryon_mass.hxx"

namespace nuX_M1 {

using namespace Loop;

extern "C" void nuX_M1_ApplySemiImplicitBackreaction(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_ApplySemiImplicitBackreaction;
  DECLARE_CCTK_PARAMETERS;

  if (!CCTK_Equals(method, "semi-implicit") || !backreact ||
      *semi_implicit_stage != 1)
    return;

  if (ngroups != 1 || nspecies != 3)
    CCTK_ERROR("nuX_M1::backreact requires ngroups=1 and nspecies=3");

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  assert(dt > 0.0);

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const CCTK_REAL mb = nuX_Utils::AverageBaryonMass(particle_mass);

  // CalcUpdate stores the accepted physical collision increment divided by
  // this implicit step's dt in source_update_delta_*. Numerical radiation
  // floors and realizability repairs are deliberately not included.
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        if (nuX_m1_mask[ijk])
          return;

        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          const CCTK_REAL delta_N = dt * source_update_delta_N[i4D];
          const CCTK_REAL delta_E = dt * source_update_delta_E[i4D];
          const CCTK_REAL delta_Fx = dt * source_update_delta_Fx[i4D];
          const CCTK_REAL delta_Fy = dt * source_update_delta_Fy[i4D];
          const CCTK_REAL delta_Fz = dt * source_update_delta_Fz[i4D];

          assert(std::isfinite(delta_N));
          assert(std::isfinite(delta_E));
          assert(std::isfinite(delta_Fx));
          assert(std::isfinite(delta_Fy));
          assert(std::isfinite(delta_Fz));

          momx[ijk] -= delta_Fx;
          momy[ijk] -= delta_Fy;
          momz[ijk] -= delta_Fz;
          tau[ijk] -= delta_E;
          DYe[ijk] +=
              -mb * ((ig == 0 ? delta_N : 0.0) -
                     (ig == 1 ? delta_N : 0.0));

          assert(std::isfinite(momx[ijk]));
          assert(std::isfinite(momy[ijk]));
          assert(std::isfinite(momz[ijk]));
          assert(std::isfinite(tau[ijk]));
          assert(std::isfinite(DYe[ijk]));
        }
      });
}

} // namespace nuX_M1
