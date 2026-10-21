#include <cmath>

#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_seed_utils.hxx"

namespace nuX_M1_Seeds {

using namespace Loop;

extern "C" void nuX_M1_Seeds_ExpandingAbsorberRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Seeds_ExpandingAbsorberRates;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const int ncomponents = radiation_components();

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        CCTK_REAL absorption = 0;
        CCTK_REAL emissivity = 0;
        const CCTK_REAL r2 = p.x * p.x + p.y * p.y + p.z * p.z;
        const CCTK_REAL dx = p.x - absorber_position_x;
        if (r2 < emitter_radius * emitter_radius) {
          // Musolino & Rezzolla Eq. (93): x^2+y^2 inside a 3D sphere.
          absorption = emitter_opacity *
                       std::exp(-16 * (p.x * p.x + p.y * p.y) /
                                (emitter_radius * emitter_radius));
          emissivity = emitter_equilibrium_energy * absorption;
        } else if (dx * dx + p.y * p.y + p.z * p.z <
                   absorber_radius * absorber_radius) {
          absorption = absorber_opacity;
        }

        for (int ig = 0; ig < ncomponents; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          // Identical independent species, with unit mean neutrino energy.
          abs_0[i4D] = abs_1[i4D] = absorption;
          eta_0[i4D] = eta_1[i4D] = emissivity;
          scat_1[i4D] = 0;
          nueave[i4D] = 1;
        }
      });
}

} // namespace nuX_M1_Seeds
