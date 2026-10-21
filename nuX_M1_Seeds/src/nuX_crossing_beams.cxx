#include <cmath>

#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_seed_utils.hxx"

namespace nuX_M1_Seeds {

using namespace Loop;

extern "C" void nuX_M1_Seeds_CrossingBeamBCs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Seeds_CrossingBeamBCs;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const int ncomponents = radiation_components();
  const CCTK_REAL nx = std::cos(crossing_beam_angle);
  const CCTK_REAL ny = std::sin(crossing_beam_angle);

  // Match the shadow/Kerr inlet: seed the adjacent interior layer as well
  // as the physical ghost cells, so a later Neumann sync retains inflow.
  const auto set_inflow =
      [=] CCTK_DEVICE(const PointDesc &p) {
        if (!(p.BI[0] == -1 || p.NI[0] < 0))
          return;
        CCTK_REAL energy = 0;
        CCTK_REAL flux_y = 0;
        for (int beam = 0; beam < crossing_beam_count; ++beam) {
          const CCTK_REAL sign = beam == 0 ? 1 : -1;
          // Slab beams enter at y=+/-offset and point toward each other.
          const CCTK_REAL transverse =
              nx * (p.y + sign * crossing_beam_offset) -
              sign * ny * (p.x - crossing_beam_inlet_x);
          if (std::abs(transverse) < beam_radius) {
            energy += static_E;
            flux_y += sign * ny * static_E;
          }
        }
        for (int ig = 0; ig < ncomponents; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = rN[i4D] = energy; // Unit mean energy.
          rFx[i4D] = nx * energy;
          rFy[i4D] = flux_y;
          rFz[i4D] = 0;
        }
      };
  grid.loop_int_device<1, 1, 1>(grid.nghostzones, set_inflow);
  grid.loop_bnd_device<1, 1, 1>(grid.nghostzones, set_inflow);
}

} // namespace nuX_M1_Seeds
