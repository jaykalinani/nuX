#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "loop_device.hxx"
#include "vect.hxx"

namespace nuX_M1 {

using namespace Loop;

extern "C" void nuX_M1_InitFluxesRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int groupspec = ngroups * nspecies;
  const int flux_comps = 5 * groupspec; 

  const auto &grid = GridDescBaseDevice(cctkGH);
  GF3D2layout layout(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                [=] CCTK_DEVICE(const PointDesc &p) {
    for (int ig = 0; ig < groupspec; ++ig) {
      const int i4D = layout.linear(p.i, p.j, p.k, ig);
      rN_rhs[i4D]  = 0.0;
      rE_rhs[i4D]  = 0.0;
      rFx_rhs[i4D] = 0.0;
      rFy_rhs[i4D] = 0.0;
      rFz_rhs[i4D] = 0.0;
    }
  });

  for (int dir = 0; dir < 3; ++dir) {
    CCTK_REAL *nu_flux_dir = (dir == 0 ? nu_flux_x : (dir == 1 ? nu_flux_y : nu_flux_z));

    grid.loop_all_device<1, 1, 1>(grid.nghostzones,
                                  [=] CCTK_DEVICE(const PointDesc &p) {
      for (int comp = 0; comp < flux_comps; ++comp) {
        const int idx = layout.linear(p.i, p.j, p.k, comp);
        nu_flux_dir[idx] = 0.0;
      }
    });
  }
}

} // namespace nuX_M1

