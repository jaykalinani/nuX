#include <cstring>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define PINDEX1D(ig, iv) ((iv) + (ig) * 5)

namespace nuX_M1 {

using namespace Loop;

extern "C" void nuX_M1_InitFluxesRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_InitFluxesRHS;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_InitFluxesRHS");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  const int GFNSIZE = (cctk_lsh[0]-1) * (cctk_lsh[1]-1) * (cctk_lsh[2]-1);
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int groupspec = ngroups * nspecies; // total components

        for (int ig = 0; ig < groupspec; ++ig) {
          const int ijk = layout2.linear(p.i, p.j, p.k);
          const int i4D = layout2.linear(p.i, p.j, p.k, ig);

          // Zero the fluxes
          for (int iv = 0; iv < 5; ++iv) {
            int comp = PINDEX1D(ig, iv);
            nu_flux_x[comp * GFNSIZE + ijk] = 0.0;
            nu_flux_y[comp * GFNSIZE + ijk] = 0.0;
            nu_flux_z[comp * GFNSIZE + ijk] = 0.0;
          }

          // Zero the RHS buffers
          rN_rhs[i4D] = 0.0;
          rE_rhs[i4D] = 0.0;
          rFx_rhs[i4D] = 0.0;
          rFy_rhs[i4D] = 0.0;
          rFz_rhs[i4D] = 0.0;

        }
      });
}

} // namespace nuX_M1
