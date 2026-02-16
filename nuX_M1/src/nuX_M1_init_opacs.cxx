#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "loop_device.hxx"
#include "vect.hxx"

namespace nuX_M1 {

using namespace Loop;

extern "C" void nuX_M1_InitOpacs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int groupspec = ngroups * nspecies;

  const auto &grid = GridDescBaseDevice(cctkGH);
  GF3D2layout layout(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout.linear(p.i, p.j, p.k, ig);
        }
      });
}

} // namespace nuX_M1
