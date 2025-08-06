#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace RatesToy {
using namespace Loop;

extern "C" void RatesToy_Init(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_RatesToy_Init;
  DECLARE_CCTK_PARAMETERS

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

			for (int ig = 0; ig < ngroups * nspecies; ++ig) {
				int const i4D = layout2.linearVec3(p.i, p.j, p.k, ig);

				rnt[i4D] = set_n;
				rJt[i4D] = set_J;
				rHt_t[i4D] = set_Ht;
				rHxt[i4D] = set_Hx;
				rHyt[i4D] = set_Hy;
				rHzt[i4D] = set_Hz;
				chit[i4D] = set_chi;
			}
		});
}

} // namespace
