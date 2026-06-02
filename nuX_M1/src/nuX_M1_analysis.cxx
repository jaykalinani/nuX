#include <algorithm>
#include <cassert>
#include <sstream>
#include <cmath>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M1_macro.hxx"
#include "nuX_utils.hxx"
#include "avg_baryon_mass.hpp"

namespace nuX_M1 {

using namespace std;
using namespace Loop;
using namespace nuX_Utils;
using namespace nuX_Rates;

extern "C" void nuX_M1_Analysis(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Analysis;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_Analysis");
  }

  // particle_mass is in MeV
  CCTK_REAL const mb = AverageBaryonMass(particle_mass);

  assert(nspecies == 3);
  assert(ngroups == 1);

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  const GF3D2<const CCTK_REAL> gf_gxx(layout_vc, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(layout_vc, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(layout_vc, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(layout_vc, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(layout_vc, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(layout_vc, gzz);

  // UTILS_LOOP3(nuX_m1_analysis, k, 0, cctk_lsh[2], j, 0, cctk_lsh[1], i, 0,
  // 						cctk_lsh[0]) {
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);

        int const ijke = layout_cc.linear(p.i, p.j, p.k, 0);
        int const ijka = layout_cc.linear(p.i, p.j, p.k, 1);
        int const ijkx = layout_cc.linear(p.i, p.j, p.k, 2);
        const CCTK_REAL gxx_cc = tensor::interp_v2c(gf_gxx, p);
        const CCTK_REAL gxy_cc = tensor::interp_v2c(gf_gxy, p);
        const CCTK_REAL gxz_cc = tensor::interp_v2c(gf_gxz, p);
        const CCTK_REAL gyy_cc = tensor::interp_v2c(gf_gyy, p);
        const CCTK_REAL gyz_cc = tensor::interp_v2c(gf_gyz, p);
        const CCTK_REAL gzz_cc = tensor::interp_v2c(gf_gzz, p);
        const CCTK_REAL volform_ijk = sqrt(nuX_Utils::metric::spatial_det(
            gxx_cc, gxy_cc, gxz_cc, gyy_cc, gyz_cc, gzz_cc));

        CCTK_REAL const nb = rho[ijk] / mb;
        ynue[ijk] = rnnu[ijke] / volform_ijk / nb;
        ynua[ijk] = rnnu[ijka] / volform_ijk / nb;
        ynux[ijk] = rnnu[ijkx] / volform_ijk / nb;

        CCTK_REAL const egas = rho[ijk] * (1 + eps[ijk]);
        CCTK_REAL const enue = rJ[ijke] / volform_ijk;
        CCTK_REAL const enua = rJ[ijka] / volform_ijk;
        CCTK_REAL const enux = rJ[ijkx] / volform_ijk;
        CCTK_REAL const etot = egas + enue + enua + enux;
        znue[ijk] = enue / etot;
        znua[ijk] = enua / etot;
        znux[ijk] = enux / etot;
      });
}

} // namespace nuX_M1
