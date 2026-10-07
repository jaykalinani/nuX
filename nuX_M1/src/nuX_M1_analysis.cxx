#include <algorithm>
#include <cassert>
#include <cmath>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_baryon_mass.hxx"
#include "nuX_utils.hxx"

namespace nuX_M1 {

using namespace std;
using namespace Loop;
using namespace nuX_Utils;

extern "C" void nuX_M1_Analysis(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Analysis;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_Analysis");
  }

  // particle_mass is in MeV
  CCTK_REAL const mb = nuX_Utils::AverageBaryonMass(particle_mass);

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
        const CCTK_REAL detg = nuX_Utils::metric::spatial_det(
            gxx_cc, gxy_cc, gxz_cc, gyy_cc, gyz_cc, gzz_cc);
        const CCTK_REAL volform_ijk =
            isfinite(detg) && detg > 0.0 ? sqrt(detg) : 0.0;

        const CCTK_REAL nb =
            isfinite(rho[ijk]) && rho[ijk] > 0.0 && isfinite(mb) && mb > 0.0
                ? rho[ijk] / mb
                : 0.0;
        const CCTK_REAL inv_baryon_volume =
            nb > 0.0 && volform_ijk > 0.0 ? 1.0 / (volform_ijk * nb) : 0.0;
        ynue[ijk] = rnnu[ijke] * inv_baryon_volume;
        ynua[ijk] = rnnu[ijka] * inv_baryon_volume;
        ynux[ijk] = rnnu[ijkx] * inv_baryon_volume;

        CCTK_REAL const egas = rho[ijk] * (1 + eps[ijk]);
        const CCTK_REAL inv_volform =
            volform_ijk > 0.0 ? 1.0 / volform_ijk : 0.0;
        CCTK_REAL const enue = rJ[ijke] * inv_volform;
        CCTK_REAL const enua = rJ[ijka] * inv_volform;
        CCTK_REAL const enux = rJ[ijkx] * inv_volform;
        CCTK_REAL const etot = egas + enue + enua + enux;
        const CCTK_REAL inv_etot =
            isfinite(etot) && etot > 0.0 ? 1.0 / etot : 0.0;
        znue[ijk] = enue * inv_etot;
        znua[ijk] = enua * inv_etot;
        znux[ijk] = enux * inv_etot;
      });
}

} // namespace nuX_M1
