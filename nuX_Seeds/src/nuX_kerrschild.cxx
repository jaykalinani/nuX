#include <cmath>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_tensor.hxx"
#include "setup_eos.hxx"

namespace nuX_Seeds {

using namespace Loop;
using namespace EOSX;
using namespace nuX_Utils;

extern "C" void nuX_Seeds_SetupHydroTest_kerrschild(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupHydroTest_kerrschild;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupHydroTest_kerrschild");

  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout3(cctkGH, {1, 0, 0});
  const GF3D2layout layout4(cctkGH, {0, 1, 0});
  const GF3D2layout layout5(cctkGH, {0, 0, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        rho[ijk] = static_rho;
        eps[ijk] = static_eps;
        velx[ijk] = static_velx;
        vely[ijk] = static_vely;
        velz[ijk] = static_velz;
        Ye[ijk] = static_ye;

        press[ijk] =
            eos_3p_ig->press_from_rho_eps_ye(rho[ijk], eps[ijk], Ye[ijk]);
      });
  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout3.linear(p.i, p.j, p.k);
        Avec_x[ijk] = 0.;
      });

  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout4.linear(p.i, p.j, p.k);
        Avec_y[ijk] = 0.;
      });

  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout5.linear(p.i, p.j, p.k);
        Avec_z[ijk] = 0.;
      });
}

extern "C" void nuX_Seeds_SetupNeutTest_kerrschild(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupNeutTest_kerrschild;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupNeutTest_kerrschild");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = 0.0;
          rFx[i4D] = 0.0;
          rFy[i4D] = 0.0;
          rFz[i4D] = 0.0;
          rN[i4D] = 0.0;
        }
      });
}

extern "C" void nuX_Seeds_KerrInflow(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_KerrInflow;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  const int nghy = cctk_nghostzones[1];
  const GF3D2<const CCTK_REAL> gf_alp(layout_vc, alp);
  const GF3D2<const CCTK_REAL> gf_betax(layout_vc, betax);
  const GF3D2<const CCTK_REAL> gf_betay(layout_vc, betay);
  const GF3D2<const CCTK_REAL> gf_betaz(layout_vc, betaz);
  const GF3D2<const CCTK_REAL> gf_gxx(layout_vc, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(layout_vc, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(layout_vc, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(layout_vc, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(layout_vc, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(layout_vc, gzz);

  // Apply faces in THC order to get deterministic edge/corner precedence:
  // -X, +X, -Y, +Y, -Z, +Z.
  //
  // This must be split into multiple passes (instead of one combined kernel):
  // Y-face copy BCs read neighboring cells, and doing that while writing all
  // faces in one kernel causes boundary read/write races on GPU.

  // -X: inject Kerr beam, vacuum outside beam window
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (p.BI[0] != -1.0)
          return;

        const bool in_beam = (fabs(p.y) <= kerr_beam_width) &&
                             (p.z >= kerr_beam_position) &&
                             (p.z <= kerr_beam_position + kerr_beam_width);
        CCTK_REAL E_new = 0.0;
        CCTK_REAL Fx_new = 0.0;
        CCTK_REAL Fy_new = 0.0;
        CCTK_REAL Fz_new = 0.0;
        CCTK_REAL N_new = 0.0;

        if (in_beam) {
          CCTK_REAL const eps_flux = 0.01;
          CCTK_REAL const alp_ijk = tensor::interp_v2c(gf_alp, p);

          CCTK_REAL const g_xx = tensor::interp_v2c(gf_gxx, p);
          CCTK_REAL const g_xy = tensor::interp_v2c(gf_gxy, p);
          CCTK_REAL const g_xz = tensor::interp_v2c(gf_gxz, p);
          CCTK_REAL const g_yy = tensor::interp_v2c(gf_gyy, p);
          CCTK_REAL const g_yz = tensor::interp_v2c(gf_gyz, p);
          CCTK_REAL const g_zz = tensor::interp_v2c(gf_gzz, p);

          CCTK_REAL const beta_u_x = tensor::interp_v2c(gf_betax, p);
          CCTK_REAL const beta_u_y = tensor::interp_v2c(gf_betay, p);
          CCTK_REAL const beta_u_z = tensor::interp_v2c(gf_betaz, p);

          CCTK_REAL const beta_d_x =
              g_xx * beta_u_x + g_xy * beta_u_y + g_xz * beta_u_z;
          CCTK_REAL const beta_d_y =
              g_xy * beta_u_x + g_yy * beta_u_y + g_yz * beta_u_z;
          CCTK_REAL const beta_d_z =
              g_xz * beta_u_x + g_yz * beta_u_y + g_zz * beta_u_z;
          CCTK_REAL const beta2 =
              beta_u_x * beta_d_x + beta_u_y * beta_d_y + beta_u_z * beta_d_z;

          CCTK_REAL const disc = beta_d_x * beta_d_x - beta2 +
                                 alp_ijk * alp_ijk * (1.0 - eps_flux);
          CCTK_REAL a = 0.0;
          if (g_xx > 0.0) {
            a = (-beta_d_x + sqrt(fmax(disc, 0.0))) / g_xx;
          }

          CCTK_REAL const detg = g_xx * (g_yy * g_zz - g_yz * g_yz) -
                                 g_xy * (g_xy * g_zz - g_yz * g_xz) +
                                 g_xz * (g_xy * g_yz - g_yy * g_xz);
          E_new = sqrt(fmax(detg, 0.0));

          CCTK_REAL const inv_alp = (alp_ijk > 0.0) ? (1.0 / alp_ijk) : 0.0;
          CCTK_REAL const F_u_x = (a + beta_u_x) * E_new * inv_alp;
          CCTK_REAL const F_u_y = beta_u_y * E_new * inv_alp;
          CCTK_REAL const F_u_z = beta_u_z * E_new * inv_alp;

          Fx_new = g_xx * F_u_x + g_xy * F_u_y + g_xz * F_u_z;
          Fy_new = g_xy * F_u_x + g_yy * F_u_y + g_yz * F_u_z;
          Fz_new = g_xz * F_u_x + g_yz * F_u_y + g_zz * F_u_z;
          N_new = 1.0;
        }

        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = E_new;
          rFx[i4D] = Fx_new;
          rFy[i4D] = Fy_new;
          rFz[i4D] = Fz_new;
          rN[i4D] = N_new;
        }
      });

  // +X: vacuum
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (p.BI[0] != 1.0)
          return;
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = 0.0;
          rFx[i4D] = 0.0;
          rFy[i4D] = 0.0;
          rFz[i4D] = 0.0;
          rN[i4D] = 0.0;
        }
      });

  // -Y: mimic THC exactly by copying every lower ghost layer from the first
  // interior Y plane (j = nghy), but only on the physical -Y boundary.
  const int jsrc_lo = nghy;
  for (int d = 0; d < nghy; ++d) {
    const int jdst = d;
    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.BI[1] != -1.0 || p.j != jdst)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4Db = layout_cc.linear(p.i, jdst, p.k, ig);
            const int i4Di = layout_cc.linear(p.i, jsrc_lo, p.k, ig);
            rE[i4Db] = rE[i4Di];
            rFx[i4Db] = rFx[i4Di];
            rFy[i4Db] = rFy[i4Di];
            rFz[i4Db] = rFz[i4Di];
            rN[i4Db] = rN[i4Di];
          }
        });
  }

  // +Y: mimic THC exactly by copying every upper ghost layer from the last
  // interior Y plane (j = cctk_lsh[1] - nghy - 1), only on physical +Y.
  const int jsrc_hi = cctk_lsh[1] - nghy - 1;
  for (int d = 0; d < nghy; ++d) {
    const int jdst = cctk_lsh[1] - nghy + d;
    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.BI[1] != 1.0 || p.j != jdst)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4Db = layout_cc.linear(p.i, jdst, p.k, ig);
            const int i4Di = layout_cc.linear(p.i, jsrc_hi, p.k, ig);
            rE[i4Db] = rE[i4Di];
            rFx[i4Db] = rFx[i4Di];
            rFy[i4Db] = rFy[i4Di];
            rFz[i4Db] = rFz[i4Di];
            rN[i4Db] = rN[i4Di];
          }
        });
  }

  // -Z: vacuum
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (p.BI[2] != -1.0)
          return;
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = 0.0;
          rFx[i4D] = 0.0;
          rFy[i4D] = 0.0;
          rFz[i4D] = 0.0;
          rN[i4D] = 0.0;
        }
      });

  // +Z: vacuum
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (p.BI[2] != 1.0)
          return;
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = 0.0;
          rFx[i4D] = 0.0;
          rFy[i4D] = 0.0;
          rFz[i4D] = 0.0;
          rN[i4D] = 0.0;
        }
      });
}

extern "C" void nuX_Seeds_KerrSchild_Mask(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_KerrSchild_Mask;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const CCTK_REAL mask_r2 = kerr_mask_radius * kerr_mask_radius;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        const bool masked = (p.x * p.x + p.y * p.y + p.z * p.z) < mask_r2;
        nuX_m1_mask[ijk] = masked ? 1.0 : 0.0;

        if (masked) {
          for (int ig = 0; ig < nspecies * ngroups; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rN[i4D] = 0.0;
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
          }
        }
      });
}

} // namespace nuX_Seeds
