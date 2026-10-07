#include <algorithm>
#include <cassert>
#include <cmath>

#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M1_closure.hxx"
#include "nuX_M1_opacity_utils.hxx"
#include "nuX_rate_units.hxx"
#include "nuX_utils.hxx"
#include "nuX_weakrates.hxx"
#include "setup_eos.hxx"

namespace nuX_M1 {

using namespace std;
using namespace Loop;
using namespace nuX_Utils;
using namespace EOSX;

extern "C" void nuX_M1_SetToEquilibriumWeakRates(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_SetToEquilibriumWeakRates;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_SetToEquilibriumWeakRates");

  assert(nspecies == 3);
  assert(ngroups == 1);

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});

  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);
  tensor::fluid_velocity_field_const fidu(layout_vc, layout_cc, alp, betax,
                                          betay, betaz, fidu_w_lorentz,
                                          fidu_velx, fidu_vely, fidu_velz);

  auto eos_3p = global_eos_3p_tab3d;
  auto weakrates = nuX_WeakRates::global_weakrates;
  if (!eos_3p)
    CCTK_ERROR("nuX_M1_SetToEquilibriumWeakRates requires a tabulated EOS");
  if (!weakrates)
    CCTK_ERROR("nuX_M1_SetToEquilibriumWeakRates requires nuX_WeakRates");

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        const CCTK_REAL rho_raw = rho[ijk];
        const CCTK_REAL temp_raw = temperature[ijk];
        const CCTK_REAL ye_raw = Ye[ijk];
        if (nuX_m1_mask[ijk] || !isfinite(rho_raw) || rho_raw <= 0.0 ||
            !isfinite(temp_raw) || temp_raw <= 0.0 || !isfinite(ye_raw) ||
            ye_raw < 0.0 || ye_raw > 1.0)
          return;

        const CCTK_REAL rho_cgs = rho_raw * rate_units::code_density_to_g_nm3 /
                                  rate_units::per_cm3_to_per_nm3;
        if (rho_cgs < equilibrium_rho_min)
          return;
        const CCTK_REAL rho_eos =
            fmin(fmax(rho_raw, eos_3p->rgrho.min), eos_3p->rgrho.max);
        const CCTK_REAL temp_eos =
            fmin(fmax(temp_raw, eos_3p->rgtemp.min), eos_3p->rgtemp.max);
        const CCTK_REAL ye_eos =
            fmin(fmax(ye_raw, eos_3p->rgye.min), eos_3p->rgye.max);

        CCTK_REAL mu_p, mu_n, mu_e;
        eos_3p->mu_pne_from_rho_temp_ye(rho_eos, temp_eos, ye_eos, mu_p, mu_n,
                                        mu_e);

        const nuX_WeakRates::EOSState weak_eos = {
            rho_eos * rate_units::code_density_to_g_nm3 /
                rate_units::per_cm3_to_per_nm3,
            temp_eos,
            ye_eos,
            particle_mass,
            mu_e,
            mu_p,
            mu_n,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0};
        const auto densities = weakrates->equilibrium_densities(weak_eos);

        CCTK_REAL nudens_0[3], nudens_1[3];
        for (int ig = 0; ig < nspecies * ngroups; ++ig) {
          nudens_0[ig] = densities.number[ig] * rate_units::per_cm3_to_per_nm3 /
                         rate_units::transport_number_density_to_nm3;
          nudens_1[ig] = densities.energy[ig] * rate_units::per_cm3_to_per_nm3 /
                         rate_units::code_energy_density_to_mev_nm3;
          fallback_equilibrium_moments(nudens_0[ig], nudens_1[ig],
                                       CCTK_REAL(0), CCTK_REAL(0));
        }

        tensor::metric<4> g_dd;
        tensor::inv_metric<4> g_uu;
        tensor::generic<CCTK_REAL, 4, 1> n_u, n_d;
        tensor::generic<CCTK_REAL, 4, 2> gamma_ud;
        geom.get_metric(p, &g_dd);
        geom.get_inv_metric(p, &g_uu);
        geom.get_normal(p, &n_u);
        geom.get_normal_form(p, &n_d);
        geom.get_space_proj(p, &gamma_ud);

        tensor::generic<CCTK_REAL, 4, 1> u_u, u_d;
        fidu.get(p, &u_u);
        tensor::contract(g_dd, u_u, &u_d);

        const CCTK_REAL detg =
            nuX_Utils::metric::spatial_det(g_dd(1, 1), g_dd(1, 2), g_dd(1, 3),
                                           g_dd(2, 2), g_dd(2, 3), g_dd(3, 3));
        if (!isfinite(detg) || detg <= 0.0 || !isfinite(fidu_w_lorentz[ijk]) ||
            fidu_w_lorentz[ijk] < 1.0)
          return;
        const CCTK_REAL volform = sqrt(detg);

        for (int ig = 0; ig < nspecies * ngroups; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          const CCTK_REAL J = volform * nudens_1[ig];

          tensor::generic<CCTK_REAL, 4, 1> H_d;
          tensor::symmetric2<CCTK_REAL, 4, 2> K_dd, rT_dd;
          for (int a = 0; a < 4; ++a) {
            H_d(a) = 0.0;
            for (int b = a; b < 4; ++b)
              K_dd(a, b) = J * (g_dd(a, b) + u_d(a) * u_d(b)) / 3.0;
          }
          assemble_rT(u_d, J, H_d, K_dd, &rT_dd);

          CCTK_REAL E = calc_J_from_rT(rT_dd, n_u);
          tensor::generic<CCTK_REAL, 4, 1> F_d;
          calc_H_from_rT(rT_dd, n_u, gamma_ud, &F_d);
          apply_floor(g_uu, &E, &F_d, rad_E_floor, rad_eps);

          rE[i4D] = E;
          unpack_F_d(F_d, &rFx[i4D], &rFy[i4D], &rFz[i4D]);
          rN[i4D] =
              max(volform * nudens_0[ig] * fidu_w_lorentz[ijk], rad_N_floor);

          assert(isfinite(rN[i4D]));
          assert(isfinite(rE[i4D]));
          assert(isfinite(rFx[i4D]));
          assert(isfinite(rFy[i4D]));
          assert(isfinite(rFz[i4D]));
        }
      });
}

} // namespace nuX_M1
