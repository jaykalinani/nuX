#ifndef NUX_LEAKAGEBASE_RATES_HXX
#define NUX_LEAKAGEBASE_RATES_HXX

#include <cctk.h>
#include <cmath>
#include <limits>

#include "m1_opacities.hpp"
#include "setup_eos.hxx"

namespace nuX_LeakageBase {

CCTK_HOST CCTK_DEVICE inline void copy_quadrature(MyQuadrature &dst,
                                                  const MyQuadrature &src) {
  dst.type = src.type;
  dst.alpha = src.alpha;
  dst.dim = src.dim;
  dst.nx = src.nx;
  dst.ny = src.ny;
  dst.nz = src.nz;
  dst.x1 = src.x1;
  dst.x2 = src.x2;
  dst.y1 = src.y1;
  dst.y2 = src.y2;
  dst.z1 = src.z1;
  dst.z2 = src.z2;
  for (int idx = 0; idx < BS_N_MAX; ++idx) {
    dst.points[idx] = src.points[idx];
    dst.w[idx] = src.w[idx];
  }
}

CCTK_HOST CCTK_DEVICE inline bool setup_equilibrium_grey_opacity_params(
    GreyOpacityParams &grey_opacity_params, const OpacityFlags &opacity_flags,
    const OpacityParams &opacity_pars, EOSX::eos_3p_tabulated3d *const eos_3p,
    const CCTK_REAL rho, const CCTK_REAL temp, const CCTK_REAL ye,
    const CCTK_REAL particle_mass) {
  grey_opacity_params = {};
  grey_opacity_params.opacity_flags = opacity_flags;
  grey_opacity_params.opacity_pars = opacity_pars;

  if (!std::isfinite(rho) || rho <= CCTK_REAL(0) || !std::isfinite(temp) ||
      temp <= CCTK_REAL(0) || !std::isfinite(ye) || ye < CCTK_REAL(0) ||
      ye > CCTK_REAL(1) || !std::isfinite(particle_mass) ||
      particle_mass <= CCTK_REAL(0))
    return false;

  const CCTK_REAL rho_eos =
      std::fmin(std::fmax(rho, eos_3p->rgrho.min), eos_3p->rgrho.max);
  const CCTK_REAL temp_eos =
      std::fmin(std::fmax(temp, eos_3p->rgtemp.min), eos_3p->rgtemp.max);
  const CCTK_REAL ye_eos =
      std::fmin(std::fmax(ye, eos_3p->rgye.min), eos_3p->rgye.max);

  grey_opacity_params.eos_pars.nb =
      rho_eos * nuX_dens_conv / (particle_mass * kBS_MeVtog);
  grey_opacity_params.eos_pars.temp = temp_eos;
  grey_opacity_params.eos_pars.ye = ye_eos;

  // NuRates' yn/yp fields are free-nucleon abundances.  Charge neutrality
  // fixes the total proton fraction to Ye, but that is not the free proton
  // fraction when the EOS contains nuclei.
  using tabulated_eos = EOSX::eos_3p_tabulated3d;
  const auto eos_state = eos_3p->interptable->interpolate<
      tabulated_eos::EV::MU_P, tabulated_eos::EV::MU_N,
      tabulated_eos::EV::MU_E, tabulated_eos::EV::XN,
      tabulated_eos::EV::XP>(std::log(rho_eos), std::log(temp_eos), ye_eos);
  if (!std::isfinite(eos_state[0]) || !std::isfinite(eos_state[1]) ||
      !std::isfinite(eos_state[2]))
    return false;
  grey_opacity_params.eos_pars.mu_p = eos_state[0];
  grey_opacity_params.eos_pars.mu_n = eos_state[1];
  grey_opacity_params.eos_pars.mu_e = eos_state[2];
  grey_opacity_params.eos_pars.yn =
      std::isfinite(eos_state[3])
          ? std::fmin(std::fmax(eos_state[3], CCTK_REAL(0)), CCTK_REAL(1))
          : CCTK_REAL(1) - ye_eos;
  grey_opacity_params.eos_pars.yp =
      std::isfinite(eos_state[4])
          ? std::fmin(std::fmax(eos_state[4], CCTK_REAL(0)), CCTK_REAL(1))
          : ye_eos;

  grey_opacity_params.distr_pars =
      NuEquilibriumParams(&grey_opacity_params.eos_pars);
  ComputeM1DensitiesEq(&grey_opacity_params.eos_pars,
                       &grey_opacity_params.distr_pars,
                       &grey_opacity_params.m1_pars);
  return true;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL
convert_number_rate(const CCTK_REAL eta_0, const CCTK_REAL species_factor) {
  return species_factor * eta_0 / nuX_ndens_conv * nuX_time_conv;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL
convert_energy_rate(const CCTK_REAL eta_1, const CCTK_REAL species_factor) {
  return species_factor * eta_1 / nuX_edens_conv * nuX_time_conv;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL
convert_number_density(const CCTK_REAL ndens_nr,
                       const CCTK_REAL species_factor) {
  return species_factor * ndens_nr / nuX_ndens_conv;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL
convert_energy_density(const CCTK_REAL edens_nr,
                       const CCTK_REAL species_factor) {
  return species_factor * edens_nr / nuX_edens_conv;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL
calc_eff_rate(CCTK_REAL const r_free, CCTK_REAL const dens,
              CCTK_REAL const kappa, CCTK_REAL const tau,
              CCTK_REAL const diff_fact) {
  const CCTK_REAL eps = std::numeric_limits<CCTK_REAL>::epsilon();
  if (!std::isfinite(r_free) || !std::isfinite(dens) ||
      !std::isfinite(kappa) || !std::isfinite(tau) ||
      !std::isfinite(diff_fact) || r_free <= 0.0 ||
      dens <= eps * r_free || kappa < eps || tau < 0.0 ||
      diff_fact <= 0.0) {
    return 0.0;
  }
  CCTK_REAL const itloss = r_free / dens;
  CCTK_REAL const lambda = 1.0 / kappa;
  CCTK_REAL const tdiff = diff_fact * lambda * tau * tau;
  const CCTK_REAL denom = 1.0 + itloss * tdiff;
  return std::isfinite(denom) && denom > 0.0 ? r_free / denom : 0.0;
}

} // namespace nuX_LeakageBase

#endif
