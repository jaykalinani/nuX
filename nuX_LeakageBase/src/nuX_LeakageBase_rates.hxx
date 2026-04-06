#ifndef NUX_LEAKAGEBASE_RATES_HXX
#define NUX_LEAKAGEBASE_RATES_HXX

#include <cctk.h>

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

CCTK_HOST CCTK_DEVICE inline void setup_equilibrium_grey_opacity_params(
    GreyOpacityParams &grey_opacity_params, const OpacityFlags &opacity_flags,
    const OpacityParams &opacity_pars, EOSX::eos_3p_tabulated3d *const eos_3p,
    const CCTK_REAL rho, const CCTK_REAL temp, const CCTK_REAL ye,
    const CCTK_REAL particle_mass) {
  grey_opacity_params = {};
  grey_opacity_params.opacity_flags = opacity_flags;
  grey_opacity_params.opacity_pars = opacity_pars;

  grey_opacity_params.eos_pars.nb =
      rho * nuX_dens_conv / (particle_mass * kBS_MeVtog);
  grey_opacity_params.eos_pars.temp = temp;
  grey_opacity_params.eos_pars.ye = ye;
  grey_opacity_params.eos_pars.yp = ye;
  grey_opacity_params.eos_pars.yn = 1.0 - ye;

  eos_3p->mu_pne_from_valid_rho_temp_ye(
      rho, temp, ye, grey_opacity_params.eos_pars.mu_p,
      grey_opacity_params.eos_pars.mu_n, grey_opacity_params.eos_pars.mu_e);

  grey_opacity_params.distr_pars =
      NuEquilibriumParams(&grey_opacity_params.eos_pars);
  ComputeM1DensitiesEq(&grey_opacity_params.eos_pars,
                       &grey_opacity_params.distr_pars,
                       &grey_opacity_params.m1_pars);
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
  if (dens <= 0.0 || kappa <= 0.0) {
    return 0.0;
  }
  CCTK_REAL const itloss = r_free / dens;
  CCTK_REAL const lambda = 1.0 / kappa;
  CCTK_REAL const tdiff = diff_fact * lambda * tau * tau;
  return r_free / (1.0 + itloss * tdiff);
}

} // namespace nuX_LeakageBase

#endif
