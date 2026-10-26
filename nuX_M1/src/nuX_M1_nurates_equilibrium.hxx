#ifndef NUX_M1_NURATES_EQUILIBRIUM_HXX
#define NUX_M1_NURATES_EQUILIBRIUM_HXX

#include "m1_opacities.hpp"
#include "nuX_M1_nurates.hxx"
#include "nuX_M1_weak_equil.hxx"
#include "nuX_rate_units.hxx"

namespace nuX_M1 {

// Use the same equilibrium prescription in the opacity evaluation and the
// population-dependent source solve. Radiation moments are physical,
// undensitized, and per flavor, as required by the NuRates interface.
template <typename EOS>
CCTK_HOST CCTK_DEVICE inline void nurates_equilibrium_targets(
    const EOS *const eos, const CCTK_REAL rho, const MyEOSParams &matter,
    const M1Quantities &radiation, const nurates_raw_opacities_t *const raw,
    const CCTK_REAL particle_mass, const CCTK_REAL proper_dt,
    const CCTK_REAL tau_trap, const CCTK_REAL tau_delta,
    CCTK_REAL *const number, CCTK_REAL *const energy) {
  const CCTK_REAL tau =
      fmin(sqrt(raw[0].abs_1_thermal *
                (raw[0].abs_1_thermal + raw[0].abs_1_nonthermal + raw[0].scat_1)),
           sqrt(raw[1].abs_1_thermal *
                (raw[1].abs_1_thermal + raw[1].abs_1_nonthermal + raw[1].scat_1))) *
      proper_dt;

  CCTK_REAL thin_number[3], thin_energy[3];
  NeutrinoDens(matter.mu_n, matter.mu_p, matter.mu_e, matter.temp,
               thin_number[0], thin_number[1], thin_number[2], thin_energy[0],
               thin_energy[1], thin_energy[2]);
  for (int sp = 0; sp < 3; ++sp)
    fallback_equilibrium_moments(thin_number[sp], thin_energy[sp], 0, 0);

  CCTK_REAL trapped_number[3]{}, trapped_energy[3]{};
  if (tau_trap >= 0 && tau > tau_trap) {
    const CCTK_REAL eps =
        eos->eps_from_rho_temp_ye(rho, matter.temp, matter.ye);
    // Species 2 represents four evolved heavy-lepton flavors. Its NuRates
    // input is per flavor, whereas the conserved energy budget is their sum.
    const CCTK_REAL eps_total =
        eps + (radiation.J[0] + radiation.J[1] + 4 * radiation.J[2]) /
                  (nuX_edens_conv * rho);
    const CCTK_REAL lepton_fraction =
        matter.ye + (radiation.n[0] - radiation.n[1]) / matter.nb;
    const CCTK_REAL nb_fm3 =
        matter.nb / nuX_Utils::rate_units::physical_number_density_fm3_to_nm3;
    CCTK_REAL temperature = matter.temp, ye = matter.ye;
    int status = BetaEquilibriumTrapped(
        rho, nb_fm3, particle_mass, eps_total, lepton_fraction, temperature,
        ye, matter.temp, matter.ye, eos);
    if (status)
      status = BetaEquilibriumTrapped(
          rho, nb_fm3, particle_mass, eps, matter.ye, temperature, ye,
          matter.temp, matter.ye, eos);

    if (status == 0) {
      CCTK_REAL mu_p, mu_n, mu_e;
      eos->mu_pne_from_rho_temp_ye(rho, temperature, ye, mu_p, mu_n, mu_e);
      NeutrinoDens(mu_n, mu_p, mu_e, temperature, trapped_number[0],
                   trapped_number[1], trapped_number[2], trapped_energy[0],
                   trapped_energy[1], trapped_energy[2]);
    }
    for (int sp = 0; sp < 3; ++sp) {
      if (status != 0) {
        trapped_number[sp] = thin_number[sp];
        trapped_energy[sp] = thin_energy[sp];
      } else {
        fallback_equilibrium_moments(trapped_number[sp], trapped_energy[sp],
                                     thin_number[sp], thin_energy[sp]);
      }
    }
  }

  for (int sp = 0; sp < 3; ++sp) {
    if (tau_trap < 0 || tau <= tau_trap) {
      number[sp] = thin_number[sp];
      energy[sp] = thin_energy[sp];
    } else if (tau > tau_trap + tau_delta) {
      number[sp] = trapped_number[sp];
      energy[sp] = trapped_energy[sp];
    } else {
      const CCTK_REAL weight = (tau - tau_trap) / tau_delta;
      number[sp] = weight * trapped_number[sp] + (1 - weight) * thin_number[sp];
      energy[sp] = weight * trapped_energy[sp] + (1 - weight) * thin_energy[sp];
    }
  }
}

} // namespace nuX_M1

#endif // NUX_M1_NURATES_EQUILIBRIUM_HXX
