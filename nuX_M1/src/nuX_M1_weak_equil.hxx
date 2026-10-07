#ifndef NUX_M1_WEAK_EQUIL_HXX
#define NUX_M1_WEAK_EQUIL_HXX

#include <cctk.h>

#include <limits>

#include "setup_eos.hxx"

namespace nuX_M1 {

// The equilibrium equations are written in EOS/code units:
//
//   Y_l       = Y_e + (n_nue - n_anue) / n_b,
//   eps_total = eps_matter + e_nu / (n_b m_b).
//
// Here n_b is in fm^-3, m_b in MeV, eps is a specific energy in units of
// c^2, and e_nu is in MeV fm^-3.  Keeping the residual in these units avoids
// mixing an EOS specific energy with a radiation energy density.
inline constexpr CCTK_REAL nu_2DNR_residual_abs_tol = 1.0e-12;
inline constexpr CCTK_REAL nu_2DNR_residual_rel_tol = 1.0e-9;
// weak_equilibrium_residual_norm returns residuals normalized by the combined
// absolute and relative tolerance, so convergence corresponds to norm <= 1.
inline constexpr CCTK_REAL nu_2DNR_residual_tol = 1.0;
inline constexpr CCTK_REAL nu_2DNR_step_tol = 1.0e-10;
inline constexpr int nu_2DNR_n_max = 100;
inline constexpr int nu_bis_n_cut_max = 12;

inline constexpr CCTK_REAL hc_mevfm = 1.23984172e3;
inline constexpr CCTK_REAL pi = 3.14159265358979323846;
inline constexpr CCTK_REAL pi2 = pi * pi;
inline constexpr CCTK_REAL pi4 = pi2 * pi2;

inline constexpr CCTK_REAL nu_n_prefactor =
    4.0 / 3.0 * pi / (hc_mevfm * hc_mevfm * hc_mevfm);
inline constexpr CCTK_REAL nu_e_prefactor =
    4.0 * pi / (hc_mevfm * hc_mevfm * hc_mevfm);

inline constexpr CCTK_REAL nu_7pi4_60 = 7.0 * pi4 / 60.0;
inline constexpr CCTK_REAL nu_7pi4_30 = 7.0 * pi4 / 30.0;
inline constexpr CCTK_REAL nu_7pi4_15 = 7.0 * pi4 / 15.0;
inline constexpr CCTK_REAL nu_14pi4_15 = 14.0 * pi4 / 15.0;

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
weak_equilibrium_residual_norm(const CCTK_REAL Yl, const CCTK_REAL eps_total,
                               const CCTK_REAL y[2]) {
  const CCTK_REAL y_scale =
      nu_2DNR_residual_abs_tol + nu_2DNR_residual_rel_tol * abs(Yl);
  const CCTK_REAL e_scale =
      nu_2DNR_residual_abs_tol + nu_2DNR_residual_rel_tol * abs(eps_total);
  return fmax(abs(y[0]) / y_scale, abs(y[1]) / e_scale);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
weak_equilibrium_finite2(const CCTK_REAL x[2]) {
  return isfinite(x[0]) && isfinite(x[1]);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
inv_jacobi(const CCTK_REAL det, const CCTK_REAL J[2][2], CCTK_REAL invJ[2][2]) {
  const CCTK_REAL inv_det = 1.0 / det;
  invJ[0][0] = J[1][1] * inv_det;
  invJ[1][1] = J[0][0] * inv_det;
  invJ[0][1] = -J[0][1] * inv_det;
  invJ[1][0] = -J[1][0] * inv_det;
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
eta_e_gradient(const CCTK_REAL rho, const CCTK_REAL T, const CCTK_REAL Ye,
               const CCTK_REAL eta, CCTK_REAL &deta_dT, CCTK_REAL &deta_dYe,
               CCTK_REAL &deps_dT, CCTK_REAL &deps_dYe,
               const EOSType *const tabeos) {
  const CCTK_REAL min_T = tabeos->rgtemp.min;
  const CCTK_REAL max_T = tabeos->rgtemp.max;
  const CCTK_REAL min_Y = tabeos->rgye.min;
  const CCTK_REAL max_Y = tabeos->rgye.max;

  if (!isfinite(rho) || !isfinite(T) || !isfinite(Ye) || !isfinite(eta) ||
      rho <= 0.0 || T <= 0.0 || !(max_T > min_T) || !(max_Y > min_Y)) {
    return 1;
  }

  // Use bound-aware centered differences.  The steps are large enough to be
  // robust for tabulated interpolation while remaining local on production
  // tables.
  const CCTK_REAL dY =
      fmax(CCTK_REAL(1.0e-6), CCTK_REAL(1.0e-4) * (max_Y - min_Y));
  const CCTK_REAL Y1 = fmax(Ye - dY, min_Y);
  const CCTK_REAL Y2 = fmin(Ye + dY, max_Y);
  if (!(Y2 > Y1)) {
    return 1;
  }

  const CCTK_REAL mu_Y1 = tabeos->mu_lepton_from_rho_temp_ye(rho, T, Y1);
  const CCTK_REAL mu_Y2 = tabeos->mu_lepton_from_rho_temp_ye(rho, T, Y2);
  const CCTK_REAL eps_Y1 = tabeos->eps_from_rho_temp_ye(rho, T, Y1);
  const CCTK_REAL eps_Y2 = tabeos->eps_from_rho_temp_ye(rho, T, Y2);
  const CCTK_REAL dmu_dYe = (mu_Y2 - mu_Y1) / (Y2 - Y1);
  deps_dYe = (eps_Y2 - eps_Y1) / (Y2 - Y1);

  const CCTK_REAL dT =
      fmax(CCTK_REAL(1.0e-6), CCTK_REAL(1.0e-4) * fmax(T, CCTK_REAL(1.0)));
  const CCTK_REAL T1 = fmax(T - dT, min_T);
  const CCTK_REAL T2 = fmin(T + dT, max_T);
  if (!(T2 > T1)) {
    return 1;
  }

  // Keep Ye fixed for the temperature derivative.  The previous
  // implementation accidentally evaluated these at Y1 and Y2.
  const CCTK_REAL mu_T1 = tabeos->mu_lepton_from_rho_temp_ye(rho, T1, Ye);
  const CCTK_REAL mu_T2 = tabeos->mu_lepton_from_rho_temp_ye(rho, T2, Ye);
  const CCTK_REAL eps_T1 = tabeos->eps_from_rho_temp_ye(rho, T1, Ye);
  const CCTK_REAL eps_T2 = tabeos->eps_from_rho_temp_ye(rho, T2, Ye);
  const CCTK_REAL dmu_dT = (mu_T2 - mu_T1) / (T2 - T1);
  deps_dT = (eps_T2 - eps_T1) / (T2 - T1);

  deta_dT = (dmu_dT - eta) / T;
  deta_dYe = dmu_dYe / T;

  return isfinite(deta_dT) && isfinite(deta_dYe) && isfinite(deps_dT) &&
                 isfinite(deps_dYe)
             ? 0
             : 1;
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
jacobi_eq_weak(const CCTK_REAL rho, const CCTK_REAL n_b,
               const CCTK_REAL particle_mass, const CCTK_REAL x[2],
               CCTK_REAL J[2][2], const EOSType *const tabeos) {
  const CCTK_REAL T = x[0];
  const CCTK_REAL Ye = x[1];
  if (!isfinite(rho) || rho <= 0.0 || !weak_equilibrium_finite2(x) ||
      T <= 0.0 || n_b <= 0.0 || particle_mass <= 0.0) {
    return 1;
  }

  const CCTK_REAL mu_l = tabeos->mu_lepton_from_rho_temp_ye(rho, T, Ye);
  const CCTK_REAL eta = mu_l / T;
  if (!isfinite(eta)) {
    return 1;
  }

  CCTK_REAL deta_dT, deta_dYe, deps_dT, deps_dYe;
  if (eta_e_gradient(rho, T, Ye, eta, deta_dT, deta_dYe, deps_dT, deps_dYe,
                     tabeos) != 0) {
    return 1;
  }

  const CCTK_REAL eta2 = eta * eta;
  const CCTK_REAL T2 = T * T;
  const CCTK_REAL T3 = T2 * T;
  const CCTK_REAL inv_nb = 1.0 / n_b;
  const CCTK_REAL inv_nb_mb = inv_nb / particle_mass;

  J[0][0] = nu_n_prefactor * inv_nb * T2 *
            (3.0 * eta * (pi2 + eta2) + T * (pi2 + 3.0 * eta2) * deta_dT);
  J[0][1] = 1.0 + nu_n_prefactor * inv_nb * T3 * (pi2 + 3.0 * eta2) * deta_dYe;

  J[1][0] = deps_dT +
            nu_e_prefactor * inv_nb_mb * T3 *
                (nu_7pi4_15 + nu_14pi4_15 + 2.0 * eta2 * (pi2 + 0.5 * eta2) +
                 eta * T * (pi2 + eta2) * deta_dT);
  J[1][1] = deps_dYe +
            nu_e_prefactor * inv_nb_mb * T3 * T * eta * (pi2 + eta2) * deta_dYe;

  return isfinite(J[0][0]) && isfinite(J[0][1]) && isfinite(J[1][0]) &&
                 isfinite(J[1][1])
             ? 0
             : 1;
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
func_eq_weak(const CCTK_REAL rho, const CCTK_REAL n_b,
             const CCTK_REAL particle_mass, const CCTK_REAL eps_total,
             const CCTK_REAL Yl, const CCTK_REAL x[2], CCTK_REAL y[2],
             const EOSType *const tabeos) {
  const CCTK_REAL T = x[0];
  const CCTK_REAL Ye = x[1];
  if (!isfinite(rho) || rho <= 0.0 || !weak_equilibrium_finite2(x) ||
      T <= 0.0 || n_b <= 0.0 || particle_mass <= 0.0 ||
      !isfinite(eps_total) || !isfinite(Yl)) {
    return 1;
  }

  const CCTK_REAL mu_l = tabeos->mu_lepton_from_rho_temp_ye(rho, T, Ye);
  const CCTK_REAL eps = tabeos->eps_from_rho_temp_ye(rho, T, Ye);
  const CCTK_REAL eta = mu_l / T;
  if (!isfinite(eta) || !isfinite(eps)) {
    return 1;
  }

  const CCTK_REAL eta2 = eta * eta;
  const CCTK_REAL T3 = T * T * T;
  const CCTK_REAL T4 = T3 * T;
  const CCTK_REAL net_nue = nu_n_prefactor * T3 * eta * (pi2 + eta2);
  const CCTK_REAL neutrino_energy =
      nu_e_prefactor * T4 *
      (nu_7pi4_60 + 0.5 * eta2 * (pi2 + 0.5 * eta2) + nu_7pi4_30);

  y[0] = Ye + net_nue / n_b - Yl;
  y[1] = eps + neutrino_energy / (n_b * particle_mass) - eps_total;
  return weak_equilibrium_finite2(y) ? 0 : 1;
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE inline int trapped_equilibrium_2DNR(
    const CCTK_REAL rho, const CCTK_REAL n_b, const CCTK_REAL particle_mass,
    const CCTK_REAL eps_total, const CCTK_REAL Yl, const CCTK_REAL x0[2],
    CCTK_REAL x1[2], const EOSType *const tabeos) {
  const CCTK_REAL min_T = tabeos->rgtemp.min;
  const CCTK_REAL max_T = tabeos->rgtemp.max;
  const CCTK_REAL min_Y = tabeos->rgye.min;
  const CCTK_REAL max_Y = tabeos->rgye.max;
  if (!isfinite(rho) || !isfinite(n_b) || !isfinite(particle_mass) ||
      !isfinite(eps_total) || !isfinite(Yl) || rho <= 0.0 || n_b <= 0.0 ||
      particle_mass <= 0.0 || !(max_T > min_T) || !(max_Y > min_Y)) {
    return 1;
  }

  x1[0] = fmin(fmax(x0[0], min_T), max_T);
  x1[1] = fmin(fmax(x0[1], min_Y), max_Y);

  CCTK_REAL y[2];
  if (func_eq_weak(rho, n_b, particle_mass, eps_total, Yl, x1, y, tabeos) !=
      0) {
    return 1;
  }
  CCTK_REAL err = weak_equilibrium_residual_norm(Yl, eps_total, y);
  if (!isfinite(err)) {
    return 1;
  }
  const CCTK_REAL eps = std::numeric_limits<CCTK_REAL>::epsilon();
  // Include the state produced by the final permitted Newton update in the
  // convergence checks.  No additional update is taken at that boundary.
  for (int iter = 0; iter <= nu_2DNR_n_max; ++iter) {
    CCTK_REAL J[2][2];
    if (jacobi_eq_weak(rho, n_b, particle_mass, x1, J, tabeos) != 0) {
      return 1;
    }

    const CCTK_REAL det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    const CCTK_REAL det_scale = abs(J[0][0] * J[1][1]) + abs(J[0][1] * J[1][0]);
    if (!isfinite(det) || !isfinite(det_scale) || !(det_scale > 0.0) ||
        abs(det) <= CCTK_REAL(64.0) * eps * det_scale) {
      return 1;
    }

    CCTK_REAL invJ[2][2];
    inv_jacobi(det, J, invJ);

    // Assess conditioning in the same dimensionless variables used by the
    // convergence tests.  The raw Jacobian mixes temperature and electron
    // fraction with lepton and specific-energy residuals, so its condition
    // number would otherwise depend on the chosen physical units.
    const CCTK_REAL residual_scale[2] = {
        nu_2DNR_residual_abs_tol + nu_2DNR_residual_rel_tol * abs(Yl),
        nu_2DNR_residual_abs_tol + nu_2DNR_residual_rel_tol * abs(eps_total)};
    const CCTK_REAL state_scale[2] = {fmax(abs(x1[0]), CCTK_REAL(1.0)),
                                      fmax(abs(x1[1]), CCTK_REAL(1.0))};
    const CCTK_REAL scaled_J[2][2] = {
        {J[0][0] * state_scale[0] / residual_scale[0],
         J[0][1] * state_scale[1] / residual_scale[0]},
        {J[1][0] * state_scale[0] / residual_scale[1],
         J[1][1] * state_scale[1] / residual_scale[1]}};
    const CCTK_REAL scaled_det =
        scaled_J[0][0] * scaled_J[1][1] - scaled_J[0][1] * scaled_J[1][0];
    const CCTK_REAL norm_J = fmax(abs(scaled_J[0][0]) + abs(scaled_J[0][1]),
                                  abs(scaled_J[1][0]) + abs(scaled_J[1][1]));
    const CCTK_REAL norm_invJ =
        fmax((abs(scaled_J[1][1]) + abs(scaled_J[0][1])) / abs(scaled_det),
             (abs(scaled_J[1][0]) + abs(scaled_J[0][0])) / abs(scaled_det));
    const CCTK_REAL condition = norm_J * norm_invJ;
    if (!isfinite(condition) ||
        condition > CCTK_REAL(1.0) / (CCTK_REAL(64.0) * eps)) {
      return 1;
    }
    CCTK_REAL dx[2] = {-(invJ[0][0] * y[0] + invJ[0][1] * y[1]),
                       -(invJ[1][0] * y[0] + invJ[1][1] * y[1])};
    if (!weak_equilibrium_finite2(dx)) {
      return 1;
    }

    // Remove only components that point out of the EOS domain.  Inward
    // components at a bound must remain active.
    if ((x1[0] <= min_T && dx[0] < 0.0) || (x1[0] >= max_T && dx[0] > 0.0)) {
      dx[0] = 0.0;
    }
    if ((x1[1] <= min_Y && dx[1] < 0.0) || (x1[1] >= max_Y && dx[1] > 0.0)) {
      dx[1] = 0.0;
    }

    // A small residual alone is not sufficient: it can occur at a poorly
    // conditioned point or after projection onto an EOS boundary.  Require
    // the corresponding projected Newton correction to be small as well.
    // This check also handles an exactly converged initial guess without
    // forcing a zero-length line-search step.
    const CCTK_REAL newton_step =
        fmax(abs(dx[0]) / fmax(abs(x1[0]), CCTK_REAL(1.0)),
             abs(dx[1]) / fmax(abs(x1[1]), CCTK_REAL(1.0)));
    if (err <= nu_2DNR_residual_tol && newton_step <= nu_2DNR_step_tol) {
      return 0;
    }
    if (iter == nu_2DNR_n_max) {
      return 1;
    }

    const CCTK_REAL base[2] = {x1[0], x1[1]};
    const CCTK_REAL err_old = err;
    bool accepted = false;
    CCTK_REAL trial_y[2];
    CCTK_REAL trial[2];
    CCTK_REAL fac = 1.0;
    for (int cut = 0; cut <= nu_bis_n_cut_max; ++cut, fac *= 0.5) {
      trial[0] = fmin(fmax(base[0] + fac * dx[0], min_T), max_T);
      trial[1] = fmin(fmax(base[1] + fac * dx[1], min_Y), max_Y);
      const CCTK_REAL step =
          fmax(abs(trial[0] - base[0]) / fmax(abs(base[0]), CCTK_REAL(1.0)),
               abs(trial[1] - base[1]) / fmax(abs(base[1]), CCTK_REAL(1.0)));
      if (step == 0.0 || func_eq_weak(rho, n_b, particle_mass, eps_total, Yl,
                                      trial, trial_y, tabeos) != 0) {
        continue;
      }
      const CCTK_REAL trial_err =
          weak_equilibrium_residual_norm(Yl, eps_total, trial_y);
      if (isfinite(trial_err) && trial_err < err_old) {
        x1[0] = trial[0];
        x1[1] = trial[1];
        y[0] = trial_y[0];
        y[1] = trial_y[1];
        err = trial_err;
        accepted = true;
        break;
      }
    }

    if (!accepted) {
      return 2;
    }
    // A small damped line-search step is not evidence that the Newton
    // correction is small.  Recompute the Jacobian and full projected Newton
    // step at the accepted state on the next iteration before declaring
    // convergence.
  }

  return 1;
}

/// Calculate hot, neutrino-trapped beta equilibrium at fixed density,
/// specific total internal energy, and total electron-lepton fraction.
template <typename EOSType>
CCTK_HOST CCTK_DEVICE inline int
BetaEquilibriumTrapped(const CCTK_REAL rho, const CCTK_REAL n_b,
                       const CCTK_REAL particle_mass, const CCTK_REAL eps_total,
                       const CCTK_REAL Yl, CCTK_REAL &T_eq, CCTK_REAL &Y_eq,
                       const CCTK_REAL T_guess, const CCTK_REAL Y_guess,
                       const EOSType *const tabeos) {
  constexpr int n_attempts = 16;
  constexpr CCTK_REAL guess_factors[n_attempts][2] = {
      {1.00, 1.00}, {0.90, 1.25}, {0.90, 1.10}, {0.90, 1.00},
      {0.90, 0.90}, {0.90, 0.75}, {0.75, 1.25}, {0.75, 1.10},
      {0.75, 1.00}, {0.75, 0.90}, {0.75, 0.75}, {0.50, 1.25},
      {0.50, 1.10}, {0.50, 1.00}, {0.50, 0.90}, {0.50, 0.75}};

  T_eq = T_guess;
  Y_eq = Y_guess;
  for (int attempt = 0; attempt < n_attempts; ++attempt) {
    const CCTK_REAL x0[2] = {guess_factors[attempt][0] * T_guess,
                             guess_factors[attempt][1] * Y_guess};
    CCTK_REAL x1[2];
    const int ierr = trapped_equilibrium_2DNR(rho, n_b, particle_mass,
                                              eps_total, Yl, x0, x1, tabeos);
    if (ierr == 0 && weak_equilibrium_finite2(x1)) {
      CCTK_REAL y[2];
      if (func_eq_weak(rho, n_b, particle_mass, eps_total, Yl, x1, y, tabeos) ==
              0 &&
          weak_equilibrium_residual_norm(Yl, eps_total, y) <=
              nu_2DNR_residual_tol) {
        T_eq = x1[0];
        Y_eq = x1[1];
        return 0;
      }
    }
  }
  return 1;
}

} // namespace nuX_M1

#endif // NUX_M1_WEAK_EQUIL_HXX
