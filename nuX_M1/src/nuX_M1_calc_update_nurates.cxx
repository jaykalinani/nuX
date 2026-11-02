#include <loop_device.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#define NUX_M1_CLOSURE_IMPLEMENTATION
#define NUX_M1_SOURCES_IMPLEMENTATION
#include "nuX_M1_nurates_source.hxx"
#include "nuX_M1_nurates_cell.hxx"

namespace nuX_M1 {

using namespace Loop;
using namespace nuX_Utils;

namespace {

// Only reached after a failed solve. Keep enough round-trip-precision input
// to replay that cell without adding diagnostic storage to successful runs.
// Repeat physical coordinates and the source interval on every line: local
// integer indices alone need not distinguish reports from separate AMR boxes.
[[noreturn]] CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_NOINLINE void
report_nurates_source_failure(
    const SourceUpdateContext &g, const MyEOSParams &matter,
    const nurates_source_options_t &options,
    const nurates_source_state_t base[3], const CCTK_REAL x,
    const CCTK_REAL y, const CCTK_REAL z, const CCTK_REAL rho,
    const CCTK_REAL volume, const CCTK_REAL macro_proper_dt) {
  printf("nuX_M1: coupled NuRates source failed at (%.17e,%.17e,%.17e) "
         "cell=(%d,%d,%d): dt=%.17e macro_proper_dt=%.17e "
         "alpha=%.17e W=%.17e volume=%.17e rho=%.17e T=%.17e Ye=%.17e\n",
         x, y, z, g.i, g.j, g.k, g.cdt, macro_proper_dt, g.alp, g.W,
         volume, rho, matter.temp, matter.ye);
  printf("nuX_M1: source-context (%d,%d,%d) at=(%.17e,%.17e,%.17e) "
         "dt=%.17e matter: nb=%.17e "
         "mu_p=%.17e mu_n=%.17e mu_e=%.17e xn=%.17e xp=%.17e\n",
         g.i, g.j, g.k, x, y, z, g.cdt,
         matter.nb, matter.mu_p, matter.mu_n, matter.mu_e,
         matter.yn, matter.yp);
  printf("nuX_M1: source-context (%d,%d,%d) at=(%.17e,%.17e,%.17e) "
         "dt=%.17e solver: closure=%d "
         "closure_eps=%.17e closure_maxiter=%d epsrel=%.17e "
         "maxiter=%d inner_maxiter=%d angular_maxiter=%d\n",
         g.i, g.j, g.k, x, y, z, g.cdt,
         int(options.closure), g.closure_epsilon,
         int(g.closure_maxiter), options.epsrel, int(options.maxiter),
         int(options.inner_maxiter), int(options.angular_maxiter));
  for (int a = 0; a < 4; ++a) {
    printf("nuX_M1: source-context (%d,%d,%d) at=(%.17e,%.17e,%.17e) "
           "dt=%.17e vector[%d]: "
           "nd=%.17e nu=%.17e ud=%.17e uu=%.17e vd=%.17e vu=%.17e\n",
           g.i, g.j, g.k, x, y, z, g.cdt,
           a, g.n_d(a), g.n_u(a), g.u_d(a), g.u_u(a),
           g.v_d(a), g.v_u(a));
    for (int b = 0; b < 4; ++b)
      printf("nuX_M1: source-context (%d,%d,%d) at=(%.17e,%.17e,%.17e) "
             "dt=%.17e tensor[%d,%d]: "
             "gdd=%.17e guu=%.17e gamma=%.17e proj=%.17e\n",
             g.i, g.j, g.k, x, y, z, g.cdt,
             a, b, g.g_dd(a, b), g.g_uu(a, b),
             g.gamma_ud(a, b), g.proj_ud(a, b));
  }
  for (int sp = 0; sp < 3; ++sp) {
    const auto &s = base[sp];
    printf("nuX_M1: source-context (%d,%d,%d) at=(%.17e,%.17e,%.17e) "
           "dt=%.17e species[%d]: "
           "N=%.17e E=%.17e Ft=%.17e Fx=%.17e Fy=%.17e Fz=%.17e chi=%.17e\n",
           g.i, g.j, g.k, x, y, z, g.cdt,
           sp, s.N, s.E, s.F_d(0), s.F_d(1),
           s.F_d(2), s.F_d(3), s.chi);
  }
  source_solver_abort();
}

} // namespace

void launch_nurates_source_increment_kernel(
    const cGH *const cctkGH, const GridDescBaseDevice &grid,
    const GF3D2layout layout_cc, const tensor::slicing_geometry_const geom,
    const tensor::fluid_velocity_field_const fidu,
    const closure_t closure_fun, const CCTK_REAL dt,
    const bool store_accepted_state) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcUpdate;
  DECLARE_CCTK_PARAMETERS;

  if (ngroups != 1 || nspecies != 3)
    CCTK_ERROR("The coupled NuRates source requires one group and three species");
  const auto eos = EOSX::global_eos_3p_tab3d;
  if (!eos)
    CCTK_ERROR("The coupled NuRates source requires a tabulated EOS");
  const CCTK_REAL macro_dt_saved = ODESolvers_GetStepDeltaTime();
  const CCTK_REAL macro_dt = macro_dt_saved > 0 ? macro_dt_saved : CCTK_DELTA_TIME;

  MyQuadrature quadrature{};
  quadrature.type = kGauleg;
  quadrature.dim = 1;
  quadrature.nx = 10;
  quadrature.ny = quadrature.nz = 1;
  quadrature.x1 = 0;
  quadrature.x2 = 1;
  GaussLegendre(&quadrature);
  const auto flags = nuX_NuRates::global_opac_flags;
  const auto parameters = nuX_NuRates::global_opac_params;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        if (nuX_m1_mask[ijk]) {
          for (int sp = 0; sp < 3; ++sp) {
            const int index = layout_cc.linear(p.i, p.j, p.k, sp);
            source_update_delta_N[index] = 0;
            source_update_delta_E[index] = 0;
            source_update_delta_Fx[index] = 0;
            source_update_delta_Fy[index] = 0;
            source_update_delta_Fz[index] = 0;
          }
          return;
        }

        tensor::metric<4> g_dd;
        tensor::inv_metric<4> g_uu;
        tensor::generic<CCTK_REAL, 4, 1> n_d, n_u, beta_u;
        tensor::generic<CCTK_REAL, 4, 2> gamma_ud;
        geom.get_metric(p, &g_dd);
        geom.get_inv_metric(p, &g_uu);
        geom.get_normal_form(p, &n_d);
        geom.get_normal(p, &n_u);
        geom.get_shift_vec(p, &beta_u);
        geom.get_space_proj(p, &gamma_ud);
        const CCTK_REAL alpha = geom.get_lapse(p);
        const CCTK_REAL W = fidu_w_lorentz[ijk];
        const CCTK_REAL volume = sqrt(nuX_Utils::metric::spatial_det(
            g_dd(1, 1), g_dd(1, 2), g_dd(1, 3), g_dd(2, 2), g_dd(2, 3), g_dd(3, 3)));
        if (!(isfinite(volume) && volume > 0 && isfinite(alpha) && alpha > 0 &&
              isfinite(W) && W >= 1))
          source_solver_abort();

        tensor::generic<CCTK_REAL, 4, 1> u_u, u_d, v_u, v_d;
        tensor::generic<CCTK_REAL, 4, 2> proj_ud;
        fidu.get(p, &u_u);
        tensor::contract(g_dd, u_u, &u_d);
        calc_proj(u_d, u_u, &proj_ud);
        pack_v_u(fidu_velx[ijk], fidu_vely[ijk], fidu_velz[ijk], &v_u);
        tensor::contract(g_dd, v_u, &v_d);

        nurates_source_state_t base[3], updated[3];
        for (int sp = 0; sp < 3; ++sp) {
          const int index = layout_cc.linear(p.i, p.j, p.k, sp);
          auto &state = base[sp];
          state.N = fmax(rN[index], rad_N_floor);
          state.E = rE[index];
          pack_F_d(beta_u(1), beta_u(2), beta_u(3), rFx[index], rFy[index],
                   rFz[index], &state.F_d);
          // Match the existing repaired provisional baseline. This repair
          // is not part of the physical collision increments below.
          apply_floor(g_uu, &state.E, &state.F_d, rad_E_floor, rad_eps);
          state.chi = chi[index];
        }

        CCTK_REAL rho_cell;
        MyEOSParams matter{};
        if (!prepare_nurates_matter(eos, rho[ijk], temperature[ijk], Ye[ijk],
                                    particle_mass, rho_cell, matter)) {
          printf("nuX_M1: invalid NuRates matter state at (%e,%e,%e)\n",
                 p.x, p.y, p.z);
          source_solver_abort();
        }
        const nurates_cell_rates_t rates{
            &quadrature, matter, flags, parameters, opacity_rate_max};
        const SourceUpdateContext geometry(
            cctkGH, p.i, p.j, p.k, 0, closure_epsilon, closure_maxiter, false,
            dt, alpha, g_dd, g_uu, n_d, n_u, gamma_ud, u_d, u_u, v_d, v_u,
            proj_ud, W, base[0].E, base[0].F_d, base[0].E, base[0].F_d, 0, 0, 0);
        nurates_source_options_t options;
        options.closure = closure_fun;
        options.maxiter = source_maxiter;
        options.inner_maxiter = source_maxiter;
        options.epsrel = source_epsrel;
        options.step_tolerance = source_epsrel;
        options.use_initial_spectrum = true;
        const CCTK_REAL eta = (matter.mu_e + matter.mu_p - matter.mu_n) / matter.temp;
        for (int sp = 0; sp < 3; ++sp) {
          options.initial_log_temperature[sp] = log(matter.temp);
          options.initial_eta[sp] = sp == 0 ? eta : sp == 1 ? -eta : 0;
          // Scattering cannot populate a genuinely empty species. Do not
          // introduce an artificial positive population for logarithmic
          // spectral coordinates when no enabled process can emit it.
          options.active_species[sp] =
              flags.use_pair || flags.use_brem ||
              (sp < 2 && flags.use_abs_em) || base[sp].N > 0 ||
              base[sp].E > 0;
        }
        if (!coupled_nurates_source(geometry, volume, base, rates, options, updated)) {
          report_nurates_source_failure(geometry, matter, options, base,
                                        p.x, p.y, p.z, rho_cell, volume,
                                        alpha * macro_dt / W);
        }
        for (int sp = 0; sp < 3; ++sp) {
          const int index = layout_cc.linear(p.i, p.j, p.k, sp);
          source_update_delta_N[index] = updated[sp].N - base[sp].N;
          source_update_delta_E[index] = updated[sp].E - base[sp].E;
          source_update_delta_Fx[index] = updated[sp].F_d(1) - base[sp].F_d(1);
          source_update_delta_Fy[index] = updated[sp].F_d(2) - base[sp].F_d(2);
          source_update_delta_Fz[index] = updated[sp].F_d(3) - base[sp].F_d(3);
          // Retain the accepted root itself: base + (root - base) can erase
          // a small, physical root in a very stiff cell. The application
          // kernel subsequently subtracts only any rejected limiter fraction.
          if (store_accepted_state) {
            rN[index] = updated[sp].N;
            rE[index] = updated[sp].E;
            unpack_F_d(updated[sp].F_d, &rFx[index], &rFy[index], &rFz[index]);
          }
          chi[index] = updated[sp].chi;
        }
      });
}

} // namespace nuX_M1
