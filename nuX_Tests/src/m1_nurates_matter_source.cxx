// Optional host integration test: unlike the algebra tests, this exercises
// the actual EOS, full NuRates callback, and trapped-equilibrium feedback.
// No diagnostic grid functions or production-kernel instrumentation.
#define NUX_M1_CLOSURE_IMPLEMENTATION
#define NUX_M1_SOURCES_IMPLEMENTATION
#include "nuX_M1_nurates_source.hxx"
#include "nuX_M1_nurates_cell.hxx"

#include <cctk_Arguments.h>

namespace nuX_M1 {
namespace {

struct counted_rates_t {
  const nurates_cell_rates_t &rates;
  int &calls;
  int &rejected;

  CCTK_HOST CCTK_DEVICE bool operator()(
      const M1Quantities &radiation, nurates_opacities_t *out) const {
    ++calls;
    const bool valid = rates(radiation, out);
    rejected += !valid;
    return valid;
  }
};

} // namespace

extern "C" void nuX_NuRates_TestMatterSource(CCTK_ARGUMENTS) {
  const auto eos = EOSX::global_eos_3p_tab3d;
  if (!eos || !CCTK_IsThornActive("nuX_NuRates"))
    CCTK_ERROR("NuRates matter-source test requires active tabulated EOSX and NuRates");

  // Recorded first-step hot-core matter state. The static conformal geometry
  // below is representative, not a claim of replaying all eight cell vertices.
  constexpr CCTK_REAL particle_mass = 930.18478708;
  CCTK_REAL rho;
  MyEOSParams matter{};
  if (!prepare_nurates_matter(eos, 7.514728e-4, 1.913318e1, 8.769979e-2,
                              particle_mass, rho, matter))
    CCTK_ERROR("Cannot prepare the NuRates matter-source test state");

  MyQuadrature quadrature{};
  quadrature.type = kGauleg;
  quadrature.dim = 1;
  quadrature.nx = 10;
  quadrature.ny = quadrature.nz = 1;
  quadrature.x1 = 0;
  quadrature.x2 = 1;
  GaussLegendre(&quadrature);

  constexpr CCTK_REAL alpha = 0.74, metric = 1.72, dt = 0.03125;
  const CCTK_REAL volume = std::pow(metric, 1.5);
  tensor::metric<4> gdd{};
  tensor::inv_metric<4> guu{};
  tensor::generic<CCTK_REAL, 4, 1> nd{}, nu{}, ud{}, uu{}, vd{}, vu{};
  tensor::generic<CCTK_REAL, 4, 2> projector{}, gamma{};
  gdd(0, 0) = -alpha * alpha;
  guu(0, 0) = -1 / (alpha * alpha);
  nd(0) = ud(0) = -alpha;
  nu(0) = uu(0) = 1 / alpha;
  for (int a = 1; a < 4; ++a) {
    gdd(a, a) = metric;
    guu(a, a) = 1 / metric;
    gamma(a, a) = 1;
  }
  calc_proj(ud, uu, &projector);
  nurates_source_state_t base[3]{}, next[3];
  for (auto &s : base) {
    s.N = 1e-10;
    s.E = 1e-15;
    s.chi = CCTK_REAL(1) / 3;
  }
  const SourceUpdateContext geometry(
      cctkGH, 0, 0, 0, 0, 1e-5, 1000, false, dt, alpha, gdd, guu, nd,
      nu, gamma, ud, uu, vd, vu, projector, 1, base[0].E, base[0].F_d,
      base[0].E, base[0].F_d, 0, 0, 0);
  const nurates_cell_rates_t rates{
      eos, &quadrature, matter, nuX_NuRates::global_opac_flags,
      nuX_NuRates::global_opac_params, rho, particle_mass, 2 * alpha * dt,
      volume, 1, 1, 3, -1, 1e-10, 1e-15};
  int calls = 0, rejected = 0;
  const counted_rates_t counted{rates, calls, rejected};
  nurates_source_options_t options;
  options.epsrel = options.step_tolerance = 1e-5;
  options.maxiter = options.inner_maxiter = 100;
  options.use_initial_spectrum = true;
  nurates_source_detail::spectral_vector x, f;
  const CCTK_REAL eta = (matter.mu_p + matter.mu_e - matter.mu_n) / matter.temp;
  for (int sp = 0; sp < 3; ++sp) {
    x(2 * sp) = options.initial_log_temperature[sp] = std::log(matter.temp);
    x(2 * sp + 1) = options.initial_eta[sp] = sp == 0 ? eta : sp == 1 ? -eta : 0;
  }
  const nurates_source_detail::spectral_source_t<counted_rates_t> source{
      geometry, volume, base, counted, options};
  const bool evaluated = source.evaluate(x, f);
  CCTK_VINFO("NuRates matter-source LTE trial: valid=%d calls=%d rejected=%d",
             evaluated, calls, rejected);
  if (evaluated) {
    CCTK_VINFO("NuRates matter-source residual: [%e,%e,%e,%e,%e,%e]",
               f(0), f(1), f(2), f(3), f(4), f(5));
    nurates_source_detail::spectral_matrix jacobian;
    const bool jacobian_valid = source.jacobian(x, f, jacobian);
    CCTK_VINFO("NuRates matter-source Jacobian: valid=%d calls=%d rejected=%d",
               jacobian_valid, calls, rejected);
  }
  const bool solved = coupled_nurates_source(
      geometry, volume, base, counted, options, next);
  CCTK_VINFO("NuRates matter-source solve: valid=%d calls=%d rejected=%d",
             solved, calls, rejected);
  if (!solved)
    CCTK_ERROR("NuRates production matter-source test failed");
  CCTK_INFO("NuRates production matter-source test passed");
}

} // namespace nuX_M1
