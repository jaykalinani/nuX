// Optional exact-input regression for a production source failure. The same
// callback and solve run on host and device; no diagnostic grid functions.
#define NUX_M1_CLOSURE_IMPLEMENTATION
#define NUX_M1_SOURCES_IMPLEMENTATION
#include "nuX_M1_nurates_source.hxx"
#include "nuX_M1_nurates_cell.hxx"

#include <AMReX_GpuMemory.H>
#include <AMReX_GpuLaunch.H>
#include <cctk_Arguments.h>

namespace nuX_M1 {
namespace {

struct replay_result {
  bool projected, rates_valid, evaluated, jacobian_valid, solved;
  CCTK_REAL number[3], energy[3], chi[3], residual[6];
  CCTK_REAL accepted[3][6]; // N, E, Fx, Fy, Fz, chi
};

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_NOINLINE replay_result replay_source(
    const cGH *const cctkGH, const EOSX::eos_3p_tabulated3d *const eos,
    const MyQuadrature &quadrature, const OpacityFlags flags,
    const OpacityParams parameters) {
  // Cell (3,3,7) from the first source update of job 3465450. Retain the
  // captured round-trip doubles, including the rounded projectors. Rebuilding
  // them from cleaner geometry would no longer replay the failed inputs.
  const CCTK_REAL vectors[6][4] = {
      {-0.803448478920372, 0, 0, 0},
      {1.2446348785721055, 2.867586511069647e-10,
       2.8675864445453695e-10, 2.9427965344176634e-9},
      {-0.8034489295513517, 0.0009210121774478071,
       0.000921012177447795, 1.655958946679565e-5},
      {1.2446355766510397, 0.0006088738571427504,
       0.0006088738571427358, 1.0950347514280765e-5},
      {-4.635476652818659e-13, 0.0009210116608795751,
       0.0009210116608795629, 1.6559580179016333e-5},
      {0, 0.0006088732288849319, 0.0006088732288849238,
       1.0947398576023367e-5}};
  const CCTK_REAL metric[4][4] = {
      {-0.6455294582794594, -3.4850805791020414e-10,
       -3.4850804982527057e-10, -3.576485948535274e-9},
      {-3.4850805791020414e-10, 1.5126492974211423,
       -1.739351433611394e-9, 4.319312066010046e-9},
      {-3.4850804982527057e-10, -1.739351433611394e-9,
       1.512649297421142, 4.319312075695958e-9},
      {-3.576485948535274e-9, 4.319312066010046e-9,
       4.319312075695958e-9, 1.5126493115412456}};
  const CCTK_REAL inverse[4][4] = {
      {-1.5491159809581998, -3.5690981890001784e-10,
       -3.569098106201742e-10, -3.6627072072773415e-9},
      {-3.5690981890001784e-10, 0.6610917690603246,
       7.601702035319188e-10, -1.887722182276762e-9},
      {-3.569098106201742e-10, 7.601702035319188e-10,
       0.6610917690603246, -1.887722186509916e-9},
      {-3.6627072072773415e-9, -1.887722182276762e-9,
       -1.887722186509916e-9, 0.6610917628892419}};
  const CCTK_REAL space_projector[4][4] = {
      {-5.220071983463775e-17, 0, 0, 0},
      {-2.3039580204914843e-10, 1, 0, 0},
      {-2.3039579670426548e-10, 0, 1, 0},
      {-2.3643853993500136e-9, 0, 0, 1}};
  const CCTK_REAL fluid_projector[4][4] = {
      {-1.1217418072550189e-6, 0.001146324522580381,
       0.0011463245225803661, 2.0610654185109686e-5},
      {-0.0004891990487531455, 1.000000560780237,
       5.607802369580822e-7, 1.008270111134833e-8},
      {-0.0004891990487531337, 5.607802369580762e-7,
       1.000000560780237, 1.0082701111348087e-8},
      {-8.798044988564185e-6, 1.008540340793791e-8,
       1.0085403407937776e-8, 1.0000000001813332}};
  tensor::metric<4> gdd{};
  tensor::inv_metric<4> guu{};
  tensor::generic<CCTK_REAL, 4, 1> nd{}, nu{}, ud{}, uu{}, vd{}, vu{};
  tensor::generic<CCTK_REAL, 4, 2> gamma{}, projector{};
  for (int a = 0; a < 4; ++a) {
    nd(a) = vectors[0][a];
    nu(a) = vectors[1][a];
    ud(a) = vectors[2][a];
    uu(a) = vectors[3][a];
    vd(a) = vectors[4][a];
    vu(a) = vectors[5][a];
    for (int b = 0; b < 4; ++b) {
      gdd(a, b) = metric[a][b];
      guu(a, b) = inverse[a][b];
      gamma(a, b) = space_projector[a][b];
      projector(a, b) = fluid_projector[a][b];
    }
  }
  MyEOSParams matter{};
  matter.nb = 1.1171770927888821e17;
  matter.temp = 9.356701812802092;
  matter.ye = 0.05235645319389099;
  matter.mu_p = -88.947931095125;
  matter.mu_n = 18.385077730150464;
  matter.mu_e = 107.28078229700728;
  matter.yn = 0.9476431366192015;
  matter.yp = 0.052356043007251156;
  nurates_source_state_t base[3]{}, projected[3], next[3];
  for (int sp = 0; sp < 3; ++sp) {
    auto &s = base[sp];
    s.N = 1e-10;
    s.E = 1e-15;
    s.F_d(0) = -2.450385764496245e-27;
    s.F_d(1) = 8.638585111745047e-19;
    s.F_d(2) = 8.63858511174498e-19;
    s.F_d(3) = 8.680171263955537e-19;
    s.chi = 0.3333341306435716;
    projected[sp] = s;
  }
  constexpr CCTK_REAL volume = 1.8604044815878096;
  const SourceUpdateContext geometry(
      cctkGH, 3, 3, 7, 0, 1e-5, 1000, false, 0.03125,
      0.803448478920372, gdd, guu, nd, nu, gamma, ud, uu, vd, vu,
      projector, 1.0000005608704579, base[0].E, base[0].F_d,
      base[0].E, base[0].F_d, 0, 0, 0);
  const nurates_cell_rates_t rates{
      eos, &quadrature, matter, flags, parameters,
      0.00029988421889534097, 930.18478708, 0.05021550176813178,
      volume, 1, 1, 3, 1e30, 1e-10, 1e-15};
  nurates_source_options_t options;
  options.epsrel = options.step_tolerance = 1e-5;
  options.maxiter = options.inner_maxiter = 100;
  options.use_initial_spectrum = true;
  nurates_source_detail::spectral_vector x, f;
  const CCTK_REAL eta = (matter.mu_e + matter.mu_p - matter.mu_n) / matter.temp;
  for (int sp = 0; sp < 3; ++sp) {
    x(2 * sp) = options.initial_log_temperature[sp] = log(matter.temp);
    x(2 * sp + 1) = options.initial_eta[sp] = sp == 0 ? eta : sp == 1 ? -eta : 0;
  }
  replay_result result{};
  M1Quantities physical{};
  result.projected = nurates_source_detail::physical_moments(
      geometry, volume, options, projected, physical);
  if (!result.projected)
    return result; // This is also an immediate failure in the actual solve.
  for (int sp = 0; sp < 3; ++sp) {
    result.number[sp] = physical.n[sp];
    result.energy[sp] = physical.J[sp];
    result.chi[sp] = projected[sp].chi;
  }
  nurates_opacities_t coefficients[3];
  result.rates_valid = rates(physical, coefficients);
  // physical_moments updates chi. Use that projected copy, exactly as the
  // production solver does, rather than the original stored closure factor.
  const nurates_source_detail::spectral_source_t<nurates_cell_rates_t> source{
      geometry, volume, projected, rates, options};
  result.evaluated = source.evaluate(x, f);
  if (result.evaluated) {
    for (int a = 0; a < 6; ++a)
      result.residual[a] = f(a);
    nurates_source_detail::spectral_matrix jacobian;
    result.jacobian_valid = source.jacobian(x, f, jacobian);
  }
  // Rejected baseline rates or an LTE trial are not necessarily fatal: the
  // complete solve retries other spectra. Always exercise that full path.
  result.solved = coupled_nurates_source(
      geometry, volume, base, rates, options, next);
  if (result.solved)
    for (int sp = 0; sp < 3; ++sp) {
      result.accepted[sp][0] = next[sp].N;
      result.accepted[sp][1] = next[sp].E;
      for (int a = 1; a < 4; ++a)
        result.accepted[sp][a + 1] = next[sp].F_d(a);
      result.accepted[sp][5] = next[sp].chi;
    }
  return result;
}

void report_replay(const char *const backend, const replay_result &r) {
  CCTK_VINFO("NuRates exact replay %s: projected=%d rates=%d LTE=%d Jacobian=%d solve=%d",
             backend, r.projected, r.rates_valid, r.evaluated,
             r.jacobian_valid, r.solved);
  if (r.projected)
    for (int sp = 0; sp < 3; ++sp)
      CCTK_VINFO("NuRates exact replay %s species %d: n=%.17e J=%.17e chi=%.17e",
                 backend, sp, r.number[sp], r.energy[sp], r.chi[sp]);
  if (r.evaluated)
    CCTK_VINFO("NuRates exact replay %s residual: [%.17e,%.17e,%.17e,%.17e,%.17e,%.17e]",
               backend, r.residual[0], r.residual[1], r.residual[2],
               r.residual[3], r.residual[4], r.residual[5]);
  if (r.solved)
    for (int sp = 0; sp < 3; ++sp) {
      const auto &s = r.accepted[sp];
      CCTK_VINFO("NuRates exact replay %s accepted %d: N=%.17e E=%.17e Fx=%.17e Fy=%.17e Fz=%.17e chi=%.17e",
                 backend, sp, s[0], s[1], s[2], s[3], s[4], s[5]);
    }
}

} // namespace

extern "C" void nuX_NuRates_ReplaySource(CCTK_ARGUMENTS) {
  const auto eos = EOSX::global_eos_3p_tab3d;
  if (!eos || !CCTK_IsThornActive("nuX_NuRates"))
    CCTK_ERROR("NuRates source replay requires active tabulated EOSX and NuRates");
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
  const auto host = replay_source(cctkGH, eos, quadrature, flags, parameters);
  report_replay("host", host);

  // GLOBAL WRAGH has no Loop tile/grid descriptor. A single test-only launch
  // needs only one temporary result, not grid storage or a manufactured tile.
  const amrex::Gpu::LaunchSafeGuard launch_region(true);
  amrex::Gpu::DeviceScalar<replay_result> storage;
  auto *const result = storage.dataPtr();
  amrex::ParallelFor(1, [=] CCTK_DEVICE(int) {
    *result = replay_source(cctkGH, eos, quadrature, flags, parameters);
  });
  amrex::Gpu::streamSynchronize();
  const auto device = storage.dataValue();
#ifdef AMREX_USE_GPU
  report_replay("device", device);
#else
  report_replay("host-parallel", device);
#endif
  if (!host.solved || !device.solved)
    CCTK_ERROR("NuRates exact source replay failed");
  CCTK_INFO("NuRates exact source replay passed");
}

} // namespace nuX_M1
