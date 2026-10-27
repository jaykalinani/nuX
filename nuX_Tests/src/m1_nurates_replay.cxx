// Optional exact-input regression for a production source failure. The same
// callback and solve run on host and device; no diagnostic grid functions.
#define NUX_M1_CLOSURE_IMPLEMENTATION
#define NUX_M1_SOURCES_IMPLEMENTATION
#include "nuX_M1_nurates_source.hxx"
#include "nuX_M1_nurates_cell.hxx"

#include <AMReX.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_ParallelDescriptor.H>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <locale>
#include <sstream>
#include <string>
#include <type_traits>

namespace nuX_M1 {
namespace {

// Test-only input: no pointers, strings, EOS objects, or reconstructed geometry.
// A copy can be broadcast to MPI ranks and captured by the one-thread launch.
struct replay_input {
  int cell[3];
  CCTK_REAL coordinates[3];
  CCTK_REAL dt, macro_proper_dt, alpha, W, volume, rho, temperature, ye;
  CCTK_REAL matter[6]; // nb, mu_p, mu_n, mu_e, xn, xp
  int closure, closure_maxiter, maxiter, inner_maxiter, angular_maxiter;
  CCTK_REAL closure_epsilon, epsrel;
  CCTK_REAL vectors[6][4]; // nd, nu, ud, uu, vd, vu
  CCTK_REAL metric[4][4], inverse[4][4];
  CCTK_REAL space_projector[4][4], fluid_projector[4][4];
  CCTK_REAL species[3][7]; // N, E, Ft, Fx, Fy, Fz, chi
  CCTK_REAL particle_mass, tau_trap, tau_delta, correction_max, rate_max;
  CCTK_REAL number_floor, energy_floor;
};
static_assert(std::is_trivially_copyable_v<replay_input>);
static_assert(sizeof(CCTK_REAL) == sizeof(double),
              "Exact replay fixtures contain round-trip double values");

struct replay_result {
  bool projected, rates_valid, evaluated, jacobian_valid, solved;
  CCTK_REAL number[3], energy[3], chi[3], residual[6];
  CCTK_REAL accepted[3][6]; // N, E, Fx, Fy, Fz, chi
};

replay_input default_replay_input() {
  replay_input input{};
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
  input.cell[0] = input.cell[1] = 3;
  input.cell[2] = 7;
  // Coordinates were not inputs to the old replay; retain an explicit neutral
  // metadata value rather than inventing a physical position for that case.
  input.dt = 0.03125;
  input.macro_proper_dt = 0.05021550176813178;
  input.alpha = 0.803448478920372;
  input.W = 1.0000005608704579;
  input.volume = 1.8604044815878096;
  input.rho = 0.00029988421889534097;
  input.temperature = 9.356701812802092;
  input.ye = 0.05235645319389099;
  input.matter[0] = 1.1171770927888821e17;
  input.matter[1] = -88.947931095125;
  input.matter[2] = 18.385077730150464;
  input.matter[3] = 107.28078229700728;
  input.matter[4] = 0.9476431366192015;
  input.matter[5] = 0.052356043007251156;
  input.closure = CLOSURE_MINERBO;
  input.closure_epsilon = input.epsrel = 1e-5;
  input.closure_maxiter = 1000;
  input.maxiter = input.inner_maxiter = 100;
  input.angular_maxiter = 32;
  input.particle_mass = 930.18478708;
  input.tau_trap = input.tau_delta = 1;
  input.correction_max = 3;
  input.rate_max = 1e30;
  input.number_floor = 1e-10;
  input.energy_floor = 1e-15;
  for (int a = 0; a < 4; ++a) {
    for (int v = 0; v < 6; ++v)
      input.vectors[v][a] = vectors[v][a];
    for (int b = 0; b < 4; ++b) {
      input.metric[a][b] = metric[a][b];
      input.inverse[a][b] = inverse[a][b];
      input.space_projector[a][b] = space_projector[a][b];
      input.fluid_projector[a][b] = fluid_projector[a][b];
    }
  }
  for (int sp = 0; sp < 3; ++sp) {
    input.species[sp][0] = 1e-10;
    input.species[sp][1] = 1e-15;
    input.species[sp][2] = -2.450385764496245e-27;
    input.species[sp][3] = 8.638585111745047e-19;
    input.species[sp][4] = 8.63858511174498e-19;
    input.species[sp][5] = 8.680171263955537e-19;
    input.species[sp][6] = 0.3333341306435716;
  }
  return input;
}

struct fixture_reader {
  std::ifstream file;
  int line_number = 0;
  std::string error;

  explicit fixture_reader(const char *const path) : file(path) {
    file.imbue(std::locale::classic());
    if (!file)
      error = "Cannot open NuRates source replay fixture";
  }

  template <typename Value>
  static bool value(std::istringstream &fields, Value &out) {
    std::string token;
    if (!(fields >> token))
      return false;
    std::istringstream number(token);
    number.imbue(std::locale::classic());
    if (!(number >> out) || !std::isfinite(static_cast<double>(out)))
      return false;
    number >> std::ws;
    return number.eof(); // Reject e.g. an integer token "2.1", not just its tail.
  }

  template <typename... Values>
  bool row(const char *const label, Values &...values) {
    if (!error.empty())
      return false;
    ++line_number;
    char line[1024] = {};
    std::string actual;
    if (file.getline(line, sizeof(line))) {
      const auto length = file.gcount() - (file.eof() ? 0 : 1);
      if (std::memchr(line, '\0', std::size_t(length))) {
        error = "Embedded NUL in fixture line " + std::to_string(line_number);
        return false;
      }
      std::istringstream fields{std::string(line, std::size_t(length))};
      fields.imbue(std::locale::classic());
      if ((fields >> actual) && actual == label &&
          (value(fields, values) && ...)) {
        fields >> std::ws;
        if (fields.eof())
          return true;
      }
    }
    error = "Invalid fixture line " + std::to_string(line_number) +
        ": expected " + label + " and " +
        std::to_string(sizeof...(Values)) + " finite values";
    return false;
  }

  bool finish() {
    if (!error.empty())
      return false;
    file >> std::ws;
    if (file.bad() || !file.eof()) {
      error = "Unexpected trailing data in NuRates source replay fixture";
      return false;
    }
    return true;
  }
};

bool valid_replay_input(const replay_input &in) {
  if (!(in.dt > 0 && in.macro_proper_dt >= 0 && in.alpha > 0 &&
        in.W >= 1 && in.volume > 0 && in.rho > 0 && in.temperature > 0 &&
        in.ye >= 0 && in.ye <= 1 && in.matter[0] > 0 &&
        in.matter[4] >= 0 && in.matter[4] <= 1 &&
        in.matter[5] >= 0 && in.matter[5] <= 1 &&
        in.closure >= CLOSURE_EDDINGTON && in.closure <= CLOSURE_THIN &&
        in.closure_epsilon > 0 && in.closure_maxiter > 0 && in.epsrel > 0 &&
        in.maxiter > 0 && in.inner_maxiter > 0 && in.angular_maxiter > 0 &&
        in.particle_mass > 0 && (in.tau_trap == -1 || in.tau_trap >= 0) &&
        in.tau_delta >= 0 && in.correction_max >= 1 &&
        (in.rate_max == -1 || in.rate_max > 0) &&
        in.number_floor >= 0 && in.energy_floor >= 0))
    return false;
  for (const auto &s : in.species)
    if (!(s[0] >= 0 && s[1] >= 0 && s[6] >= CCTK_REAL(1) / 3 && s[6] <= 1))
      return false;
  return true;
}

bool read_replay_input(const char *const path, replay_input &in,
                       std::string &error) {
  fixture_reader r(path);
  bool ok = r.row("nux_nurates_source_replay_v1") &&
      r.row("cell", in.cell[0], in.cell[1], in.cell[2]) &&
      r.row("coordinates", in.coordinates[0], in.coordinates[1], in.coordinates[2]) &&
      r.row("scalars", in.dt, in.macro_proper_dt, in.alpha, in.W, in.volume,
            in.rho, in.temperature, in.ye) &&
      r.row("matter", in.matter[0], in.matter[1], in.matter[2],
            in.matter[3], in.matter[4], in.matter[5]) &&
      r.row("solver", in.closure, in.closure_epsilon, in.closure_maxiter,
            in.epsrel, in.maxiter, in.inner_maxiter, in.angular_maxiter);
  const char *const vector_labels[] = {"nd", "nu", "ud", "uu", "vd", "vu"};
  for (int v = 0; ok && v < 6; ++v)
    ok = r.row(vector_labels[v], in.vectors[v][0], in.vectors[v][1],
               in.vectors[v][2], in.vectors[v][3]);
  CCTK_REAL (*const matrices[])[4] = {
      in.metric, in.inverse, in.space_projector, in.fluid_projector};
  const char *const matrix_labels[] = {"gdd", "guu", "gamma", "proj"};
  for (int m = 0; ok && m < 4; ++m)
    for (int a = 0; ok && a < 4; ++a) {
      const auto label = std::string(matrix_labels[m]) + std::to_string(a);
      const auto row = matrices[m][a];
      ok = r.row(label.c_str(), row[0], row[1], row[2], row[3]);
    }
  for (int sp = 0; ok && sp < 3; ++sp) {
    const auto label = "species" + std::to_string(sp);
    const auto s = in.species[sp];
    ok = r.row(label.c_str(), s[0], s[1], s[2], s[3], s[4], s[5], s[6]);
  }
  ok = ok && r.row("rate_options", in.particle_mass, in.tau_trap, in.tau_delta,
                    in.correction_max, in.rate_max, in.number_floor,
                    in.energy_floor) && r.finish();
  error = r.error;
  if (ok && !valid_replay_input(in)) {
    error = "Out-of-range matter, solver, moment, or rate option in replay fixture";
    ok = false;
  }
  return ok;
}

replay_input load_replay_input(const char *const path) {
  replay_input input = default_replay_input();
  int loaded = 1;
  char message[256] = {};
  const int root = amrex::ParallelDescriptor::IOProcessorNumber();
  if (amrex::ParallelDescriptor::IOProcessor() && path[0]) {
    try {
      std::string error;
      loaded = read_replay_input(path, input, error);
      if (!loaded)
        std::snprintf(message, sizeof(message), "%s", error.c_str());
    } catch (const std::exception &e) {
      loaded = 0;
      std::snprintf(message, sizeof(message), "%s", e.what());
    } catch (...) {
      loaded = 0;
      std::snprintf(message, sizeof(message), "%s", "Unknown fixture-read exception");
    }
  }
  // All ranks learn of a malformed/missing root-side file before any abort.
  amrex::ParallelDescriptor::Bcast(&loaded, 1, root);
  amrex::ParallelDescriptor::Bcast(message, sizeof(message), root);
  if (!loaded)
    CCTK_VERROR("NuRates source replay fixture '%s': %s", path, message);
  amrex::ParallelDescriptor::Bcast(
      reinterpret_cast<unsigned char *>(&input), sizeof(input), root);
  return input;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_NOINLINE replay_result replay_source(
    const cGH *const cctkGH, const EOSX::eos_3p_tabulated3d *const eos,
    const MyQuadrature &quadrature, const OpacityFlags flags,
    const OpacityParams parameters, const replay_input &input) {
  tensor::metric<4> gdd{};
  tensor::inv_metric<4> guu{};
  tensor::generic<CCTK_REAL, 4, 1> nd{}, nu{}, ud{}, uu{}, vd{}, vu{};
  tensor::generic<CCTK_REAL, 4, 2> gamma{}, projector{};
  for (int a = 0; a < 4; ++a) {
    nd(a) = input.vectors[0][a];
    nu(a) = input.vectors[1][a];
    ud(a) = input.vectors[2][a];
    uu(a) = input.vectors[3][a];
    vd(a) = input.vectors[4][a];
    vu(a) = input.vectors[5][a];
    for (int b = 0; b < 4; ++b) {
      gdd(a, b) = input.metric[a][b];
      guu(a, b) = input.inverse[a][b];
      gamma(a, b) = input.space_projector[a][b];
      projector(a, b) = input.fluid_projector[a][b];
    }
  }
  MyEOSParams matter{};
  matter.nb = input.matter[0];
  matter.temp = input.temperature;
  matter.ye = input.ye;
  matter.mu_p = input.matter[1];
  matter.mu_n = input.matter[2];
  matter.mu_e = input.matter[3];
  matter.yn = input.matter[4];
  matter.yp = input.matter[5];
  nurates_source_state_t base[3]{}, projected[3], next[3];
  for (int sp = 0; sp < 3; ++sp) {
    auto &s = base[sp];
    s.N = input.species[sp][0];
    s.E = input.species[sp][1];
    for (int a = 0; a < 4; ++a)
      s.F_d(a) = input.species[sp][a + 2];
    s.chi = input.species[sp][6];
    projected[sp] = s;
  }
  const CCTK_REAL volume = input.volume;
  const SourceUpdateContext geometry(
      cctkGH, input.cell[0], input.cell[1], input.cell[2], 0,
      input.closure_epsilon, input.closure_maxiter, false, input.dt,
      input.alpha, gdd, guu, nd, nu, gamma, ud, uu, vd, vu,
      projector, input.W, base[0].E, base[0].F_d,
      base[0].E, base[0].F_d, 0, 0, 0);
  const nurates_cell_rates_t rates{
      eos, &quadrature, matter, flags, parameters,
      input.rho, input.particle_mass, input.macro_proper_dt,
      volume, input.tau_trap, input.tau_delta, input.correction_max,
      input.rate_max, input.number_floor, input.energy_floor};
  nurates_source_options_t options;
  options.closure = static_cast<closure_t>(input.closure);
  options.epsrel = options.step_tolerance = input.epsrel;
  options.maxiter = input.maxiter;
  options.inner_maxiter = input.inner_maxiter;
  options.angular_maxiter = input.angular_maxiter;
  options.use_initial_spectrum = true;
  nurates_source_detail::spectral_vector x, f;
  const CCTK_REAL eta = (matter.mu_e + matter.mu_p - matter.mu_n) / matter.temp;
  for (int sp = 0; sp < 3; ++sp) {
    options.active_species[sp] = flags.use_pair || flags.use_brem ||
        (sp < 2 && flags.use_abs_em) || base[sp].N > 0 || base[sp].E > 0;
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
  DECLARE_CCTK_PARAMETERS;
  const auto eos = EOSX::global_eos_3p_tab3d;
  if (!eos || !CCTK_IsThornActive("nuX_NuRates") || !amrex::Initialized())
    CCTK_ERROR("NuRates source replay requires active tabulated EOSX, NuRates, and an initialized driver");
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
  const auto input = load_replay_input(nurates_source_replay_file);
  CCTK_VINFO("NuRates exact replay fixture: %s; cell=(%d,%d,%d), xyz=(%.17e,%.17e,%.17e), dt=%.17e, macro_proper_dt=%.17e",
             nurates_source_replay_file[0] ? nurates_source_replay_file :
                 "built-in job3465450 (coordinates unavailable)",
             input.cell[0], input.cell[1], input.cell[2], input.coordinates[0],
             input.coordinates[1], input.coordinates[2], input.dt,
             input.macro_proper_dt);
  CCTK_VINFO("NuRates exact replay reactions: beta=%d pair=%d brem=%d inelastic=%d elastic=%d; brem_implementation=%d",
             flags.use_abs_em, flags.use_pair, flags.use_brem,
             flags.use_inelastic_scatt, flags.use_iso,
             static_cast<int>(parameters.brem_implementation));
  CCTK_VINFO("NuRates exact replay corrections: dU=%d dm_eff=%d WM_abs=%d WM_scat=%d decay=%d NN_medium=%d neglect_blocking=%d beta_fallback=%d beta_nb_threshold=%.17e",
             parameters.use_dU, parameters.use_dm_eff, parameters.use_WM_ab,
             parameters.use_WM_sc, parameters.use_decay,
             parameters.use_NN_medium_corr, parameters.neglect_blocking,
             parameters.use_beta_low_density_fallback,
             parameters.beta_low_density_nb_threshold);
  CCTK_VINFO("NuRates exact replay rate options: particle_mass=%.17e tau_trap=%.17e tau_delta=%.17e correction_max=%.17e rate_max=%.17e N_floor=%.17e E_floor=%.17e",
             input.particle_mass, input.tau_trap, input.tau_delta,
             input.correction_max, input.rate_max, input.number_floor,
             input.energy_floor);
  const auto host = replay_source(cctkGH, eos, quadrature, flags, parameters, input);
  report_replay("host", host);

  // GLOBAL WRAGH has no Loop tile/grid descriptor. A single test-only launch
  // needs only one temporary result, not grid storage or a manufactured tile.
  const amrex::Gpu::LaunchSafeGuard launch_region(true);
  amrex::Gpu::DeviceScalar<replay_result> storage;
  auto *const result = storage.dataPtr();
  amrex::ParallelFor(1, [=] CCTK_DEVICE(int) {
    *result = replay_source(cctkGH, eos, quadrature, flags, parameters, input);
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
