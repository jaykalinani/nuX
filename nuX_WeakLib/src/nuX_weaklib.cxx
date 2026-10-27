#include <algorithm>
#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_ParallelDescriptor.H>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "nuX_M1_closure.hxx"
#include "nuX_rate_units.hxx"
#include "nuX_utils.hxx"
#include "setup_eos.hxx"
#include "wli_io_opacity.H"
#include "nuX_weaklib.hxx"
#include "nuX_weaklib_io.hxx"

namespace nuX_WeakLib {
using namespace Loop;
using namespace nuX_Utils;

namespace {
// No device-owning object is constructed at static initialization or destroyed
// after AMReX finalization. Host source arrays are discarded after upload.
struct Tables {
  TableView view;
  std::array<amrex::Gpu::DeviceVector<double>, 9> axes; // physical E, node indices, weights, 2x(rho,T,Ye)
  wli::ResidentTable<4> emab[2];
  wli::ResidentTable<5> iso[2];
};
std::unique_ptr<Tables> tables;

void require(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error("nuX_WeakLib: " + message);
}

void upload(amrex::Gpu::DeviceVector<double> &dest, const std::vector<double> &src) {
  dest.resize(src.size());
  amrex::Gpu::htod_memcpy(dest.data(), src.data(), src.size() * sizeof(double));
}

AxisView axis(amrex::Gpu::DeviceVector<double> &dest,
              const wli::io::HostAxis &src, bool logarithmic) {
  require(src.kind == (logarithmic ? wli::AxisKind::Log : wli::AxisKind::Linear),
          "unexpected interpolation kind for " + src.name);
  require(src.points.size() >= 2, "axis needs at least two nodes");
  std::vector<double> nodes = src.points;
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    require(std::isfinite(nodes[i]) && (!logarithmic || nodes[i] > 0) &&
                (i == 0 || nodes[i] > src.points[i - 1]),
            "non-finite/non-monotone axis " + src.name);
    if (logarithmic) nodes[i] = std::log10(nodes[i]);
  }
  require(coordinate_axis_is_uniform(nodes.data(), int(nodes.size())),
          "WLI requires uniformly spaced interpolation coordinates for " + src.name);
  upload(dest, nodes);
  return {dest.data(), int(nodes.size()), src.points.front(), src.points.back()};
}

ThermoView thermo(Tables &owner, const wli::io::HostOpacityCommon &common, int start) {
  // The tables' physical shape and unit convention are rho,T,Ye. Do not guess
  // when an unsupported table layout or coordinate convention is supplied.
  require(common.tsIndices == std::array<int, 3>{1, 2, 3} &&
              common.axes[0].name == "Density" && common.axes[1].name == "Temperature" &&
              common.axes[2].name == "Electron Fraction", "unsupported thermodynamic axis order");
  return {axis(owner.axes[start], common.axes[0], true),
          axis(owner.axes[start + 1], common.axes[1], true),
          axis(owner.axes[start + 2], common.axes[2], false)};
}

void energy_grid(Tables &owner, const wli::io::HostOpacityCommon &common,
                 std::vector<double> &reference) {
  const auto &nodes = common.energyGrid.points;
  require(common.energyGrid.kind == wli::AxisKind::Log &&
              physical_energy_nodes_are_valid(nodes.data(), nodes.size()),
          "expected a declared-log MeV energy grid with 2..INT_MAX finite positive increasing nodes");
  if (!reference.empty()) {
    require(nodes == reference, "enabled channels must have identical energy nodes");
    return;
  }
  reference = nodes;
  std::vector<double> node_coordinates(nodes.size()), weights(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    node_coordinates[i] = double(i); // exact-node VALUE access, not physical log(E)
    // Composite trapezoidal quadrature on physical E, not log E; no implied
    // opacity extrapolation beyond the first/last tabulated node.
    const auto left = i == 0 ? 0 : i - 1;
    const auto right = i + 1 == nodes.size() ? i : i + 1;
    weights[i] = (nodes[right] - nodes[left]) / 2;
  }
  upload(owner.axes[0], nodes);
  upload(owner.axes[1], node_coordinates); // REPLACES the old log-energy buffer
  upload(owner.axes[2], weights);
  owner.view.energies = owner.axes[0].data();
  owner.view.energy_node_coordinates = owner.axes[1].data();
  owner.view.weights = owner.axes[2].data();
  owner.view.ne = int(nodes.size());
}

template <int D>
void table(wli::ResidentTable<D> &dest, const std::array<int, D> &shape,
           const std::vector<double> &values) {
  wli::TableMeta<D> meta{};
  std::size_t size = 1;
  for (int d = 0; d < D; ++d) {
    require(shape[d] >= 1, "invalid table extent");
    meta.n[d] = shape[d];
    size *= std::size_t(shape[d]);
  }
  require(values.size() == size &&
              std::all_of(values.begin(), values.end(), [](double x) { return std::isfinite(x); }),
          "invalid table data size or non-finite log-stored values");
  dest.upload(meta, values.data(), values.size());
}

void finish_errors(const amrex::Gpu::DeviceScalar<int> &error) {
  // AMReX's device Abort/assert can be compiled out under NDEBUG. A transient
  // integer flag plus synchronized host error works in release GPU builds too.
  amrex::Gpu::synchronize();
  const int status = error.dataValue();
  if (status)
    CCTK_VERROR("nuX_WeakLib rejected a grid point (status %d): 1 invalid state, "
                "2 outside opacity table, 3 insufficient spectral coverage, "
                "4 invalid kernel, 5 invalid moment, 6 outside EOSX, 7 invalid geometry. "
                "No coefficient or thermodynamic state has been clamped.", status);
}
} // namespace

extern "C" void nuX_WeakLib_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_WeakLib_ParamCheck;
  DECLARE_CCTK_PARAMETERS;
  if (!CCTK_Equals(rates_lib, "WeakLib")) return;
  if (nspecies != 3 || ngroups != 1)
    CCTK_PARAMWARN("WeakLib requires exactly three grey transport species");
  if (!include_emab && !include_iso)
    CCTK_PARAMWARN("WeakLib needs at least one enabled channel");
  if (include_pair || include_nes || include_brem || include_nuclei_ec)
    CCTK_PARAMWARN("WeakLib Pair/NES/Brem and EC_table closures are not implemented; "
                   "enabling their interpolators alone would not define conservative collision terms");
  if ((include_emab && !*emab_table) || (include_iso && !*iso_table))
    CCTK_PARAMWARN("WeakLib requires a nonempty file path for every enabled channel");
  if (opacity_tau_trap >= 0 || opacity_corr_fac_max != 1)
    CCTK_PARAMWARN("WeakLib requires opacity_tau_trap=-1 and opacity_corr_fac_max=1; "
                   "its rates are direct equilibrium-spectrum table integrals");
}

extern "C" void nuX_WeakLib_Load(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_WeakLib_Load;
  DECLARE_CCTK_PARAMETERS;
  if (!CCTK_Equals(rates_lib, "WeakLib")) return;
  if (tables) CCTK_ERROR("WeakLib tables were already loaded");
  if (!EOSX::global_eos_3p_tab3d)
    CCTK_ERROR("nuX_WeakLib requires EOSX's tabulated evolution EOS");
  try {
    require(nspecies == 3 && ngroups == 1, "only three-species grey transport is supported");
    require(include_emab || include_iso, "no interactions enabled");
    require(!include_pair && !include_nes && !include_brem && !include_nuclei_ec,
            "Pair/NES/Brem/EC_table coupling is not implemented");
    require(opacity_tau_trap < 0 && opacity_corr_fac_max == 1,
            "disable the other providers' trapped-equilibrium and energy-squared corrections");
    auto owner = std::make_unique<Tables>();
    owner->view.tail_tolerance = spectral_tail_tolerance;
    std::vector<double> energies;
    if (include_emab) {
      host_io::check_emab_file_collective(emab_table);
      const auto host = wli::io::read_emab_table(emab_table);
      if (host.hasECTable)
        CCTK_INFO("WeakLib EmAb file contains an optional EC_table; include_nuclei_ec=no, so only its main EmAb arrays are used");
      require(host.nOpacities == 2 && host.values.size() == 2 && host.offset.size() == 2,
              "EmAb must contain exactly nue and anti-nue");
      energy_grid(*owner, host.common, energies);
      owner->view.emab_thermo = thermo(*owner, host.common, 3);
      const auto t = owner->view.emab_thermo;
      require(host.nPoints == std::array<int, 4>{owner->view.ne, t.rho.n, t.temp.n, t.ye.n},
              "EmAb shape disagrees with its coordinate axes");
      for (int s = 0; s < 2; ++s) {
        require(std::isfinite(host.offset[s]), "non-finite EmAb offset");
        table<4>(owner->emab[s], host.nPoints, host.values[s]);
        owner->view.emab[s] = owner->emab[s].view();
        owner->view.emab_offsets[s] = host.offset[s];
      }
      owner->view.include_emab = true;
    }
    if (include_iso) {
      const auto metadata = host_io::check_iso_file_collective(iso_table, iso_table_convention);
      CCTK_VINFO("WeakLib Iso convention=%s, Units=['%s','%s']; raw Legendre expansion "
                 "coefficients interpreted as cm^-1 MeV^-2; user-declared convention, "
                 "EOS/model provenance not verified",
                 host_io::iso_convention_name(metadata.convention),
                 metadata.units[0].data(), metadata.units[1].data());
      if (metadata.convention == host_io::IsoConvention::fullweakphysics_5836de98)
        CCTK_VINFO("WeakLib Iso native metadata: weak_magnetism_corr=%d ion_ion_corr=%d "
                   "many_body_corr=%d ga_strange=%.17g (recorded, not changed)",
                   metadata.corrections[0], metadata.corrections[1],
                   metadata.corrections[2], metadata.ga_strange);
      const auto host = wli::io::read_scat_iso_table(iso_table);
      require(host.nOpacities == 2 && host.nMoments == 2 &&
                  host.values.size() == 2 && host.offset.size() == 4,
              "Iso must contain neutrino/antineutrino tables and two Legendre moments");
      energy_grid(*owner, host.common, energies);
      owner->view.iso_thermo = thermo(*owner, host.common, 6);
      const auto t = owner->view.iso_thermo;
      require(host.nPoints == std::array<int, 5>{owner->view.ne, 2, t.rho.n, t.temp.n, t.ye.n},
              "Iso shape disagrees with its coordinate axes");
      for (int s = 0; s < 2; ++s) {
        table<5>(owner->iso[s], host.nPoints, host.values[s]);
        owner->view.iso[s] = owner->iso[s].view();
        for (int m = 0; m < 2; ++m) {
          const double offset = wli::IsoOffset(host.offset.data(), 2, 2, s, m);
          require(std::isfinite(offset), "non-finite Iso offset");
          owner->view.iso_offsets[s][m] = offset;
        }
      }
      owner->view.include_iso = true;
    }
    amrex::Gpu::synchronize();
    if (run_self_tests) self_test();
    tables = std::move(owner);
    CCTK_VINFO("WeakLib tables resident: EmAb=%d Iso=%d, %d energy nodes; EOS remains EOSX",
               int(include_emab), int(include_iso), tables->view.ne);
  } catch (const std::exception &e) {
    CCTK_VERROR("WeakLib table setup failed: %s", e.what());
  }
}

extern "C" void nuX_WeakLib_Free(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_WeakLib_Free;
  if (tables) {
    amrex::Gpu::synchronize();
    tables.reset();
  }
}

extern "C" void nuX_WeakLib_CalcOpacity(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_WeakLib_CalcOpacity;
  DECLARE_CCTK_PARAMETERS;
  if (!tables || !EOSX::global_eos_3p_tab3d)
    CCTK_ERROR("WeakLib provider/EOSX is not initialized; select the provider before setup");
  const amrex::Gpu::LaunchSafeGuard launch_region(true);
  const auto view = tables->view;
  const auto eos = EOSX::global_eos_3p_tab3d;
  const bool rest_subtracted = CCTK_Equals(chemical_potentials, "rest-mass-subtracted");
  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout(cctkGH, {1, 1, 1});
  amrex::Gpu::DeviceScalar<int> error(0);
  auto *error_ptr = error.dataPtr();
  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const int ijk = layout.linear(p.i, p.j, p.k);
    if (nuX_m1_mask[ijk]) {
      for (int s = 0; s < 3; ++s) {
        const int i = layout.linear(p.i, p.j, p.k, s);
        abs_0[i] = abs_1[i] = eta_0[i] = eta_1[i] = scat_1[i] = nueave[i] = 0;
      }
      return;
    }
    const auto d = rho[ijk], t = temperature[ijk], y = Ye[ijk];
    if (!(std::isfinite(d) && d > 0 && std::isfinite(t) && t > 0 && std::isfinite(y) && y >= 0 && y <= 1 &&
          d >= eos->rgrho.min && d <= eos->rgrho.max &&
          t >= eos->rgtemp.min && t <= eos->rgtemp.max &&
          y >= eos->rgye.min && y <= eos->rgye.max)) {
      amrex::Gpu::Atomic::Max(error_ptr, 6);
      return;
    }
    CCTK_REAL mup, mun, mue;
    eos->mu_pne_from_rho_temp_ye(d, t, y, mup, mun, mue);
    Rates rates;
    const auto status = evaluate(view, d * rate_units::code_density_to_g_nm3 /
        rate_units::per_cm3_to_per_nm3, t, y,
        neutrino_chemical_potential(mup, mun, mue, rest_subtracted), rates);
    if (status != Status::success) {
      amrex::Gpu::Atomic::Max(error_ptr, int(status));
      return;
    }
    for (int s = 0; s < 3; ++s) {
      const int i = layout.linear(p.i, p.j, p.k, s);
      constexpr double inverse_cm_to_code = rate_units::per_cm_to_per_nm * rate_units::code_length_to_nm;
      constexpr double number_cgs_to_transport = rate_units::per_cm3_to_per_nm3 /
          rate_units::transport_number_density_to_nm3;
      constexpr double energy_cgs_to_code = rate_units::per_cm3_to_per_nm3 /
          rate_units::code_energy_density_to_mev_nm3;
      // Use the same length convention for emission and absorption. Convert
      // the spectral emission integrals directly, retaining finite emission
      // even if the LTE population rounded to zero before multiplication by
      // a large opacity.
      abs_0[i] = rates.absorption_number[s] * inverse_cm_to_code;
      abs_1[i] = rates.absorption_energy[s] * inverse_cm_to_code;
      scat_1[i] = rates.scattering[s] * inverse_cm_to_code;
      eta_0[i] = rates.emission_number[s] *
          (inverse_cm_to_code * number_cgs_to_transport);
      eta_1[i] = rates.emission_energy[s] *
          (inverse_cm_to_code * energy_cgs_to_code);
      nueave[i] = rates.mean_energy[s] *
          rate_units::transport_number_density_to_nm3 /
          rate_units::code_energy_density_to_mev_nm3; // legacy transport unit, not MeV
      if (!nonnegative_finite(abs_0[i]) || !nonnegative_finite(abs_1[i]) ||
          !nonnegative_finite(scat_1[i]) || !nonnegative_finite(eta_0[i]) ||
          !nonnegative_finite(eta_1[i]) || !nonnegative_finite(nueave[i]))
        amrex::Gpu::Atomic::Max(error_ptr, int(Status::invalid_moment));
    }
  });
  finish_errors(error);
}

extern "C" void nuX_WeakLib_SetEquilibrium(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_WeakLib_SetEquilibrium;
  DECLARE_CCTK_PARAMETERS;
  if (!tables || !EOSX::global_eos_3p_tab3d) CCTK_ERROR("WeakLib/EOSX not initialized");
  const amrex::Gpu::LaunchSafeGuard launch_region(true);
  const auto view = tables->view;
  const auto eos = EOSX::global_eos_3p_tab3d;
  const bool rest_subtracted = CCTK_Equals(chemical_potentials, "rest-mass-subtracted");
  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout cc(cctkGH, {1, 1, 1}), vc(cctkGH, {0, 0, 0});
  const tensor::slicing_geometry_const geom(vc, cc, alp, betax, betay, betaz,
      gxx, gxy, gxz, gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz);
  const tensor::fluid_velocity_field_const fidu(vc, cc, alp, betax, betay, betaz,
      fidu_w_lorentz, fidu_velx, fidu_vely, fidu_velz);
  amrex::Gpu::DeviceScalar<int> error(0);
  auto *error_ptr = error.dataPtr();
  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    const int ijk = cc.linear(p.i, p.j, p.k);
    if (nuX_m1_mask[ijk]) return;
    const auto d = rho[ijk], t = temperature[ijk], y = Ye[ijk];
    const double dcgs = d * rate_units::code_density_to_g_nm3 / rate_units::per_cm3_to_per_nm3;
    if (std::isfinite(dcgs) && dcgs < equilibrium_rho_min) return;
    if (!(std::isfinite(d) && d > 0 && std::isfinite(t) && t > 0 && std::isfinite(y) && y >= 0 && y <= 1 &&
          d >= eos->rgrho.min && d <= eos->rgrho.max &&
          t >= eos->rgtemp.min && t <= eos->rgtemp.max &&
          y >= eos->rgye.min && y <= eos->rgye.max)) {
      amrex::Gpu::Atomic::Max(error_ptr, 6);
      return;
    }
    CCTK_REAL mup, mun, mue;
    eos->mu_pne_from_rho_temp_ye(d, t, y, mup, mun, mue);
    Rates rates;
    const auto status = evaluate(view, dcgs, t, y,
        neutrino_chemical_potential(mup, mun, mue, rest_subtracted), rates);
    if (status != Status::success) {
      amrex::Gpu::Atomic::Max(error_ptr, int(status));
      return;
    }
    tensor::metric<4> gdd;
    tensor::inv_metric<4> guu;
    tensor::generic<CCTK_REAL, 4, 1> nu, ud, uu;
    tensor::generic<CCTK_REAL, 4, 2> projector;
    geom.get_metric(p, &gdd);
    geom.get_inv_metric(p, &guu);
    geom.get_normal(p, &nu);
    geom.get_space_proj(p, &projector);
    fidu.get(p, &uu);
    tensor::contract(gdd, uu, &ud);
    const CCTK_REAL det = metric::spatial_det(gdd(1,1), gdd(1,2), gdd(1,3),
                                              gdd(2,2), gdd(2,3), gdd(3,3));
    if (!(std::isfinite(det) && det > 0 && std::isfinite(fidu_w_lorentz[ijk]) && fidu_w_lorentz[ijk] >= 1)) {
      amrex::Gpu::Atomic::Max(error_ptr, 7);
      return;
    }
    const CCTK_REAL volume = std::sqrt(det);
    for (int s = 0; s < 3; ++s) {
      const int i = cc.linear(p.i, p.j, p.k, s);
      const CCTK_REAL n = rates.number[s] * rate_units::per_cm3_to_per_nm3 /
          rate_units::transport_number_density_to_nm3;
      const CCTK_REAL j = volume * rates.energy[s] * rate_units::per_cm3_to_per_nm3 /
          rate_units::code_energy_density_to_mev_nm3;
      tensor::generic<CCTK_REAL, 4, 1> h, f;
      tensor::symmetric2<CCTK_REAL, 4, 2> k, stress;
      for (int a = 0; a < 4; ++a) {
        h(a) = 0;
        for (int b = a; b < 4; ++b) k(a,b) = j * (gdd(a,b) + ud(a) * ud(b)) / 3;
      }
      nuX_M1::assemble_rT(ud, j, h, k, &stress);
      CCTK_REAL e = nuX_M1::calc_J_from_rT(stress, nu);
      nuX_M1::calc_H_from_rT(stress, nu, projector, &f);
      nuX_M1::apply_floor(guu, &e, &f, rad_E_floor, rad_eps);
      rE[i] = e;
      nuX_M1::unpack_F_d(f, &rFx[i], &rFy[i], &rFz[i]);
      rN[i] = std::fmax(volume * n * fidu_w_lorentz[ijk], CCTK_REAL(rad_N_floor));
    }
  });
  finish_errors(error);
}
} // namespace nuX_WeakLib
