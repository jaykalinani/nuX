// Opt-in HOST-ONLY actual-HDF fixture sources. Not an executed test result.
#include "nuX_weaklib_io.hxx"
#include "nuX_weaklib.hxx"
#include "wli_io_opacity.H"
#include "wli_opacity_emab_iso.H"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
#include <unistd.h>

#include <AMReX_ParallelDescriptor.H>
#include "cctk.h"

namespace nuX_WeakLib::host_io {
namespace {
namespace pd = amrex::ParallelDescriptor;
constexpr const char *explicit_name = "explicit-kernel-units";
constexpr const char *native_name = "weaklib-fullweakphysics-5836de98";
constexpr const char *group = "/Scat_Iso_Kernels";

enum class Fault {
  none, mixed_units, unknown_units, literal_inverse_cm, native_punctuation, missing_units,
  long_explicit_units, overwide_units, nonuniform_energy,
  unit_count, unit_rank, variable_units, native_width, missing_flag0,
  missing_flag1, missing_flag2, flag_value, flag_float, flag_shape,
  ga_nan, ga_shape, ga_float, moments, count_float, flat_offsets,
  nan_offset, species_shape, species_float, energy_count, log_flag_float
};

std::string long_explicit_label() {
  // Exactly 256 non-NUL bytes, with internal rather than trailing whitespace.
  // This exercises the entire accepted fixed-width label plus its 257th-byte
  // audit-record terminator, not merely a label just beyond the old 32 limit.
  return "Per" + std::string(233, ' ') + "Centimeter Per MeV^2";
}

void write_int(H5::H5File &f, const std::string &path, const std::vector<int> &values,
               const std::vector<hsize_t> &dims) {
  const H5::DataSpace space(int(dims.size()), dims.data());
  auto ds = f.createDataSet(path, H5::PredType::NATIVE_INT, space);
  ds.write(values.data(), H5::PredType::NATIVE_INT);
}

void write_double(H5::H5File &f, const std::string &path, const std::vector<double> &values,
                  const std::vector<hsize_t> &dims, bool as_float = false) {
  const H5::DataSpace space(int(dims.size()), dims.data());
  auto ds = f.createDataSet(path, as_float ? H5::PredType::NATIVE_FLOAT :
                                          H5::PredType::NATIVE_DOUBLE, space);
  ds.write(values.data(), H5::PredType::NATIVE_DOUBLE);
}

void write_strings(H5::H5File &f, const std::string &path,
                   const std::vector<std::string> &values,
                   const std::vector<hsize_t> &dims, bool variable = false,
                   std::size_t width = 32) {
  const H5::DataSpace space(int(dims.size()), dims.data());
  const H5::StrType type(H5::PredType::C_S1, variable ? H5T_VARIABLE : width);
  auto ds = f.createDataSet(path, type, space);
  if (variable) {
    std::vector<const char *> pointers;
    for (const auto &value : values) pointers.push_back(value.c_str());
    ds.write(pointers.data(), type);
  } else {
    std::vector<char> data(values.size() * width, '\0');
    for (std::size_t i = 0; i < values.size(); ++i)
      std::copy_n(values[i].data(), std::min(width, values[i].size()), data.data() + i * width);
    ds.write(data.data(), type);
  }
}

void write_fixture(const std::string &path, bool native, Fault fault) {
  // Only paths within our mkdtemp-owned private directory; never truncate a
  // preexisting path. No external table is opened for writing.
  H5::H5File f(path, H5F_ACC_EXCL);
  f.createGroup("/EnergyGrid");
  write_strings(f, "/EnergyGrid/Name", {"Energy"}, {1});
  write_strings(f, "/EnergyGrid/Unit", {"MeV"}, {1});
  write_int(f, "/EnergyGrid/nPoints", {fault == Fault::energy_count ? 4 : 3}, {1});
  write_int(f, "/EnergyGrid/LogInterp", {1}, {1});
  write_double(f, "/EnergyGrid/Values",
               fault == Fault::nonuniform_energy ? std::vector<double>{1, 2, 9} :
                                                  std::vector<double>{1, 2, 4}, {3});
  f.createGroup("/ThermoState");
  write_strings(f, "/ThermoState/Names", {"Density", "Temperature", "Electron Fraction"}, {3});
  write_strings(f, "/ThermoState/Units", {"g/cm3", "K", "None"}, {3});
  if (fault == Fault::log_flag_float)
    write_double(f, "/ThermoState/LogInterp", {1, 1, 0}, {3});
  else
    write_int(f, "/ThermoState/LogInterp", {1, 1, 0}, {3});
  write_double(f, "/ThermoState/Density", {1e8, 1e9}, {2});
  write_double(f, "/ThermoState/Temperature", {1e10, 1e11}, {2});
  write_double(f, "/ThermoState/Electron Fraction", {0.1, 0.4}, {2});
  write_int(f, "/ThermoState/iRho", {1}, {1});
  write_int(f, "/ThermoState/iT", {2}, {1});
  write_int(f, "/ThermoState/iYe", {3}, {1});
  f.createGroup(group);
  const std::string g = group;
  write_int(f, g + "/nOpacities", {2}, {1});
  if (fault == Fault::count_float)
    write_double(f, g + "/nMoments", {2.5}, {1});
  else
    write_int(f, g + "/nMoments", {fault == Fault::moments ? 1 : 2}, {1});
  if (fault != Fault::missing_units) {
    const std::string label = native ? "Per Centimeter" : "Per Centimeter Per MeV^2";
    std::vector<std::string> labels{label, label};
    if (fault == Fault::mixed_units) labels[1] = "Per Centimeter Per MeV^2";
    if (fault == Fault::unknown_units) labels[0] = "undocumented";
    if (fault == Fault::literal_inverse_cm) labels[0] = "1/cm";
    if (fault == Fault::native_punctuation) labels[0] = "Per*(Centimeter)";
    if (fault == Fault::long_explicit_units) labels = {long_explicit_label(), long_explicit_label()};
    if (fault == Fault::unit_count) labels.resize(1);
    write_strings(f, g + "/Units", labels,
                  fault == Fault::unit_rank ? std::vector<hsize_t>{1, 2} :
                      std::vector<hsize_t>{hsize_t(labels.size())},
                  fault == Fault::variable_units,
                  fault == Fault::native_width ? 16 : fault == Fault::long_explicit_units ? 256 :
                      fault == Fault::overwide_units ? 257 : 32);
  }
  std::vector<double> offsets{1, 2, 3, 4}; // s + 2*m, NOT s*2 + m
  if (fault == Fault::nan_offset) offsets[0] = std::numeric_limits<double>::quiet_NaN();
  write_double(f, g + "/Offsets", offsets,
               fault == Fault::flat_offsets ? std::vector<hsize_t>{4} :
                                             std::vector<hsize_t>{2, 2});
  for (int s = 0; s < 2; ++s) {
    std::vector<double> values(48);
    for (std::size_t i = 0; i < values.size(); ++i) {
      const int m = int(i / 3) % 2;
      const double phi = m == 0 ? (s == 0 ? 2 : 5) : (s == 0 ? -0.5 : 1.5);
      // Keep data finite even for the deliberately NaN offset dataset; do not
      // perform log arithmetic on NaN if the application enables FP traps.
      values[i] = std::log10(phi + (1 + s + 2 * m));
    }
    write_double(f, g + (s == 0 ? "/Electron Neutrino" : "/Electron Antineutrino"), values,
                 s == 1 && fault == Fault::species_shape ? std::vector<hsize_t>{2, 2, 2, 3, 2} :
                                                          std::vector<hsize_t>{2, 2, 2, 2, 3},
                 s == 1 && fault == Fault::species_float);
  }
  if (native) {
    const char *flags[] = {"weak_magnetism_corr", "ion_ion_corr", "many_body_corr"};
    for (int i = 0; i < 3; ++i) {
      if ((fault == Fault::missing_flag0 && i == 0) ||
          (fault == Fault::missing_flag1 && i == 1) ||
          (fault == Fault::missing_flag2 && i == 2)) continue;
      const auto p = g + "/" + flags[i];
      if (i == 0 && fault == Fault::flag_float) write_double(f, p, {0.5}, {1});
      else if (i == 0 && fault == Fault::flag_shape) write_int(f, p, {1, 1}, {2});
      else write_int(f, p, {i == 0 && fault == Fault::flag_value ? 2 : (i == 1 ? 0 : 1)}, {1});
    }
    const double ga = fault == Fault::ga_nan ? std::numeric_limits<double>::quiet_NaN() : -0.04;
    write_double(f, g + "/ga_strange", fault == Fault::ga_shape ? std::vector<double>{ga, ga} :
                                                                             std::vector<double>{ga},
                 fault == Fault::ga_shape ? std::vector<hsize_t>{2} : std::vector<hsize_t>{1},
                 fault == Fault::ga_float);
  }
}

// Explicit cleanup is collective-safe; destructor is best-effort on exception.
// Never recursively remove a directory or delete an unowned/external pathname.
struct Scratch {
  std::array<char, 64> directory{};
  bool owned = false;
  std::vector<std::string> files;
  ~Scratch() {
    if (!owned) return;
    for (const auto &path : files) std::remove(path.c_str());
    ::rmdir(directory.data());
  }
  void cleanup() {
    for (const auto &path : files)
      if (std::remove(path.c_str()) && errno != ENOENT)
        throw std::runtime_error("cannot remove owned fixture " + path);
    if (::rmdir(directory.data()))
      throw std::runtime_error("cannot remove owned fixture directory " + std::string(directory.data()));
    files.clear();
    owned = false;
  }
};

bool near(double got, double expected) {
  return std::isfinite(got) && std::abs(got - expected) <= 1e-12 * std::abs(expected);
}

int check_coefficients(const wli::io::HostScatIsoTable &table) {
  if (table.nPoints != std::array<int, 5>{3, 2, 2, 2, 2} ||
      table.nOpacities != 2 || table.nMoments != 2 || table.offset.size() != 4 ||
      table.values.size() != 2 || table.values[0].size() != 48 || table.values[1].size() != 48 ||
      table.common.energyGrid.points.size() != 3)
    return 1;
  int failures = 0;
  const double coordinates[] = {0, 1, 2};
  const double logd[] = {8, 9}, logt[] = {10, 11}, ye[] = {0.1, 0.4};
  TableView v;
  v.ne = 3;
  v.energy_node_coordinates = coordinates;
  v.iso_thermo = {{logd, 2, 1e8, 1e9}, {logt, 2, 1e10, 1e11}, {ye, 2, 0.1, 0.4}};
  for (int s = 0; s < 2; ++s) {
    v.iso[s] = {table.values[s].data(), {3, 2, 2, 2, 2}};
    for (int m = 0; m < 2; ++m)
      v.iso_offsets[s][m] = wli::IsoOffset(table.offset.data(), 2, 2, s, m);
  }
  for (int ie = 0; ie < 3; ++ie) { // exact energy nodes, interior MATTER coordinates
    const double energy = table.common.energyGrid.points[ie];
    for (int s = 0; s < 2; ++s) {
      double phi[2];
      for (int m = 0; m < 2; ++m) {
        phi[m] = iso_at_energy_node(v, s, m, ie, 8.5, 10.5, 0.25);
        const double expected = m == 0 ? (s == 0 ? 2 : 5) : (s == 0 ? -0.5 : 1.5);
        failures += !near(phi[m], expected);
      }
      // Independent angular integral: 2*pi*E^2*int(1-mu)*(phi0+phi1*mu)dmu.
      const double pi = std::acos(-1.0);
      const double expected = 4 * pi * energy * energy * (s == 0 ? 13.0 / 6 : 9.0 / 2);
      failures += !near(4 * pi * energy * energy * (phi[0] - phi[1] / 3), expected);
    }
  }
  return failures;
}
} // namespace

void io_self_test() {
  Scratch scratch;
  collective_root_action("Iso fixture directory creation failed", [&] {
    constexpr char pattern[] = "/tmp/nux-weaklib-iso-XXXXXX";
    std::copy_n(pattern, sizeof(pattern), scratch.directory.begin());
    if (!::mkdtemp(scratch.directory.data()))
      throw std::runtime_error("mkdtemp for Iso fixtures failed");
    scratch.owned = true;
  });
  pd::Bcast(scratch.directory.data(), int(scratch.directory.size()), pd::IOProcessorNumber());
  int sequence = 0, failures = 0;
  const auto fixture = [&](bool native, Fault fault) {
    const auto path = std::string(scratch.directory.data()) + "/case-" + std::to_string(sequence++) + ".h5";
    collective_root_action("Iso fixture creation failed", [&] {
      scratch.files.push_back(path);
      write_fixture(path, native, fault);
    });
    return path;
  };
  const auto rejected = [&](const std::string &path, const std::string &convention,
                             const char *diagnostic) {
    bool failed = false;
    try {
      check_iso_file_collective(path, convention);
    } catch (const std::exception &e) {
      failed = std::string(e.what()).find(diagnostic) != std::string::npos;
    }
    failures += !failed;
  };

  const auto explicit_path = fixture(false, Fault::none);
  const auto native_path = fixture(true, Fault::none);
  const auto long_explicit_path = fixture(false, Fault::long_explicit_units);
  const auto nonuniform_path = fixture(true, Fault::nonuniform_energy);
  const auto explicit_metadata = check_iso_file_collective(explicit_path, explicit_name);
  const auto native_metadata = check_iso_file_collective(native_path, native_name);
  const auto long_metadata = check_iso_file_collective(long_explicit_path, explicit_name);
  check_iso_file_collective(nonuniform_path, native_name);
  const auto long_label = long_explicit_label();
  failures += long_label.size() != 256;
  failures += long_metadata.convention != IsoConvention::explicit_kernel_units;
  for (const auto &unit : long_metadata.units) {
    failures += !std::equal(long_label.begin(), long_label.end(), unit.begin());
    failures += unit[256] != '\0'; // no truncation and room for the terminator
  }
  failures += explicit_metadata.convention != IsoConvention::explicit_kernel_units;
  failures += explicit_metadata.corrections != std::array<int, 3>{-1, -1, -1};
  failures += native_metadata.convention != IsoConvention::fullweakphysics_5836de98;
  failures += native_metadata.corrections != std::array<int, 3>{1, 0, 1};
  failures += native_metadata.ga_strange != -0.04;
  for (const auto &unit : native_metadata.units)
    failures += std::strcmp(unit.data(), "Per Centimeter") != 0;
  // Actual public reader, root-only HDF I/O and collective host-array broadcast.
  const auto explicit_table = wli::io::read_scat_iso_table(explicit_path);
  const auto native_table = wli::io::read_scat_iso_table(native_path);
  const auto long_table = wli::io::read_scat_iso_table(long_explicit_path);
  const auto nonuniform_table = wli::io::read_scat_iso_table(nonuniform_path);
  failures += nonuniform_table.common.energyGrid.points != std::vector<double>{1, 2, 9};
  failures += !physical_energy_nodes_are_valid(nonuniform_table.common.energyGrid.points.data(), 3);
  failures += check_coefficients(nonuniform_table);
  failures += explicit_table.nPoints != native_table.nPoints;
  failures += explicit_table.offset != native_table.offset;
  failures += explicit_table.values != native_table.values; // exact, no format rescaling
  failures += check_coefficients(explicit_table) + check_coefficients(native_table);
  failures += long_table.offset != explicit_table.offset;
  failures += long_table.values != explicit_table.values;
  failures += check_coefficients(long_table);
  rejected(native_path, explicit_name, "iso_table_convention");
  rejected(explicit_path, native_name, "iso_table_convention");
  rejected(native_path, "unknown-convention", "iso_table_convention");
  // An invalid enum must also fail the local API, without accidentally choosing
  // either numeric interpretation. This action/exception is root-only.
  bool invalid_enum_failed = false;
  try {
    collective_root_action("Iso invalid enum test", [&] {
      H5::H5File file(native_path, H5F_ACC_RDONLY);
      validate_iso_metadata(file, static_cast<IsoConvention>(99));
    });
  } catch (const std::exception &e) {
    invalid_enum_failed = std::string(e.what()).find("enum") != std::string::npos;
  }
  failures += !invalid_enum_failed;

  const std::pair<Fault, const char *> malformed[] = {
      {Fault::mixed_units, "/Units"}, {Fault::unknown_units, "/Units"},
      {Fault::native_punctuation, "/Units"},
      {Fault::literal_inverse_cm, "/Units"}, {Fault::missing_units, "/Units"},
      {Fault::unit_count, "/Units"}, {Fault::unit_rank, "/Units"},
      {Fault::variable_units, "/Units"}, {Fault::native_width, "/Units"},
      {Fault::missing_flag0, "weak_magnetism_corr"}, {Fault::missing_flag1, "ion_ion_corr"},
      {Fault::missing_flag2, "many_body_corr"}, {Fault::flag_value, "weak_magnetism_corr"},
      {Fault::flag_float, "weak_magnetism_corr"}, {Fault::flag_shape, "weak_magnetism_corr"},
      {Fault::ga_nan, "ga_strange"}, {Fault::ga_shape, "ga_strange"},
      {Fault::ga_float, "ga_strange"}, {Fault::moments, "/nMoments"},
      {Fault::count_float, "/nMoments"}, {Fault::flat_offsets, "/Offsets"},
      {Fault::nan_offset, "/Offsets"}, {Fault::species_shape, "Electron Antineutrino"},
      {Fault::species_float, "Electron Antineutrino"}, {Fault::energy_count, "/EnergyGrid/nPoints"},
      {Fault::log_flag_float, "/ThermoState/LogInterp"}};
  for (const auto &[fault, diagnostic] : malformed)
    rejected(fixture(true, fault), native_name, diagnostic);
  rejected(fixture(false, Fault::literal_inverse_cm), explicit_name, "/Units");
  rejected(fixture(false, Fault::overwide_units), explicit_name, "/Units");
  rejected(fixture(false, Fault::unknown_units), explicit_name, "/Units");
  // These failures happen only on root, but every rank must catch the status.
  rejected(std::string(scratch.directory.data()) + "/missing.h5", native_name, "validation failed");
  bool creation_failed = false;
  try {
    collective_root_action("expected Iso fixture creation failure", [&] {
      H5::H5File file(std::string(scratch.directory.data()) + "/absent/fail.h5", H5F_ACC_EXCL);
    });
  } catch (const std::exception &e) {
    creation_failed = std::string(e.what()).find("expected Iso fixture creation failure") !=
                      std::string::npos;
  }
  failures += !creation_failed;
  pd::ReduceIntMax(failures);
  pd::Barrier(); // all readers finished before root removes only its own fixtures
  collective_root_action("Iso fixture cleanup failed", [&] { scratch.cleanup(); });
  if (failures)
    throw std::runtime_error("nuX_WeakLib: actual-HDF Iso fixture checks failed: " +
                             std::to_string(failures));
  CCTK_INFO("WeakLib host Iso actual-HDF convention/schema/normalization and "
            "collective-failure fixture checks passed");
}
} // namespace nuX_WeakLib::host_io
