#include "nuX_weaklib_io.hxx"

#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include <AMReX_ParallelDescriptor.H>

namespace nuX_WeakLib::host_io {
namespace {
void require(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error("nuX_WeakLib: " + message);
}

H5::DataSet dataset(const H5::H5File &file, const std::string &path) {
  require(file.nameExists(path), "missing dataset " + path);
  return file.openDataSet(path);
}

std::vector<hsize_t> shape(const H5::DataSet &ds, const std::string &path) {
  const auto space = ds.getSpace();
  const int rank = space.getSimpleExtentNdims();
  require(rank >= 0 && rank <= 5, "unsupported HDF dataset rank for " + path);
  std::vector<hsize_t> dims(rank);
  if (rank) space.getSimpleExtentDims(dims.data());
  return dims;
}

void require_shape(const H5::DataSet &ds, const std::string &path,
                   std::initializer_list<hsize_t> expected) {
  require(shape(ds, path) == std::vector<hsize_t>(expected), "wrong shape for " + path);
}

void require_double(const H5::DataSet &ds, const std::string &path) {
  require(ds.getTypeClass() == H5T_FLOAT && ds.getDataType().getSize() == 8 &&
              ds.getFloatType().getPrecision() == 64,
          path + " must have double-precision HDF floating type");
}

std::vector<std::string> strings(const H5::H5File &file, const std::string &path,
                               hsize_t count, std::size_t required_width = 0) {
  const auto ds = dataset(file, path);
  require_shape(ds, path, {count});
  require(ds.getTypeClass() == H5T_STRING, path + " must have string HDF type");
  const auto type = ds.getStrType();
  require(!type.isVariableStr(), path + " must use WeakLib fixed-width strings");
  const auto width = type.getSize();
  require(width > 0 && width <= 256 && (!required_width || width == required_width),
          "invalid fixed string width for " + path);
  std::vector<char> data(std::size_t(count) * width);
  ds.read(data.data(), type);
  std::vector<std::string> result;
  for (hsize_t i = 0; i < count; ++i) {
    std::string value(data.data() + std::size_t(i) * width, width);
    const auto end = value.find('\0');
    if (end != std::string::npos) value.resize(end);
    const auto last = value.find_last_not_of(' ');
    value.resize(last == std::string::npos ? 0 : last + 1);
    result.push_back(value);
  }
  return result;
}

std::string normalized(std::string value, bool unit_punctuation = true) {
  value.erase(std::remove_if(value.begin(), value.end(), [=](unsigned char c) {
    return std::isspace(c) ||
           (unit_punctuation && (c == '^' || c == '*' || c == '(' || c == ')'));
  }), value.end());
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return char(std::tolower(c));
  });
  return value;
}

bool one_of(const std::string &value, std::initializer_list<const char *> choices) {
  return std::any_of(choices.begin(), choices.end(), [&](const char *x) { return value == x; });
}

std::vector<long long> integers(const H5::H5File &file, const std::string &path,
                                hsize_t count) {
  const auto ds = dataset(file, path);
  require_shape(ds, path, {count});
  require(ds.getTypeClass() == H5T_INTEGER && ds.getDataType().getSize() <= 8,
          path + " must have integer HDF type");
  std::vector<long long> result(count);
  ds.read(result.data(), H5::PredType::NATIVE_LLONG);
  return result;
}

long long integer(const H5::H5File &file, const std::string &path) {
  return integers(file, path, 1)[0];
}

void common_units(const H5::H5File &file) {
  require(normalized(strings(file, "/EnergyGrid/Unit", 1)[0]) == "mev",
          "/EnergyGrid/Unit must be MeV");
  const auto thermo = strings(file, "/ThermoState/Units", 3);
  require(one_of(normalized(thermo[0]), {"g/cm3", "gram/cm3", "gramspercm3", "gcm-3"}) &&
              one_of(normalized(thermo[1]), {"k", "kelvin"}) &&
              one_of(normalized(thermo[2]), {"", "none", "1", "dimensionless"}),
          "/ThermoState/Units must be [g/cm3, K, dimensionless]");
}

hsize_t axis_extent(const H5::H5File &file, const std::string &path) {
  const auto ds = dataset(file, path);
  require_double(ds, path);
  const auto dims = shape(ds, path);
  require(dims.size() == 1 && dims[0] >= 2 && dims[0] <= INT_MAX,
          path + " must be a 1D axis with 2..INT_MAX nodes");
  return dims[0];
}

void standard_thermo_layout(const H5::H5File &file) {
  const auto names = strings(file, "/ThermoState/Names", 3);
  require(names == std::vector<std::string>{"Density", "Temperature", "Electron Fraction"},
          "unsupported /ThermoState/Names or axis order");
  require(integers(file, "/ThermoState/LogInterp", 3) == std::vector<long long>{1, 1, 0},
          "/ThermoState/LogInterp must be [1,1,0]");
  for (int i = 0; i < 3; ++i) {
    const char *roles[] = {"iRho", "iT", "iYe"};
    require(integer(file, std::string("/ThermoState/") + roles[i]) == i + 1,
            std::string("unsupported /ThermoState/") + roles[i]);
  }
}
} // namespace

IsoConvention parse_iso_convention(const std::string &name) {
  if (name == "explicit-kernel-units") return IsoConvention::explicit_kernel_units;
  if (name == "weaklib-fullweakphysics-5836de98") return IsoConvention::fullweakphysics_5836de98;
  throw std::runtime_error("nuX_WeakLib: unknown iso_table_convention '" + name + "'");
}

const char *iso_convention_name(IsoConvention convention) {
  switch (convention) {
  case IsoConvention::explicit_kernel_units: return "explicit-kernel-units";
  case IsoConvention::fullweakphysics_5836de98: return "weaklib-fullweakphysics-5836de98";
  }
  throw std::runtime_error("nuX_WeakLib: invalid Iso convention enum");
}

IsoMetadata validate_iso_metadata(const H5::H5File &file, IsoConvention convention) {
  iso_convention_name(convention); // reject invalid direct-API enum values too
  const bool native = convention == IsoConvention::fullweakphysics_5836de98;
  common_units(file);
  const std::string group = "/Scat_Iso_Kernels";
  IsoMetadata metadata;
  metadata.convention = convention;
  const auto units = strings(file, group + "/Units", 2, native ? 32 : 0);
  for (int s = 0; s < 2; ++s) {
    const auto unit = normalized(units[s], !native);
    require(native ? unit == "percentimeter"
                   : one_of(unit, {"1/cmmev2", "1/cm/mev2", "cm-1mev-2",
                                   "1percmpermev2", "percentimeterpermev2"}),
            group + "/Units: label '" + units[s] + "' contradicts iso_table_convention='" +
                iso_convention_name(convention) + "'; native Per Centimeter requires an explicit "
                "weaklib-fullweakphysics-5836de98 user declaration, not dimensional inference");
    std::copy(units[s].begin(), units[s].end(), metadata.units[s].begin());
  }
  require(integer(file, group + "/nOpacities") == 2, group + "/nOpacities must be 2");
  require(integer(file, group + "/nMoments") == 2, group + "/nMoments must be 2");
  const auto offset = dataset(file, group + "/Offsets");
  require_shape(offset, group + "/Offsets", {2, 2});
  require_double(offset, group + "/Offsets");
  std::array<double, 4> offsets{};
  offset.read(offsets.data(), H5::PredType::NATIVE_DOUBLE);
  require(std::all_of(offsets.begin(), offsets.end(), [](double x) { return std::isfinite(x); }),
          group + "/Offsets contains a non-finite value");

  // Check the actual on-disk coordinate roles/types/extents BEFORE WLI reads
  // arrays. Numeric energy ranges and MATTER-axis uniformity remain checked by
  // the upload path; exact energy-node access does not require uniform log(E).
  const auto energy_name = strings(file, "/EnergyGrid/Name", 1);
  require(!energy_name[0].empty(), "/EnergyGrid/Name must not be empty");
  const auto energy_log = integer(file, "/EnergyGrid/LogInterp");
  require(energy_log == 0 || energy_log == 1, "/EnergyGrid/LogInterp must be 0 or 1");
  const auto ne = axis_extent(file, "/EnergyGrid/Values");
  require(integer(file, "/EnergyGrid/nPoints") == static_cast<long long>(ne),
          "/EnergyGrid/nPoints disagrees with /EnergyGrid/Values");
  standard_thermo_layout(file);
  const auto nd = axis_extent(file, "/ThermoState/Density");
  const auto nt = axis_extent(file, "/ThermoState/Temperature");
  const auto ny = axis_extent(file, "/ThermoState/Electron Fraction");
  const std::vector<hsize_t> expected{ny, nt, nd, 2, ne}; // C-order is reversed Fortran
  std::size_t count = 1;
  for (const auto n : expected) {
    require(n <= INT_MAX / count, group + " exceeds WLI's INT_MAX-element broadcast contract");
    count *= n;
  }
  for (const char *species : {"Electron Neutrino", "Electron Antineutrino"}) {
    const auto path = group + "/" + species;
    const auto ds = dataset(file, path);
    require(shape(ds, path) == expected, "wrong species/axis shape for " + path);
    require_double(ds, path);
  }
  if (native) {
    const char *flags[] = {"weak_magnetism_corr", "ion_ion_corr", "many_body_corr"};
    for (int i = 0; i < 3; ++i) {
      const auto path = group + "/" + flags[i];
      const auto value = integer(file, path);
      require(value == 0 || value == 1, path + " must be 0 or 1");
      metadata.corrections[i] = int(value);
    }
    const auto path = group + "/ga_strange";
    const auto ds = dataset(file, path);
    require_shape(ds, path, {1});
    require_double(ds, path);
    ds.read(&metadata.ga_strange, H5::PredType::NATIVE_DOUBLE);
    require(std::isfinite(metadata.ga_strange), path + " must be finite");
  }
  return metadata;
}

void collective_root_action(const std::string &context,
                            const std::function<void()> &action) {
  const int root = amrex::ParallelDescriptor::IOProcessorNumber();
  int bad = 0;
  std::array<char, 1024> message{};
  if (amrex::ParallelDescriptor::IOProcessor()) {
    try {
      action();
    } catch (const H5::Exception &e) {
      bad = 1;
      std::snprintf(message.data(), message.size(), "%s", e.getDetailMsg().c_str());
    } catch (const std::exception &e) {
      bad = 1;
      std::snprintf(message.data(), message.size(), "%s", e.what());
    } catch (...) {
      bad = 1;
      std::snprintf(message.data(), message.size(), "%s", "unknown host I/O exception");
    }
  }
  amrex::ParallelDescriptor::Bcast(&bad, 1, root);
  amrex::ParallelDescriptor::Bcast(message.data(), int(message.size()), root);
  if (bad) require(false, context + ": " + message.data());
}

IsoMetadata check_iso_file_collective(const std::string &path,
                                     const std::string &convention) {
  IsoMetadata metadata;
  collective_root_action("Iso unit/schema validation failed for " + path, [&] {
    const auto selected = parse_iso_convention(convention);
    H5::H5File file(path, H5F_ACC_RDONLY);
    metadata = validate_iso_metadata(file, selected);
  });
  const int root = amrex::ParallelDescriptor::IOProcessorNumber();
  int selected = int(metadata.convention);
  amrex::ParallelDescriptor::Bcast(&selected, 1, root);
  metadata.convention = IsoConvention(selected);
  for (auto &unit : metadata.units)
    amrex::ParallelDescriptor::Bcast(unit.data(), int(unit.size()), root);
  amrex::ParallelDescriptor::Bcast(metadata.corrections.data(), 3, root);
  amrex::ParallelDescriptor::Bcast(&metadata.ga_strange, 1, root);
  return metadata;
}

void check_emab_file_collective(const std::string &path) {
  collective_root_action("EmAb unit/schema validation failed for " + path, [&] {
    H5::H5File file(path, H5F_ACC_RDONLY);
    common_units(file);
    const std::string group = file.nameExists("/EmAb") ? "/EmAb" : "/EmAb_CorrectedAbsorption";
    for (const auto &unit : strings(file, group + "/Units", 2))
      require(one_of(normalized(unit), {"1/cm", "cm-1", "1percm", "percentimeter"}),
              "unsupported kernel unit '" + unit + "' in " + group);
    std::vector<hsize_t> previous;
    for (const char *name : {"Electron Neutrino", "Electron Antineutrino"}) {
      const auto path = group + "/" + name;
      const auto dims = shape(dataset(file, path), path);
      require(dims.size() == 4, "wrong table rank for " + path);
      require(previous.empty() || previous == dims, "species table shapes disagree in " + group);
      previous = dims;
    }
  });
}

void check_two_energy_file_collective(const std::string &path,
                                      const std::string &channel) {
  collective_root_action(channel + " unit/schema validation failed for " + path, [&] {
    require(channel == "Pair" || channel == "NES" || channel == "Brem",
            "unsupported two-energy channel");
    H5::H5File file(path, H5F_ACC_RDONLY);
    common_units(file);
    standard_thermo_layout(file);
    const bool brem = channel == "Brem";
    const std::string group = "/Scat_" + channel + "_Kernels";
    const auto unit = normalized(strings(file, group + "/Units", 1)[0]);
    require(brem ? one_of(unit, {"1/mev", "mev-1", "permev"}) :
        one_of(unit, {"1/cmmev3", "1/cm/mev3", "cm-1mev-3", "percentimeterpermev3"}),
        "unsupported physical unit in " + group + "/Units");
    const int moments = brem ? 1 : 4;
    require(integer(file, group + "/nOpacities") == 1 &&
                integer(file, group + "/nMoments") == moments,
            "unsupported species/moment layout in " + group);
    const auto ne = axis_extent(file, "/EnergyGrid/Values");
    const auto nt = axis_extent(file, "/ThermoState/Temperature");
    const auto fourth = brem ? axis_extent(file, "/ThermoState/Density") : nt;
    const auto fifth = brem ? nt : axis_extent(file, "/EtaGrid/Values");
    if (!brem) {
      const auto etaunit = normalized(strings(file, "/EtaGrid/Unit", 1)[0]);
      require(one_of(etaunit, {"dimensionless", "none", "1"}),
              "/EtaGrid must be dimensionless electron degeneracy");
      require(integer(file, "/EtaGrid/LogInterp") == 1,
              "/EtaGrid requires logarithmic interpolation");
    }
    if (channel == "NES")
      require(integer(file, group + "/NPS") == 1,
              "NES requires NPS=1 to include both electron and positron scattering");
    const std::string values = group + (brem ? "/S_sigma" : "/Kernels");
    const auto data = dataset(file, values);
    require_shape(data, values, {fifth, fourth, hsize_t(moments), ne, ne});
    require_double(data, values);
    std::size_t count = 1;
    for (const auto extent : shape(data, values)) {
      require(extent <= INT_MAX / count,
              values + " exceeds WLI's INT_MAX-element broadcast contract");
      count *= extent;
    }
    const auto offsets = dataset(file, group + "/Offsets");
    require_shape(offsets, group + "/Offsets", {hsize_t(moments), 1});
    require_double(offsets, group + "/Offsets");
  });
}
} // namespace nuX_WeakLib::host_io
