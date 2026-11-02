#ifndef NUX_WEAKLIB_IO_HXX
#define NUX_WEAKLIB_IO_HXX

// HOST ONLY: never include this header from the device math header or capture
// its metadata in TableView. This adapter does not rescale any table values.
#include <array>
#include <functional>
#include <string>

#include <H5Cpp.h>

namespace nuX_WeakLib::host_io {

enum class IsoConvention { explicit_kernel_units, fullweakphysics_5836de98 };

struct IsoMetadata {
  IsoConvention convention = IsoConvention::explicit_kernel_units;
  // Preserve the old explicit-label gate's <=256-byte fixed-string allowance.
  // The pinned native format separately requires exactly 32-byte strings.
  std::array<std::array<char, 257>, 2> units{};
  std::array<int, 3> corrections{-1, -1, -1}; // native-only scientific metadata
  double ga_strange = 0; // meaningful only for the native convention
};

IsoConvention parse_iso_convention(const std::string &name);
const char *iso_convention_name(IsoConvention convention);
// Root-local, throws with a dataset path. Also usable by host-only tests.
IsoMetadata validate_iso_metadata(const H5::H5File &file, IsoConvention convention);
// Every rank calls these. Only the I/O root opens the file; errors are broadcast
// before any WLI array-load collective. Unknown declarations fail collectively.
IsoMetadata check_iso_file_collective(const std::string &path,
                                     const std::string &convention);
void check_emab_file_collective(const std::string &path);
void check_two_energy_file_collective(const std::string &path,
                                      const std::string &channel);

// Shared host-only handshake, also used for opt-in fixture creation/cleanup.
// action executes on the I/O root ONLY; it must not contain MPI collectives.
void collective_root_action(const std::string &context,
                            const std::function<void()> &action);
void io_self_test();

} // namespace nuX_WeakLib::host_io
#endif
