#ifndef NUX_M1_DIAGNOSTICS_HXX
#define NUX_M1_DIAGNOSTICS_HXX

#include "cctk.h"

namespace nuX_M1 {

// Empty specializations keep diagnostic grid-function pointers out of the
// production device kernels.  Besides avoiding unnecessary kernel arguments,
// this sidesteps NVCC's restriction on first captures inside if constexpr in
// extended device lambdas.
template <bool store> struct EquilibriumDiagnostics {
  EquilibriumDiagnostics(CCTK_REAL *, CCTK_REAL *, CCTK_REAL *) {}

  CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE void clear(const int) const {}
  CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE void invalid(const int) const {}

  template <typename Layout>
  CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE void
  record_all(const Layout &, const int, const int, const int, const int,
             const int, const CCTK_REAL, const CCTK_REAL) const {}
};

template <> struct EquilibriumDiagnostics<true> {
  CCTK_REAL *const status;
  CCTK_REAL *const lepton_residual;
  CCTK_REAL *const energy_residual;

  EquilibriumDiagnostics(CCTK_REAL *const status,
                         CCTK_REAL *const lepton_residual,
                         CCTK_REAL *const energy_residual)
      : status(status), lepton_residual(lepton_residual),
        energy_residual(energy_residual) {}

  CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE void clear(const int index) const {
    status[index] = 0.0;
    lepton_residual[index] = 0.0;
    energy_residual[index] = 0.0;
  }

  CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE void invalid(const int index) const {
    status[index] = -2.0;
  }

  template <typename Layout>
  CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE void
  record_all(const Layout &layout, const int i, const int j, const int k,
             const int ng, const int equilibrium_status,
             const CCTK_REAL lepton_error, const CCTK_REAL energy_error) const {
    for (int ig = 0; ig < ng; ++ig) {
      const int index = layout.linear(i, j, k, ig);
      status[index] = equilibrium_status;
      lepton_residual[index] = lepton_error;
      energy_residual[index] = energy_error;
    }
  }
};

} // namespace nuX_M1

#endif
