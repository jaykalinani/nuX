#ifndef NUX_M1_COLLISION_PROVIDER_HXX
#define NUX_M1_COLLISION_PROVIDER_HXX

#include "cctk.h"
#include "cctk_Functions.h"

namespace nuX_M1 {
// Host query once per launch, never once per cell. Optional provider absence
// leaves every existing rate provider's floor-aware path unchanged.
inline bool coupled_collision_active() {
  return CCTK_IsFunctionAliased("nuX_M1_CoupledCollisionActive") &&
         nuX_M1_CoupledCollisionActive() != 0;
}

// Host-call interface for optional coupled collision providers. Existing stage
// workspace is passed explicitly; no extra GF or provider-dependent M1 build
// dependency is needed. Device kernels capture the pointer values, never this
// stack object's address. The caller's schedule owns all read/write metadata.
struct collision_update_buffers_t {
  CCTK_REAL *N, *E, *Fx, *Fy, *Fz, *chi;
  CCTK_REAL *delta_N, *delta_E, *delta_Fx, *delta_Fy, *delta_Fz;
};
} // namespace nuX_M1

#endif
