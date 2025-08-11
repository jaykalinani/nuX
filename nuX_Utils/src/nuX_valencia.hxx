#ifndef NUX_VALENCIA_HXX
#define NUX_VALENCIA_HXX

#include <cctk.h>

namespace nuX_Utils {
namespace valencia {

//! Compute u^μ from (α, β^i, W, v^i) using ADM relations:
//!   n^μ = (1/α, -β^i/α),  u^μ = W ( n^μ + (0, v^i) )
//! ⇒ u^0 = W/α,  u^i = W ( v^i - β^i/α )
CCTK_HOST CCTK_DEVICE inline
void uvel(CCTK_REAL alp,
          CCTK_REAL betax,
          CCTK_REAL betay,
          CCTK_REAL betaz,
          CCTK_REAL w_lorentz,
          CCTK_REAL velx,
          CCTK_REAL vely,
          CCTK_REAL velz,
          CCTK_REAL u[4]) {
  const CCTK_REAL inva = 1.0 / alp;
  u[0] = w_lorentz * inva;
  u[1] = w_lorentz * (velx - betax * inva);
  u[2] = w_lorentz * (vely - betay * inva);
  u[3] = w_lorentz * (velz - betaz * inva);
}

} // namespace valencia
} // namespace nuX_Utils

#endif

