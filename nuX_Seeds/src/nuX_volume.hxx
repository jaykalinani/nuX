#ifndef NUX_VOLUME_HXX
#define NUX_VOLUME_HXX

#include <cctk.h>

namespace nuX_Seeds {
//! Misc utilities for manipulating metric tensors
/*!
 *  NOTE: these are low level routines, kept header-inline and device-safe.
 */
namespace nuX_Seeds_volume {

//! Get normal to space hypersurface: n^μ = (1/α, -β^i/α)
static inline CCTK_HOST CCTK_DEVICE CCTK_REAL volume_f(CCTK_REAL R, CCTK_REAL xp, CCTK_REAL yp,
                             CCTK_REAL zp, CCTK_REAL dx, CCTK_REAL dy,
                             CCTK_REAL dz) {
  const int NPOINTS = 10;
  int inside = 0, count = 0;
  for (int ii = 0; ii < NPOINTS; ++ii) {
    CCTK_REAL const myx = (xp - dx / 2.) + (ii + 0.5) * (dx / NPOINTS);
    for (int jj = 0; jj < NPOINTS; ++jj) {
      CCTK_REAL const myy = (yp - dy / 2.) + (jj + 0.5) * (dy / NPOINTS);
      for (int kk = 0; kk < NPOINTS; ++kk) {
        CCTK_REAL const myz = (zp - dz / 2.) + (kk + 0.5) * (dz / NPOINTS);
        ++count;
        if (myx * myx + myy * myy + myz * myz < R * R)
          ++inside;
      }
    }
  }
  return static_cast<CCTK_REAL>(inside) / static_cast<CCTK_REAL>(count);
}
} // namespace nuX_Seeds_volume
} // namespace nuX_Seeds
#endif