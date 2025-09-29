// nuX_M1_init_fluxes_rhs.cxx — correct, ghost-aware zeroing
#include <loop_device.hxx>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include <vec.hxx>

#define PINDEX1D(ig, iv) ((iv) + (ig) * 5)

namespace nuX_M1 {

using namespace Loop;

CCTK_DEVICE CCTK_HOST inline vect<int, 3> face_centering(int dir) {
  return {dir == 0 ? 0 : 1, dir == 1 ? 0 : 1, dir == 2 ? 0 : 1};
}

CCTK_DEVICE CCTK_HOST inline int face_stride(int dir, const int *lsh) {
  return (dir == 0)   ? (lsh[0] + 1) * lsh[1] * lsh[2]
         : (dir == 1) ? lsh[0] * (lsh[1] + 1) * lsh[2]
                      : lsh[0] * lsh[1] * (lsh[2] + 1);
}

template <int dir> inline void ZeroFaceFlux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_InitFluxesRHS;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_fc(cctkGH, face_centering(dir));

  CCTK_REAL *__restrict__ nu_flux_dir = (dir == 0   ? nu_flux_x
                                         : dir == 1 ? nu_flux_y
                                                    : nu_flux_z);

  const int groupspec = ngroups * nspecies;

  grid.loop_all_device<(dir == 0 ? 0 : 1), (dir == 1 ? 0 : 1),
                       (dir == 2 ? 0 : 1)>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < groupspec; ++ig) {
          for (int iv = 0; iv < 5; ++iv) {
            const int comp = PINDEX1D(ig, iv);
            const ptrdiff_t idx = layout_fc.linear(p.i, p.j, p.k, comp);
            nu_flux_dir[idx] = 0.0;
          }
        }
      });
}

inline void ZeroCellRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_InitFluxesRHS;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const int groupspec = ngroups * nspecies;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rN_rhs[i4D] = 0.0;
          rFx_rhs[i4D] = 0.0;
          rFy_rhs[i4D] = 0.0;
          rFz_rhs[i4D] = 0.0;
          rE_rhs[i4D] = 0.0;
        }
      });
}

extern "C" void nuX_M1_InitFluxesRHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_InitFluxesRHS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_M1_InitFluxesRHS");

  ZeroFaceFlux<0>(cctkGH); // x-faces (Nx+1, Ny,   Nz) + ghosts
  ZeroFaceFlux<1>(cctkGH); // y-faces (Nx,   Ny+1, Nz) + ghosts
  ZeroFaceFlux<2>(cctkGH); // z-faces (Nx,   Ny,   Nz+1) + ghosts
  ZeroCellRHS(cctkGH);
}

} // namespace nuX_M1
