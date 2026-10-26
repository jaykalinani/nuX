#include <cassert>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "nuX_M1_closure.hxx"
#include "nuX_baryon_mass.hxx"

namespace nuX_M1 {

// RK42-IMEX2L keeps the solved state, then recovers
// g = (Z - base) / (gamma*dt). The base is copied after this PostStep, so a
// repair done here is part of the state and not of g. A zero-rate solve then
// leaves g = 0. IMEX42L rebuilds the step from g and includes a floor in that
// recovered tendency. RK42-IMEX does not call this PostStep before its solve.
extern "C" void nuX_M1_ProjectBeforeImplicit(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_ProjectBeforeImplicit;
  DECLARE_CCTK_PARAMETERS;

  if (!ODESolvers_InPreImplicitStage())
    return;

  // method is the shared ODESolvers::method parameter.
  if (!CCTK_Equals(method, "RK42-IMEX2L"))
    return;

  if (verbose)
    CCTK_INFO("nuX_M1_ProjectBeforeImplicit");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);
  const CCTK_REAL mb = nuX_Utils::AverageBaryonMass(particle_mass);
  const bool couple = backreact;

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        if (nuX_m1_mask[ijk])
          return;

        tensor::generic<CCTK_REAL, 4, 1> beta_u;
        geom.get_shift_vec(p, &beta_u);
        const CCTK_REAL betax_ijk = beta_u(1);
        const CCTK_REAL betay_ijk = beta_u(2);
        const CCTK_REAL betaz_ijk = beta_u(3);
        tensor::inv_metric<4> g_uu;
        geom.get_inv_metric(p, &g_uu);

        const int groupspec = ngroups * nspecies;
        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          const CCTK_REAL Eold = rE[i4D];
          const CCTK_REAL Fxold = rFx[i4D];
          const CCTK_REAL Fyold = rFy[i4D];
          const CCTK_REAL Fzold = rFz[i4D];
          const CCTK_REAL Nold = rN[i4D];

          CCTK_REAL E = Eold;
          tensor::generic<CCTK_REAL, 4, 1> F_d;
          pack_F_d(betax_ijk, betay_ijk, betaz_ijk, Fxold, Fyold, Fzold, &F_d);
          apply_floor(g_uu, &E, &F_d, rad_E_floor, rad_eps);
          const CCTK_REAL N = std::max(Nold, rad_N_floor);

          // The implicit solve applies this same repair to its result. Coupling
          // the repair here, and leaving a zero increment there when the rates
          // are zero, counts the hydro update once.
          if (couple) {
            const CCTK_REAL delta_E = E - Eold;
            const CCTK_REAL delta_Fx = F_d(1) - Fxold;
            const CCTK_REAL delta_Fy = F_d(2) - Fyold;
            const CCTK_REAL delta_Fz = F_d(3) - Fzold;
            const CCTK_REAL delta_N = N - Nold;
            const CCTK_REAL delta_DYe =
                -mb * ((ig == 0 ? delta_N : 0.0) - (ig == 1 ? delta_N : 0.0));
            momx[ijk] -= delta_Fx;
            momy[ijk] -= delta_Fy;
            momz[ijk] -= delta_Fz;
            tau[ijk] -= delta_E;
            DYe[ijk] += delta_DYe;
            // netabs and netheat are analysis-only and are not valid before
            // the first implicit solve, so this repair does not update them.
            assert(isfinite(momx[ijk]));
            assert(isfinite(momy[ijk]));
            assert(isfinite(momz[ijk]));
            assert(isfinite(tau[ijk]));
            assert(isfinite(DYe[ijk]));
          }

          rE[i4D] = E;
          unpack_F_d(F_d, &rFx[i4D], &rFy[i4D], &rFz[i4D]);
          rN[i4D] = N;
          assert(isfinite(rE[i4D]));
          assert(isfinite(rN[i4D]));
          assert(isfinite(rFx[i4D]));
          assert(isfinite(rFy[i4D]));
          assert(isfinite(rFz[i4D]));
        }
      });
}

} // namespace nuX_M1
