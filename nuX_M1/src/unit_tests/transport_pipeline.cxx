#include <cmath>

#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "test.hxx"

namespace nuX_M1_Tests {

namespace {

using namespace Loop;

#define PINDEX1D(ig, iv) ((iv) + (ig) * 5)

CCTK_HOST CCTK_DEVICE inline CCTK_REAL test_number_density(const int ig) {
  return CCTK_REAL(0.25) + CCTK_REAL(0.01) * ig;
}

CCTK_HOST CCTK_DEVICE inline CCTK_REAL test_energy_density(const int ig) {
  return CCTK_REAL(1.0) + CCTK_REAL(0.05) * ig;
}

bool close_to_zero(const CCTK_REAL value, const CCTK_REAL atol) {
  return std::isfinite(value) && std::fabs(value) <= atol;
}

void require_finite(const CCTK_REAL value, const char *quantity) {
  if (!std::isfinite(value)) {
    CCTK_VERROR("nuX_M1 pipeline test failed: %s is not finite", quantity);
  }
}

void require_zero(const CCTK_REAL value, const CCTK_REAL atol,
                  const char *quantity) {
  if (!close_to_zero(value, atol)) {
    CCTK_VERROR("nuX_M1 pipeline test failed: %s expected 0 got %.17g",
                quantity, static_cast<double>(value));
  }
}

template <typename F>
void for_each_interior(const cGH *cctkGH, F &&f) {
  const GF3D2layout layout(cctkGH, {1, 1, 1});
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    f(layout, p.i, p.j, p.k);
  });
}

} // namespace

} // namespace nuX_M1_Tests

extern "C" void nuX_M1_TestPipelineSetup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GridDescBaseDevice grid(cctkGH);
  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});
  const int groupspec = ngroups * nspecies;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout.linear(p.i, p.j, p.k);
        nuX_m1_mask[ijk] = 0.0;

        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout.linear(p.i, p.j, p.k, ig);
          const CCTK_REAL n = nuX_M1_Tests::test_number_density(ig);
          const CCTK_REAL e = nuX_M1_Tests::test_energy_density(ig);

          rN[i4D] = n;
          rE[i4D] = e;
          rFx[i4D] = 0.0;
          rFy[i4D] = 0.0;
          rFz[i4D] = 0.0;

          rN_p[i4D] = n;
          rE_p[i4D] = e;
          rFx_p[i4D] = 0.0;
          rFy_p[i4D] = 0.0;
          rFz_p[i4D] = 0.0;
        }
      });
}

extern "C" void nuX_M1_TestPipelineZeroOpacity(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GridDescBaseDevice grid(cctkGH);
  const Loop::GF3D2layout layout(cctkGH, {1, 1, 1});
  const int groupspec = ngroups * nspecies;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < groupspec; ++ig) {
          const int i4D = layout.linear(p.i, p.j, p.k, ig);
          abs_0[i4D] = 0.0;
          abs_1[i4D] = 0.0;
          eta_0[i4D] = 0.0;
          eta_1[i4D] = 0.0;
          scat_1[i4D] = 0.0;
          nueave[i4D] = 0.0;
        }
      });
}

extern "C" void nuX_M1_TestPipelineAfterCalcFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace nuX_M1_Tests;
  const CCTK_REAL atol = unit_test_pipeline_atol;
  const int groupspec = ngroups * nspecies;

  for_each_interior(cctkGH,
                    [&](const Loop::GF3D2layout &layout, const int i,
                        const int j, const int k) {
    for (int ig = 0; ig < groupspec; ++ig) {
      const int idx = layout.linear(i, j, k, ig);
      require_close(nu_cmax_x[idx], CCTK_REAL(1), atol,
                    "nu_cmax_x for uniform Minkowski state");
      require_close(nu_cmax_y[idx], CCTK_REAL(1), atol,
                    "nu_cmax_y for uniform Minkowski state");
      require_close(nu_cmax_z[idx], CCTK_REAL(1), atol,
                    "nu_cmax_z for uniform Minkowski state");

      const int comp_n = PINDEX1D(ig, 0);
      const int comp_e = PINDEX1D(ig, 4);
      require_zero(nu_flux_x[layout.linear(i, j, k, comp_n)], atol,
                   "x number flux for zero-flux uniform state");
      require_zero(nu_flux_y[layout.linear(i, j, k, comp_n)], atol,
                   "y number flux for zero-flux uniform state");
      require_zero(nu_flux_z[layout.linear(i, j, k, comp_n)], atol,
                   "z number flux for zero-flux uniform state");
      require_zero(nu_flux_x[layout.linear(i, j, k, comp_e)], atol,
                   "x energy flux for zero-flux uniform state");
      require_zero(nu_flux_y[layout.linear(i, j, k, comp_e)], atol,
                   "y energy flux for zero-flux uniform state");
      require_zero(nu_flux_z[layout.linear(i, j, k, comp_e)], atol,
                   "z energy flux for zero-flux uniform state");

      for (int iv = 0; iv < 5; ++iv) {
        const int comp = PINDEX1D(ig, iv);
        require_finite(nu_flux_x[layout.linear(i, j, k, comp)], "nu_flux_x");
        require_finite(nu_flux_y[layout.linear(i, j, k, comp)], "nu_flux_y");
        require_finite(nu_flux_z[layout.linear(i, j, k, comp)], "nu_flux_z");
      }
    }
                    });
}

extern "C" void nuX_M1_TestPipelineAfterUpdateRHSFromFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace nuX_M1_Tests;
  const CCTK_REAL atol = unit_test_pipeline_atol;
  const int groupspec = ngroups * nspecies;

  for_each_interior(cctkGH,
                    [&](const Loop::GF3D2layout &layout, const int i,
                        const int j, const int k) {
    for (int ig = 0; ig < groupspec; ++ig) {
      const int idx = layout.linear(i, j, k, ig);
      require_zero(rN_rhs[idx], atol,
                   "rN_rhs after uniform UpdateRHSFromFluxes");
      require_zero(rE_rhs[idx], atol,
                   "rE_rhs after uniform UpdateRHSFromFluxes");
      require_zero(rFx_rhs[idx], atol,
                   "rFx_rhs after uniform UpdateRHSFromFluxes");
      require_zero(rFy_rhs[idx], atol,
                   "rFy_rhs after uniform UpdateRHSFromFluxes");
      require_zero(rFz_rhs[idx], atol,
                   "rFz_rhs after uniform UpdateRHSFromFluxes");
    }
                    });
}

extern "C" void nuX_M1_TestPipelineAfterCalcUpdate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace nuX_M1_Tests;
  const CCTK_REAL atol = unit_test_pipeline_atol;
  const int groupspec = ngroups * nspecies;

  for_each_interior(cctkGH,
                    [&](const Loop::GF3D2layout &layout, const int i,
                        const int j, const int k) {
    for (int ig = 0; ig < groupspec; ++ig) {
      const int idx = layout.linear(i, j, k, ig);
      require_close(rN[idx], test_number_density(ig), atol,
                    "rN after zero-source CalcUpdate");
      require_close(rE[idx], test_energy_density(ig), atol,
                    "rE after zero-source CalcUpdate");
      require_zero(rFx[idx], atol, "rFx after zero-source CalcUpdate");
      require_zero(rFy[idx], atol, "rFy after zero-source CalcUpdate");
      require_zero(rFz[idx], atol, "rFz after zero-source CalcUpdate");

      require_zero(source_update_delta_N[idx], atol,
                   "source_update_delta_N for zero source");
      require_zero(source_update_delta_E[idx], atol,
                   "source_update_delta_E for zero source");
      require_zero(source_update_delta_Fx[idx], atol,
                   "source_update_delta_Fx for zero source");
      require_zero(source_update_delta_Fy[idx], atol,
                   "source_update_delta_Fy for zero source");
      require_zero(source_update_delta_Fz[idx], atol,
                   "source_update_delta_Fz for zero source");

      require_close(source_update_status[idx], CCTK_REAL(0), atol,
                    "source_update_status for zero source");
      require_close(closure_update_status[idx], CCTK_REAL(0), atol,
                    "closure_update_status for zero source");
    }
                    });
}
