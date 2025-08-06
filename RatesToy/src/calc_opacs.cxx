#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "setup_eos.hxx"

#include "m1_opacities.hpp"

namespace RatesToy {
using namespace Loop;
using namespace EOSX;

extern "C" void RatesToy_Calc(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_RatesToy_Calc;
  DECLARE_CCTK_PARAMETERS

  // Set ccc layout
  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  // Init structs for nurates calls
  MyQuadrature my_quad = {.type   = kGauleg,
                          .alpha  = -42.,
                          .dim    = 1,
                          .nx     = 10,
                          .ny     = 1,
                          .nz     = 1,
                          .x1     = 0.,
                          .x2     = 1.,
                          .y1     = -42.,
                          .y2     = -42.,
                          .z1     = -42.,
                          .z2     = -42.,
                          .points = {0},
                          .w      = {0}};
  GaussLegendre(&my_quad);
  printf("Quadratures generated.\n");

  // Opacity flags (activate all reactions)
  OpacityFlags opacity_flags = opacity_flags_default_all;

  // Opacity parameters (corrections all switched off)
  OpacityParams opacity_pars = opacity_params_default_none;
    
  // Init EOS
  auto eos_3p = global_eos_3p_tab3d;

  printf("Starting Loop.\n");
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

      // Init GreyOpacs struct
      printf("Init Gray Opac Params.\n");
      GreyOpacityParams my_grey_opacity_params;

      // my_grey_opacity_params.opacity_flags = {
      //     .use_abs_em          = opacity_flags(0),
      //     .use_pair            = opacity_flags(1),
      //     .use_brem            = opacity_flags(2),
      //     .use_inelastic_scatt = opacity_flags(3),
      //     .use_iso             = opacity_flags(4)};
      my_grey_opacity_params.opacity_flags = opacity_flags;

      // Opacity parameters (corrections all switched off)
      // my_grey_opacity_params.opacity_pars = {
      //     .use_dU             = opacity_pars(0),
      //     .use_dm_eff         = opacity_pars(1),
      //     .use_WM_ab          = opacity_pars(2),
      //     .use_WM_sc          = opacity_pars(3),
      //     .use_decay          = opacity_pars(4),
      //     .use_BRT_brem       = opacity_pars(5),
      //     .use_NN_medium_corr = opacity_pars(6),
      //     .neglect_blocking   = opacity_pars(7)};
      my_grey_opacity_params.opacity_pars = opacity_pars;

      const int ijk = layout2.linear(p.i, p.j, p.k);

      printf("Loading Thermo Data.\n");
      // Load Thermodynamic Data
      CCTK_REAL rhoL = rho[ijk];
      CCTK_REAL tempL = temperature[ijk];
      CCTK_REAL yeL = Ye[ijk];
      my_grey_opacity_params.eos_pars.nb   = rhoL / kBS_Mu * 1e-21;
      my_grey_opacity_params.eos_pars.temp = tempL;
      my_grey_opacity_params.eos_pars.yp   = yeL;
      my_grey_opacity_params.eos_pars.yn   = 1.0 - yeL;
      CCTK_REAL mu_pL, mu_nL, mu_eL;
      eos_3p->mu_pne_from_valid_rho_temp_ye(rhoL, tempL, yeL, &mu_pL, &mu_nL, &mu_eL);
      my_grey_opacity_params.eos_pars.mu_p = mu_pL;
      my_grey_opacity_params.eos_pars.mu_n = mu_nL;
      my_grey_opacity_params.eos_pars.mu_e = mu_eL;

      // Load M1 Data
      printf("Loading M1 Data.\n");
			// for (int ig = 0; ig < ngroups * nspecies; ++ig) {
			for (int ig = 0; ig < 3; ++ig) {
				// int const i4D = CCTK_VectGFIndex3D(cctkGH, p.i, p.j, p.k, ig);
        const int i4D = layout2.linearVec3(p.i, p.j, p.k, ig);
        // printf("At i4D = %i, ijk = %i, point = (%i, %i, %i), group %i.\n", i4D, ijk, p.i, p.j, p.k, ig);
        my_grey_opacity_params.m1_pars.n[ig]   = rnt[i4D];
        my_grey_opacity_params.m1_pars.J[ig]   = rJt[i4D];
        my_grey_opacity_params.m1_pars.H[ig][0]   = rHt_t[i4D];
        my_grey_opacity_params.m1_pars.H[ig][1]   = rHxt[i4D];
        my_grey_opacity_params.m1_pars.H[ig][2]   = rHyt[i4D];
        my_grey_opacity_params.m1_pars.H[ig][3]   = rHzt[i4D];
      }

      // Distribution parameters
      my_grey_opacity_params.distr_pars = CalculateDistrParamsFromM1(
          &my_grey_opacity_params.m1_pars,
          &my_grey_opacity_params.eos_pars);
          
      printf("Done Loading Gray Opac Params.\n");
      // Set up quadrature on GPU
      MyQuadrature gpu_quad;
      gpu_quad.nx = my_quad.nx; 
      for (int idx = 0; idx < gpu_quad.nx; idx++)
      {
          gpu_quad.w[idx]      = my_quad.w[idx];
          gpu_quad.points[idx] = my_quad.points[idx];
      }

      M1Opacities coeffs = ComputeM1Opacities(&gpu_quad, &gpu_quad,
                                              &my_grey_opacity_params);
            
      // Scatter emissivities, opacities
			for (int ig = 0; ig < ngroups * nspecies; ++ig) {
				// int const i4D = CCTK_VectGFIndex3D(cctkGH, p.i, p.j, p.k, ig);
        const int i4D = layout2.linearVec3(p.i, p.j, p.k, ig);
        printf("At idx (%i, %i, %i), rho T Ye = (%e, %e, %e), group %i; abs0[%i] = %e, abs1[%i] = %e\n", 
            p.i, p.j, p.k, rhoL, tempL, yeL, ig, coeffs.kappa_0_a[ig], coeffs.kappa_a[ig]);

        abs_0t[i4D] = coeffs.kappa_0_a[ig];
        abs_1t[i4D] = coeffs.kappa_a[ig];
        eta_0t[i4D] = coeffs.eta_0[ig];
        eta_1t[i4D] = coeffs.eta[ig];
        scat_1t[i4D] = coeffs.kappa_s[ig];
			}
		});
}

} // namespace
