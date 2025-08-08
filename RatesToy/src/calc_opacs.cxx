#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "setup_eos.hxx"
#include "m1_opacities.hpp"

namespace RatesToy {
using namespace Loop;
using namespace EOSX;
using namespace nuX_Rates;

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

  // Opacity flags
  OpacityFlags opacity_flags = global_opac_flags;

  // Opacity parameters (corrections all switched off)
  OpacityParams opacity_pars = global_opac_params;
    
  // Init EOS
  auto eos_3p = global_eos_3p_tab3d;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

      // Init GreyOpacs struct
      GreyOpacityParams my_grey_opacity_params;

      // Neutrino reactions
      my_grey_opacity_params.opacity_flags = opacity_flags;
      if (debug_ratestoy) { 
        printf("beta = %i\n", my_grey_opacity_params.opacity_flags.use_abs_em);
        printf("pair = %i\n", my_grey_opacity_params.opacity_flags.use_pair);
        printf("brems = %i\n", my_grey_opacity_params.opacity_flags.use_brem);
        printf("inelastic = %i\n", my_grey_opacity_params.opacity_flags.use_inelastic_scatt);
        printf("elastic = %i\n", my_grey_opacity_params.opacity_flags.use_iso);
      }

      // Opacity parameters
      my_grey_opacity_params.opacity_pars = opacity_pars;
      if (debug_ratestoy) { 
        printf("dU = %i\n", my_grey_opacity_params.opacity_pars.use_dU);
        printf("dmeff = %i\n", my_grey_opacity_params.opacity_pars.use_dm_eff);
        printf("WMab = %i\n", my_grey_opacity_params.opacity_pars.use_WM_ab);
        printf("WMsc = %i\n", my_grey_opacity_params.opacity_pars.use_WM_sc);
        printf("decay = %i\n", my_grey_opacity_params.opacity_pars.use_decay);
        printf("BRT brem = %i\n", my_grey_opacity_params.opacity_pars.use_BRT_brem);
        printf("NN medium = %i\n", my_grey_opacity_params.opacity_pars.use_NN_medium_corr);
        printf("neglect block = %i\n", my_grey_opacity_params.opacity_pars.neglect_blocking);
      }

      const int ijk = layout2.linear(p.i, p.j, p.k);

      // Convert Thermodynamic Data to nurates
      CCTK_REAL rhoL = rho[ijk];
      CCTK_REAL tempL = temperature[ijk];
      CCTK_REAL yeL = Ye[ijk];
      my_grey_opacity_params.eos_pars.nb   = rhoL * nuX_dens_conv / (particle_mass * kBS_MeVtog); // CU to nm^-3
      my_grey_opacity_params.eos_pars.temp = tempL;
      my_grey_opacity_params.eos_pars.yp   = yeL;
      my_grey_opacity_params.eos_pars.yn   = 1.0 - yeL;

      CCTK_REAL mu_pL, mu_nL, mu_eL;
      eos_3p->mu_pne_from_valid_rho_temp_ye(rhoL, tempL, yeL, mu_pL, mu_nL, mu_eL);
      my_grey_opacity_params.eos_pars.mu_p = mu_pL;
      my_grey_opacity_params.eos_pars.mu_n = mu_nL;
      my_grey_opacity_params.eos_pars.mu_e = mu_eL;

      // Convert M1 Data to nurates
			for (int ig = 0; ig < ngroups * nspecies; ++ig) {
        const int i4D = layout2.linearVec3(p.i, p.j, p.k, ig);
        CCTK_REAL in_fac = 1.0;
        if (ig == 2)
          CCTK_REAL in_fac = 1.0/4.0; // Heavy neutrinos account for 4 species
        my_grey_opacity_params.m1_pars.n[ig]   = rnt[i4D] * nuX_ndens_conv * in_fac; // fm^-3 to nm^-3
        my_grey_opacity_params.m1_pars.J[ig]   = rJt[i4D] * nuX_edens_conv * in_fac; // CU to MeV nm^-3
        my_grey_opacity_params.m1_pars.chi[ig] = chit[i4D];
      }

      // Distribution parameters
      my_grey_opacity_params.distr_pars = CalculateDistrParamsFromM1(
          &my_grey_opacity_params.m1_pars,
          &my_grey_opacity_params.eos_pars);
          
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
            
      // Convert emissivities, opacities from nurates
			for (int ig = 0; ig < ngroups * nspecies; ++ig) {
        const int i4D = layout2.linearVec3(p.i, p.j, p.k, ig);

        CCTK_REAL out_fac = 1.0;
        if (ig == 2)
          CCTK_REAL out_fac = 4.0; // Heavy neutrinos account for 4 species

        abs_0t[i4D] = coeffs.kappa_0_a[ig] * nuX_length_conv;
        abs_1t[i4D] = coeffs.kappa_a[ig] * nuX_length_conv;
        scat_1t[i4D] = coeffs.kappa_s[ig] * nuX_length_conv;
        eta_0t[i4D] = coeffs.eta_0[ig] / nuX_ndens_conv * nuX_time_conv * out_fac;
        eta_1t[i4D] = coeffs.eta[ig] / nuX_edens_conv * nuX_time_conv * out_fac;
			}
		});
}

} // namespace
