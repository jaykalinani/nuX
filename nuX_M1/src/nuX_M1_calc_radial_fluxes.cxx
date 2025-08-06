//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2020, David Radice <david.radice@psu.edu>
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <cassert>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <loop_device.hxx>
#include "nuX_utils.hxx"
#include "nuX_M1_closure.hxx"

using namespace utils;
using namespace Loop;
using namespace std;

extern "C" void nuX_M1_CalcRadialFluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CalcRadialFluxes;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_CalcRadialFluxes");
  }

	const GridDescBaseDevice grid(cctkGH);
	const GF3D2layout layout2(cctkGH, {1, 1, 1});

  tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
                                      gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz,
                                      kzz, volform);
  tensor::fluid_velocity_field_const fidu(alp, betax, betay, betaz,
                                          fidu_w_lorentz, fidu_velx, fidu_vely,
                                          fidu_velz);

	// UTILS_LOOP3(thc_m1_calc_radial_fluxes, k, 0, cctk_lsh[2], j, 0, cctk_lsh[1],
	// 						i, 0, cctk_lsh[0]) {
	grid.loop_all_device<1, 1, 1>(
		grid.nghostzones,
		[=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

		// int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
		const int ijk = layout2.linear(p.i, p.j, p.k);

		if (nuX_m1_mask[ijk]) {
			for (int ig = 0; ig < nspecies * ngroups; ++ig) {
			        int const i4D = layout2.linear(p.i, p.j, p.k, ig);
				radial_flux_0[i4D] = 0;
				radial_flux_1[i4D] = 0;
			}
		}
		else {

			tensor::metric<4> g_dd;
			geom.get_metric(ijk, &g_dd);

			tensor::inv_metric<4> g_uu;
			geom.get_inv_metric(ijk, &g_uu);

			tensor::generic<CCTK_REAL, 4, 1> beta_u;
			geom.get_shift_vec(ijk, &beta_u);

			tensor::generic<CCTK_REAL, 4, 1> r_d;
			r_d(0) = 0.0;
			r_d(1) = p.x; // x[ijk];
			r_d(2) = p.y; // y[ijk];
			r_d(3) = p.z; // z[ijk];
			CCTK_REAL const rr = sqrt(tensor::dot(g_uu, r_d, r_d));
			CCTK_REAL const irr = 1.0 / rr;

			tensor::generic<CCTK_REAL, 4, 1> u_u;
			fidu.get(ijk, &u_u);

			tensor::generic<CCTK_REAL, 4, 1> F_d;
			tensor::generic<CCTK_REAL, 4, 1> F_u;
			tensor::generic<CCTK_REAL, 4, 1> H_d;
			tensor::generic<CCTK_REAL, 4, 1> H_u;
			tensor::generic<CCTK_REAL, 4, 1> fnu_u;

			for (int ig = 0; ig < ngroups * nspecies; ++ig) {
				int const i4D = layout2.linear(p.i, p.j, p.k, ig);
				pack_F_d(betax[ijk], betay[ijk], betaz[ijk], rFx[i4D], rFy[i4D],
								 rFz[i4D], &F_d);
				pack_H_d(rHt[i4D], rHx[i4D], rHy[i4D], rHz[i4D], &H_d);

				tensor::contract(g_uu, F_d, &F_u);
				tensor::contract(g_uu, H_d, &H_u);

				radial_flux_1[i4D] = 0.0;
				for (int a = 1; a < 4; ++a) {
					radial_flux_1[i4D] +=
							r_d(a) * irr * calc_E_flux(alp[ijk], beta_u, rE[i4D], F_u, a);
				}

				assemble_fnu(u_u, rJ[i4D], H_u, &fnu_u);
				CCTK_REAL const Gamma = alp[ijk] * fnu_u(0);
				// Note that nnu is densitized here
				CCTK_REAL const nnu = rN[i4D] / Gamma;

				radial_flux_0[i4D] = alp[ijk] * irr * nnu * tensor::dot(fnu_u, r_d);
			}

		} // else
	}); // grid loop
}
