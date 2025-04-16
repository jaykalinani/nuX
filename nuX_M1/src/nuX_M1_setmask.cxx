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

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M1_macro.h"
#include "utils.hh"
#include <loop_device.hxx>

using namespace Loop;

extern "C" void nuX_M1_SetMask(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_SetMask");
  }

	// UTILS_LOOP3(thc_m1_setmask, k, 0, cctk_lsh[2], j, 0, cctk_lsh[1], i, 0,
	// 						cctk_lsh[0]) {
	const GridDescBaseDevice grid(cctkGH);
	const GF3D2layout layout2(cctkGH, {1, 1, 1});

	grid.loop_all_device<1, 1, 1>(
		grid.nghostzones,
		[=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

		const int ijk = layout2.linear(p.i, p.j, p.k);
		/*
		int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
		if (hydro_excision && hydro_excision_mask[ijk]) {
			thc_m1_mask[ijk] = 1;
			for (int ig = 0; ig < nspecies * ngroups; ++ig) {
				int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);
				rN[i4D] = 0;
				rE[i4D] = 0;
				rFx[i4D] = 0;
				rFy[i4D] = 0;
				rFz[i4D] = 0;
			}
		}
		else {
			thc_m1_mask[ijk] = 0;
		}
		*/
		nuX_m1_mask[ijk] = 0;
	});
}
