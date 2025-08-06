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

#include <cstring>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_utils.hxx"
#include <loop_device.hxx>

using namespace Loop;

extern "C" void nuX_M1_Reset(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Reset;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_Reset");
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout2(cctkGH, {1, 1, 1});

  // size_t siz = UTILS_GFSIZE(cctkGH) * nspecies * sizeof(CCTK_REAL);
  // size_t isiz = UTILS_GFSIZE(cctkGH) * sizeof(CCTK_INT);

  // UTILS_LOOP3(nuX_m1_analysis, k, 0, cctk_lsh[2], j, 0, cctk_lsh[1], i, 0,
  //            cctk_lsh[0]) {
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

			for (int ig = 0; ig < ngroups * nspecies; ++ig) {
			        int const i4D = layout2.linear(p.i, p.j, p.k, ig);
				rE[i4D] = rad_E_floor;
				rN[i4D] = rad_N_floor;
				rFx[i4D] = 0.0;
				rFy[i4D] = 0.0;
				rFz[i4D] = 0.0;
				nuX_m1_mask[i4D] = 0.0;
			}
		});
}
