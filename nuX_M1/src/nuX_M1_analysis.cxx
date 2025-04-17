//  Templated Hydrodynamics Code: an hydro code built on top of HRSCCore
//  Copyright (C) 2021, David Radice <david.radice@psu.edu>
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

#include <algorithm>
#include <cassert>
#include <sstream>
#include <cmath>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_M1_macro.hxx"
#include <loop_device.hxx>

#include "utils.hh"

using namespace std;
using namespace Loop;

extern "C" void nuX_M1_Analysis(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_Analysis;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_Analysis");
  }

	// particle_mass is in MeV
	// normfact (?) * cgs2cactusMass * mev_to_erg / (clight*clight)
  CCTK_REAL const mb = (1e50) * (5.0278543128934301e-34) * (1.60217733e-6) * particle_mass / (2.99792458e10 * 2.99792458e10); // AverageBaryonMass();

  assert(nspecies == 3);
  assert(ngroups == 1);

	const GridDescBaseDevice grid(cctkGH);
	const GF3D2layout layout2(cctkGH, {1, 1, 1});

	// UTILS_LOOP3(nuX_m1_analysis, k, 0, cctk_lsh[2], j, 0, cctk_lsh[1], i, 0,
	// 						cctk_lsh[0]) {
	grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

		// int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
		const int ijk = layout2.linear(p.i, p.j, p.k);

		int const ijke = CCTK_VectGFIndex3D(cctkGH, p.i, p.j, p.k, 0);
		int const ijka = CCTK_VectGFIndex3D(cctkGH, p.i, p.j, p.k, 1);
		int const ijkx = CCTK_VectGFIndex3D(cctkGH, p.i, p.j, p.k, 2);

		CCTK_REAL const nb = rho[ijk] / mb;
		ynue[ijk] = rnnu[ijke] / volform[ijk] / nb;
		ynua[ijk] = rnnu[ijka] / volform[ijk] / nb;
		ynux[ijk] = rnnu[ijkx] / volform[ijk] / nb;

		CCTK_REAL const egas = rho[ijk] * (1 + eps[ijk]);
		CCTK_REAL const enue = rJ[ijke] / volform[ijk];
		CCTK_REAL const enua = rJ[ijka] / volform[ijk];
		CCTK_REAL const enux = rJ[ijkx] / volform[ijk];
		CCTK_REAL const etot = egas + enue + enua + enux;
		znue[ijk] = enue / etot;
		znua[ijk] = enua / etot;
		znux[ijk] = enux / etot;

	});
}
