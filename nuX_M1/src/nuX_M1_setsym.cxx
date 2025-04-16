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

#include <cstdio>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "Symmetry.h"

extern "C" void nuX_M1_SetSym(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_SetSym;
  DECLARE_CCTK_PARAMETERS

  char vname[BUFSIZ];

  int const sym_scal[3] = {1, 1, 1};
  int const sym_vecx[3] = {-1, 1, 1};
  int const sym_vecy[3] = {1, -1, 1};
  int const sym_vecz[3] = {1, 1, -1};

  for (int ig = 0; ig < nspecies; ++ig) {
    std::snprintf(vname, BUFSIZ, "nuX_M1::rN[%d]", ig);
    SetCartSymVN(cctkGH, sym_scal, vname);

    std::snprintf(vname, BUFSIZ, "nuX_M1::rE[%d]", ig);
    SetCartSymVN(cctkGH, sym_scal, vname);

    std::snprintf(vname, BUFSIZ, "nuX_M1::rFx[%d]", ig);
    SetCartSymVN(cctkGH, sym_vecx, vname);

    std::snprintf(vname, BUFSIZ, "nuX_M1::rFy[%d]", ig);
    SetCartSymVN(cctkGH, sym_vecy, vname);

    std::snprintf(vname, BUFSIZ, "nuX_M1::rFz[%d]", ig);
    SetCartSymVN(cctkGH, sym_vecz, vname);
  }
}
