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

extern "C" void nuX_M1_CopyLevels(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_CopyLevels;
  DECLARE_CCTK_PARAMETERS

  if (verbose) {
    CCTK_INFO("nuX_M1_CopyLevels");
  }

  size_t siz = UTILS_GFSIZE(cctkGH) * nspecies * sizeof(CCTK_REAL);

  std::memcpy(rN_p, rN, siz);
  std::memcpy(rE_p, rE, siz);
  std::memcpy(rFx_p, rFx, siz);
  std::memcpy(rFy_p, rFy, siz);
  std::memcpy(rFz_p, rFz, siz);

  if (timelevels > 2) {
    std::memcpy(rN_p_p, rN, siz);
    std::memcpy(rE_p_p, rE, siz);
    std::memcpy(rFx_p_p, rFx, siz);
    std::memcpy(rFy_p_p, rFy, siz);
    std::memcpy(rFz_p_p, rFz, siz);
  }
}
