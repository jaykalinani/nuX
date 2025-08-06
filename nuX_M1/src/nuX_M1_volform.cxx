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

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

//#include "nuX_printer.hh" // To be removed


using namespace nuX::m1;

extern "C" void nuX_M1_InitVolform(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTSX_nuX_M1_InitVolform;
    DECLARE_CCTK_PARAMETERS;

    if (verbose) {
        CCTK_INFO("nuX_M1_InitVolform");
    }

    {
        grid.loop_all_device<1, 1, 1>(
            grid.nghostzones,
            [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            // TODO: layout2 dependence
            if (nuX_m1_mask[p.I]) {
                
                // TODO: Unclear whether the following call will work
                volform[p.I] = 0;
                continue;
            }
         const smat<GF3D2<const CCTK_REAL>, dim> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
         const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
            return calc_avg_v2c(gf_g(i, j), p, dir);
         });

         /* determinant of spatial metric */
         const CCTK_REAL detg_avg = calc_det(g_avg);
         volform[p.I] = sqrt(detg_avg);
         assert(isfinite(volform[p.I]));
           
        });
    }

    // Done with printing
   //  nuX::Printer::stop();

}
