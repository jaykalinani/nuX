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

#include "finite_difference.h"
#include "nuX_M1_closure.hxx"
#include "nuX_M1_macro.hxx"
#include "utils.hh"

#define FDORDER 2

using namespace utils;

// NOTE. Unlike the rest of the M1 code, here spatial indices go over 0,1,2
extern "C" void nuX_M1_CalcGRSources(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    if (verbose) {
        CCTK_INFO("nuX_M1_CalcGRSources");
    }

    // Grid data
    CCTK_REAL const idelta[3] = {
        1.0/CCTK_DELTA_SPACE(0),
        1.0/CCTK_DELTA_SPACE(1),
        1.0/CCTK_DELTA_SPACE(2),
    };

    // Slicing geometry
    utils::tensor::slicing_geometry_const geom(alp, betax, betay, betaz, gxx, gxy, gxz,
            gyy, gyz, gzz, kxx, kxy, kxz, kyy, kyz, kzz, volform);

#pragma omp parallel
    {
        UTILS_LOOP3(nuX_m1_grsource,
                k, NUX_M1_NGHOST, cctk_lsh[2]-NUX_M1_NGHOST,
                j, NUX_M1_NGHOST, cctk_lsh[1]-NUX_M1_NGHOST,
                i, NUX_M1_NGHOST, cctk_lsh[0]-NUX_M1_NGHOST) {
            int const ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);

            // Skip masked points
            if (nuX_m1_mask[ijk]) {
                continue;
            }

            // Geometry on the slice
            utils::tensor::inv_metric<3> g_uu;
            geom.get_inv_metric(ijk, &g_uu);
            utils::tensor::symmetric2<CCTK_REAL, 3, 2> K_dd;
            geom.get_extr_curv(ijk, &K_dd);

            // Lapse derivatives
            utils::tensor::generic<CCTK_REAL, 3, 1> dalp_d;
            for (int a = 0; a < 3; ++a) {
                dalp_d(a) = idelta[a] * cdiff_1(cctkGH, alp, i, j, k, a, FDORDER);
            }
            // Shift derivatives
            utils::tensor::generic<CCTK_REAL, 3, 2> dbeta_du;
            for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                CCTK_REAL const * betab = geom.get_shift_comp(b);
                dbeta_du(a,b) = idelta[a] *
                    cdiff_1(cctkGH, betab, i, j, k, a, FDORDER);
            }
            // Metric derivatives
            utils::tensor::symmetric2<CCTK_REAL, 3, 3> dg_ddd;
            for(int a = 0; a < 3; ++a)
            for(int b = 0; b < 3; ++b)
            for(int c = b; c < 3; ++c) {
                CCTK_REAL const * gbc = geom.get_space_metric_comp(b, c);
                dg_ddd(a,b,c) = idelta[a] *
                    cdiff_1(cctkGH, gbc, i, j, k, a, FDORDER);
            }

            for (int ig = 0; ig < ngroups*nspecies; ++ig) {
                int const i4D = CCTK_VectGFIndex3D(cctkGH, i, j, k, ig);

                // Radiation quantities
                tensor::generic<CCTK_REAL, 3, 1> F_d;
                pack_F_d(rFx[i4D], rFy[i4D], rFz[i4D], &F_d);
                tensor::symmetric2<CCTK_REAL, 3, 2> P_dd;
                pack_P_dd(rPxx[i4D], rPxy[i4D], rPxz[i4D],
                          rPyy[i4D], rPyz[i4D], rPzz[i4D], &P_dd);

                // Contravariant radiation pressure
                tensor::symmetric2<CCTK_REAL, 3, 2> P_uu;
                tensor::contract2(g_uu, P_dd, &P_uu);

                // Radiation flux sources
                tensor::generic<CCTK_REAL *, 3, 1> rF_rhs;
                rF_rhs(0) = &rFx_rhs[i4D];
                rF_rhs(1) = &rFy_rhs[i4D];
                rF_rhs(2) = &rFz_rhs[i4D];

                // Compute radiation energy sources
                // Note that everything is already densitized
                rE_rhs[i4D] += alp[ijk]*tensor::dot(P_uu, K_dd) -
                    tensor::dot(g_uu, F_d, dalp_d);

                // Compute the radiation flux sources
                for (int a = 0; a < 3; ++a) {
                    *rF_rhs(a) -= rE[i4D]*dalp_d(a);
                    for (int b = 0; b < 3; ++b) {
                        *rF_rhs(a) += F_d(b)*dbeta_du(a,b);
                    }
                    for (int b = 0; b < 3; ++b)
                    for (int c = 0; c < 3; ++c) {
                        *rF_rhs(a) += alp[ijk]/2. * P_uu(b,c) * dg_ddd(a,b,c);
                    }
                }
            }

        } UTILS_ENDLOOP3(nuX_m1_grsource);
    }
}
