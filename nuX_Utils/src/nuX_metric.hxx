#ifndef NUX_METRIC_HXX
#define NUX_METRIC_HXX

#include <cctk.h>

namespace nuX_Utils {
//! Misc utilities for manipulating metric tensors
/*!
 *  NOTE: these are low level routines, kept header-inline and device-safe.
 */
namespace metric {

//! Get normal to space hypersurface: n^μ = (1/α, -β^i/α)
CCTK_HOST CCTK_DEVICE inline
void normal(CCTK_REAL const alp,
            CCTK_REAL const betax,
            CCTK_REAL const betay,
            CCTK_REAL const betaz,
            CCTK_REAL n[4]) {
  // Future-directed unit normal (contravariant)
  n[0] = 1.0 / alp;
  n[1] = -betax / alp;
  n[2] = -betay / alp;
  n[3] = -betaz / alp;
}

//! Construct a spatial metric matrix γ_ij in row-major (xx,xy,xz, yx,yy,yz, zx,zy,zz)
CCTK_HOST CCTK_DEVICE inline
void space(CCTK_REAL const gxx, CCTK_REAL const gxy, CCTK_REAL const gxz,
           CCTK_REAL const gyy, CCTK_REAL const gyz, CCTK_REAL const gzz,
           CCTK_REAL gamma[9]) {
  gamma[0] = gxx; gamma[1] = gxy; gamma[2] = gxz;
  gamma[3] = gxy; gamma[4] = gyy; gamma[5] = gyz;
  gamma[6] = gxz; gamma[7] = gyz; gamma[8] = gzz;
}

//! det(γ_ij) from components
CCTK_HOST CCTK_DEVICE inline
CCTK_REAL spatial_det(CCTK_REAL const gxx, CCTK_REAL const gxy, CCTK_REAL const gxz,
                      CCTK_REAL const gyy, CCTK_REAL const gyz, CCTK_REAL const gzz) {
  // Standard symmetric 3x3 determinant
  return gxx*gyy*gzz + 2.0*gxy*gxz*gyz - gxx*gyz*gyz - gyy*gxz*gxz - gzz*gxy*gxy;
}

//! det(γ_ij) from 3x3 array (row-major)
CCTK_HOST CCTK_DEVICE inline
CCTK_REAL spatial_det(CCTK_REAL const gamma[9]) {
  const CCTK_REAL gxx = gamma[0], gxy = gamma[1], gxz = gamma[2];
  const CCTK_REAL gyy = gamma[4], gyz = gamma[5];
  const CCTK_REAL gzz = gamma[8];
  return spatial_det(gxx,gxy,gxz,gyy,gyz,gzz);
}

//! γ^{ij} via cofactors / det (component form)
CCTK_HOST CCTK_DEVICE inline
void spatial_inv(CCTK_REAL const det,
                 CCTK_REAL const gxx, CCTK_REAL const gxy, CCTK_REAL const gxz,
                 CCTK_REAL const gyy, CCTK_REAL const gyz, CCTK_REAL const gzz,
                 CCTK_REAL *uxx, CCTK_REAL *uxy, CCTK_REAL *uxz,
                 CCTK_REAL *uyy, CCTK_REAL *uyz, CCTK_REAL *uzz) {
  const CCTK_REAL cxx =  gyy*gzz - gyz*gyz;
  const CCTK_REAL cxy = -(gxy*gzz - gxz*gyz);
  const CCTK_REAL cxz =  gxy*gyz - gxz*gyy;
  const CCTK_REAL cyy =  gxx*gzz - gxz*gxz;
  const CCTK_REAL cyz = -(gxx*gyz - gxy*gxz);
  const CCTK_REAL czz =  gxx*gyy - gxy*gxy;

  const CCTK_REAL invdet = 1.0 / det;
  *uxx = cxx * invdet;
  *uxy = cxy * invdet;
  *uxz = cxz * invdet;
  *uyy = cyy * invdet;
  *uyz = cyz * invdet;
  *uzz = czz * invdet;
}

//! γ^{ij} from 3x3 array (row-major) to 3x3 array (row-major)
CCTK_HOST CCTK_DEVICE inline
void spatial_inv(CCTK_REAL const det,
                 CCTK_REAL const gamma[9],
                 CCTK_REAL ugamma[9]) {
  CCTK_REAL uxx,uxy,uxz,uyy,uyz,uzz;
  spatial_inv(det,
              gamma[0],gamma[1],gamma[2],
              gamma[4],gamma[5],gamma[8],
              &uxx,&uxy,&uxz,&uyy,&uyz,&uzz);
  ugamma[0]=uxx; ugamma[1]=uxy; ugamma[2]=uxz;
  ugamma[3]=uxy; ugamma[4]=uyy; ugamma[5]=uyz;
  ugamma[6]=uxz; ugamma[7]=uyz; ugamma[8]=uzz;
}

//! Spacetime metric g_{μν} from ADM (β_i = γ_{ij} β^j)
CCTK_HOST CCTK_DEVICE inline
void spacetime(CCTK_REAL const alp,
               CCTK_REAL const betax, CCTK_REAL const betay, CCTK_REAL const betaz,
               CCTK_REAL const gxx,  CCTK_REAL const gxy,  CCTK_REAL const gxz,
               CCTK_REAL const gyy,  CCTK_REAL const gyz,  CCTK_REAL const gzz,
               CCTK_REAL g[16]) {
  // Build γ_ij
  CCTK_REAL gamma[9];
  space(gxx,gxy,gxz,gyy,gyz,gzz,gamma);

  // β_i = γ_ij β^j
  const CCTK_REAL bix = gamma[0]*betax + gamma[1]*betay + gamma[2]*betaz;
  const CCTK_REAL biy = gamma[3]*betax + gamma[4]*betay + gamma[5]*betaz;
  const CCTK_REAL biz = gamma[6]*betax + gamma[7]*betay + gamma[8]*betaz;

  // β_i β^i
  const CCTK_REAL beta_sq = betax*bix + betay*biy + betaz*biz;

  // Fill g_{μν} (row-major 4x4)
  g[0] = -alp*alp + beta_sq;  g[1] = bix;              g[2] = biy;              g[3] = biz;
  g[4] = bix;                 g[5] = gxx;              g[6] = gxy;              g[7] = gxz;
  g[8] = biy;                 g[9] = gxy;              g[10]= gyy;              g[11]= gyz;
  g[12]= biz;                 g[13]= gxz;              g[14]= gyz;              g[15]= gzz;
}

//! Inverse spacetime metric g^{μν}
CCTK_HOST CCTK_DEVICE inline
void spacetime_upper(CCTK_REAL const alp,
                     CCTK_REAL const betax, CCTK_REAL const betay, CCTK_REAL const betaz,
                     CCTK_REAL const gxx,  CCTK_REAL const gxy,  CCTK_REAL const gxz,
                     CCTK_REAL const gyy,  CCTK_REAL const gyz,  CCTK_REAL const gzz,
                     CCTK_REAL u[16]) {
  // γ^{ij}
  const CCTK_REAL det = spatial_det(gxx,gxy,gxz,gyy,gyz,gzz);
  CCTK_REAL uxx,uxy,uxz,uyy,uyz,uzz;
  spatial_inv(det,gxx,gxy,gxz,gyy,gyz,gzz,&uxx,&uxy,&uxz,&uyy,&uyz,&uzz);

  const CCTK_REAL inv_a2 = 1.0/(alp*alp);

  // g^{00}, g^{0i}
  u[0] = -inv_a2;
  u[1] =  betax*inv_a2;
  u[2] =  betay*inv_a2;
  u[3] =  betaz*inv_a2;

  // g^{i0} = g^{0i}
  u[4] = u[1];
  u[8] = u[2];
  u[12]= u[3];

  // g^{ij} = γ^{ij} - β^i β^j / α^2
  u[5]  = uxx - betax*betax*inv_a2;  u[6]  = uxy - betax*betay*inv_a2;  u[7]  = uxz - betax*betaz*inv_a2;
  u[9]  = uxy - betay*betax*inv_a2;  u[10] = uyy - betay*betay*inv_a2;  u[11] = uyz - betay*betaz*inv_a2;
  u[13] = uxz - betaz*betax*inv_a2;  u[14] = uyz - betaz*betay*inv_a2;  u[15] = uzz - betaz*betaz*inv_a2;
}

} // namespace metric
} // namespace nuX_Utils

#endif

