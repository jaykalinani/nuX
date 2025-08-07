#ifndef NUX_TENSOR_HXX
#define NUX_TENSOR_HXX

#include <cctk.h>
#include <cstddef>
#include <algorithm>
#include <cmath>

namespace nuX_Utils {
namespace tensor {

//------------------------------------------------------------------------------
// Kronecker delta
//------------------------------------------------------------------------------
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
inline CCTK_REAL delta(int a, int b) { return a == b ? CCTK_REAL(1) : CCTK_REAL(0); }

//------------------------------------------------------------------------------
// generic<T, ndim, rank>: rank-1/2/3 dense tensors (row-major)
//------------------------------------------------------------------------------
template<typename T, int ndim, int rank> class generic;

// rank-1
template<typename T, int ndim_>
class generic<T, ndim_, 1> {
public:
  using data_t = T;
  static int const rank = 1;
  static int const ndim = ndim_;
  static int const ndof = ndim;

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  static int multiindex(int a) { return a; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& operator[](int a) { return m_data[a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& operator[](int a) const { return m_data[a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& operator()(int a) { return m_data[a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& operator()(int a) const { return m_data[a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& at(int a) { return m_data[a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& at(int a) const { return m_data[a]; }

private:
  T m_data[ndim];
};

// rank-2
template<typename T, int ndim_>
class generic<T, ndim_, 2> {
public:
  using data_t = T;
  static int const rank = 2;
  static int const ndim = ndim_;
  static int const ndof = ndim*ndim;

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  static int multiindex(int a, int b) { return a*ndim + b; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& operator[](int a) { return m_data[a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& operator[](int a) const { return m_data[a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& operator()(int a, int b) { return m_data[multiindex(a,b)]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& operator()(int a, int b) const { return m_data[multiindex(a,b)]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& at(int a, int b) { return m_data[multiindex(a,b)]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& at(int a, int b) const { return m_data[multiindex(a,b)]; }

private:
  T m_data[ndof];
};

// rank-3
template<typename T, int ndim_>
class generic<T, ndim_, 3> {
public:
  using data_t = T;
  static int const rank = 3;
  static int const ndim = ndim_;
  static int const ndof = ndim*ndim*ndim;

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  static int multiindex(int a, int b, int c) { return a*ndim*ndim + b*ndim + c; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& operator[](int a) { return m_data[a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& operator[](int a) const { return m_data[a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& operator()(int a, int b, int c) { return m_data[multiindex(a,b,c)]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& operator()(int a, int b, int c) const { return m_data[multiindex(a,b,c)]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& at(int a, int b, int c) { return m_data[multiindex(a,b,c)]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& at(int a, int b, int c) const { return m_data[multiindex(a,b,c)]; }

private:
  T m_data[ndof];
};

//------------------------------------------------------------------------------
// symmetric2<T,ndim,2>: store only i<=j components in row-major packed layout
//------------------------------------------------------------------------------
template<typename T, int ndim, int rank> class symmetric2;

template<typename T, int ndim_>
class symmetric2<T, ndim_, 2> {
public:
  using data_t = T;
  static int const rank = 2;
  static int const ndim = ndim_;
  static int const ndof = (ndim*(ndim+1))/2;

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  static int multiindex(int a, int b) {
    if (b < a) std::swap(a,b);
    int const offset = ndim*a - (a*(a-1))/2;
    return offset + (b - a);
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& operator[](int a) { return m_data[a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& operator[](int a) const { return m_data[a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& operator()(int a, int b) { return m_data[multiindex(a,b)]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& operator()(int a, int b) const { return m_data[multiindex(a,b)]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T& at(int a, int b) { return m_data[multiindex(a,b)]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  T const& at(int a, int b) const { return m_data[multiindex(a,b)]; }

private:
  T m_data[ndof];
};

//------------------------------------------------------------------------------
// Basic contractions / dots (host+device safe, simple loops)
//------------------------------------------------------------------------------
template<typename T, int N>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
void contract(symmetric2<T,N,2> const& g, generic<T,N,1> const& va, generic<T,N,1>* vb) {
  for (int a=0;a<N;++a) {
    (*vb)[a] = T(0);
    for (int b=0;b<N;++b) (*vb)[a] += g(a,b) * va(b);
  }
}

template<typename T, int N>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
void contract(symmetric2<T,N,2> const& g, symmetric2<T,N,2> const& A, generic<T,N,2>* B) {
  for (int a=0;a<N;++a)
  for (int b=0;b<N;++b) {
    (*B)(a,b) = T(0);
    for (int c=0;c<N;++c) (*B)(a,b) += g(a,c) * A(c,b);
  }
}

template<typename T, int N>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
void contract2(symmetric2<T,N,2> const& g, symmetric2<T,N,2> const& A, symmetric2<T,N,2>* B) {
  for (int a=0;a<N;++a)
  for (int b=a;b<N;++b) {
    B->at(a,b) = T(0);
    for (int c=0;c<N;++c)
    for (int d=0;d<N;++d) B->at(a,b) += g(a,c)*g(b,d)*A(c,d);
  }
}

template<typename T, int N>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
T dot(generic<T,N,1> const& va, generic<T,N,1> const& vb) {
  T s(0); for (int a=0;a<N;++a) s += va(a)*vb(a); return s;
}

template<typename T, int N>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
T dot(symmetric2<T,N,2> const& A, symmetric2<T,N,2> const& B) {
  T s(0); for (int a=0;a<N;++a) for (int b=0;b<N;++b) s += A(a,b)*B(a,b); return s;
}

template<typename T, int N>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
T dot(symmetric2<T,N,2> const& g, generic<T,N,1> const& va, generic<T,N,1> const& vb) {
  T s(0); for (int a=0;a<N;++a) for (int b=0;b<N;++b) s += g(a,b)*va(a)*vb(b); return s;
}

template<typename T, int N>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
T dot(symmetric2<T,N,2> const& g, symmetric2<T,N,2> const& A, symmetric2<T,N,2> const& B) {
  T s(0);
  for (int a=0;a<N;++a) for (int b=0;b<N;++b)
  for (int c=0;c<N;++c) for (int d=0;d<N;++d) s += g(a,c)*g(b,d)*A(c,d)*B(a,b);
  return s;
}

//------------------------------------------------------------------------------
// 3x3 symmetric determinant + inverse (ADM spatial metric)
//------------------------------------------------------------------------------
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
inline CCTK_REAL spatial_det(CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz,
                             CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz) {
  // det = gxx*(gyy*gzz - gyz^2) - gxy*(gxy*gzz - gxz*gyz) + gxz*(gxy*gyz - gxz*gyy)
  return gxx*(gyy*gzz - gyz*gyz) - gxy*(gxy*gzz - gxz*gyz) + gxz*(gxy*gyz - gxz*gyy);
}

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
inline void spatial_inv(CCTK_REAL det,
                        CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz,
                        CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz,
                        CCTK_REAL* uxx, CCTK_REAL* uxy, CCTK_REAL* uxz,
                        CCTK_REAL* uyy, CCTK_REAL* uyz, CCTK_REAL* uzz) {
  CCTK_REAL const A11 = gyy*gzz - gyz*gyz;
  CCTK_REAL const A12 = gxz*gyz - gxy*gzz;
  CCTK_REAL const A13 = gxy*gyz - gxz*gyy;
  CCTK_REAL const A22 = gxx*gzz - gxz*gxz;
  CCTK_REAL const A23 = gxy*gxz - gxx*gyz;
  CCTK_REAL const A33 = gxx*gyy - gxy*gxy;
  CCTK_REAL const invdet = CCTK_REAL(1)/det;

  *uxx = A11*invdet;
  *uxy = A12*invdet;
  *uxz = A13*invdet;
  *uyy = A22*invdet;
  *uyz = A23*invdet;
  *uzz = A33*invdet;
}

//------------------------------------------------------------------------------
// metric / inv_metric specializations
//------------------------------------------------------------------------------
template<int ndim> class metric;
template<int ndim> class inv_metric;

// Spatial metric (lower)
template<>
class metric<3> : public symmetric2<CCTK_REAL,3,2> {
public:
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  CCTK_REAL det() const {
    return spatial_det(this->operator, this->operator, this->operator,
                       this->operator, this->operator, this->operator);
  }
};

// Spatial inverse metric (upper)
template<>
class inv_metric<3> : public symmetric2<CCTK_REAL,3,2> {
public:
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void from_metric(CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz,
                   CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz) {
    CCTK_REAL const det = spatial_det(gxx,gxy,gxz,gyy,gyz,gzz);
    CCTK_REAL uxx,uxy,uxz,uyy,uyz,uzz;
    spatial_inv(det,gxx,gxy,gxz,gyy,gyz,gzz,&uxx,&uxy,&uxz,&uyy,&uyz,&uzz);
    this->operator=uxx; this->operator=uxy; this->operator=uxz;
    this->operator=uyy; this->operator=uyz; this->operator=uzz;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void from_metric(CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz,
                   CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz,
                   CCTK_REAL det) {
    CCTK_REAL uxx,uxy,uxz,uyy,uyz,uzz;
    spatial_inv(det,gxx,gxy,gxz,gyy,gyz,gzz,&uxx,&uxy,&uxz,&uyy,&uyz,&uzz);
    this->operator=uxx; this->operator=uxy; this->operator=uxz;
    this->operator=uyy; this->operator=uyz; this->operator=uzz;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void from_metric(metric<3> const& g) {
    from_metric(g[0],g[1],g[2],g[3],g[4],g[5]);
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void from_metric_det(metric<3> const& g, CCTK_REAL det) {
    from_metric(g[0],g[1],g[2],g[3],g[4],g[5],det);
  }
};

// Spacetime metric (lower) from ADM (for completeness, you can skip if unused)
template<>
class metric<4> : public symmetric2<CCTK_REAL,4,2> {
public:
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void from_adm(CCTK_REAL alp, CCTK_REAL betax, CCTK_REAL betay, CCTK_REAL betaz,
                CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz,
                CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz) {
    // Build full g_ab (lower) using ADM: g_00 = -alpha^2 + beta_i beta^i; g_0i = beta_i; g_ij = gamma_ij
    // Here we need beta_i = gamma_ij beta^j
    inv_metric<3> gamma_uu; gamma_uu.from_metric(gxx,gxy,gxz,gyy,gyz,gzz); // for convenience if needed later
    // Compute beta_i
    CCTK_REAL beta_u[3] = {betax,betay,betaz};
    CCTK_REAL beta_d[3];
    // beta_i = gamma_ij beta^j
    beta_d[0] = gxx*beta_u[0] + gxy*beta_u[1] + gxz*beta_u[2];
    beta_d[1] = gxy*beta_u[0] + gyy*beta_u[1] + gyz*beta_u[2];
    beta_d[2] = gxz*beta_u[0] + gyz*beta_u[1] + gzz*beta_u[2];
    CCTK_REAL beta2 = beta_u[0]*beta_d[0] + beta_u[1]*beta_d[1] + beta_u[2]*beta_d[2];

    // Fill packed symmetric 4x4 (i<=j)
    auto set = [&](int a, int b, CCTK_REAL v){ this->at(a,b)=v; };
    set(0,0, -alp*alp + beta2);
    set(0,1, beta_d[0]); set(0,2, beta_d[1]); set(0,3, beta_d[2]);
    set(1,1, gxx); set(1,2, gxy); set(1,3, gxz);
    set(2,2, gyy); set(2,3, gyz);
    set(3,3, gzz);
  }
};

// Spacetime inverse metric (upper) from ADM identities
template<>
class inv_metric<4> : public symmetric2<CCTK_REAL,4,2> {
public:
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void from_adm(CCTK_REAL alp, CCTK_REAL betax, CCTK_REAL betay, CCTK_REAL betaz,
                CCTK_REAL gxx, CCTK_REAL gxy, CCTK_REAL gxz,
                CCTK_REAL gyy, CCTK_REAL gyz, CCTK_REAL gzz) {
    // gamma^{ij}
    inv_metric<3> gamma_uu; gamma_uu.from_metric(gxx,gxy,gxz,gyy,gyz,gzz);
    CCTK_REAL const ialpha2 = CCTK_REAL(1)/(alp*alp);

    // g^{00} = -1/alpha^2
    this->at(0,0) = -ialpha2;
    // g^{0i} = beta^i/alpha^2
    this->at(0,1) = betax * ialpha2;
    this->at(0,2) = betay * ialpha2;
    this->at(0,3) = betaz * ialpha2;
    // g^{ij} = gamma^{ij} - beta^i beta^j / alpha^2
    this->at(1,1) = gamma_uu(0,0) - betax*betax*ialpha2;
    this->at(1,2) = gamma_uu(0,1) - betax*betay*ialpha2;
    this->at(1,3) = gamma_uu(0,2) - betax*betaz*ialpha2;
    this->at(2,2) = gamma_uu(1,1) - betay*betay*ialpha2;
    this->at(2,3) = gamma_uu(1,2) - betay*betaz*ialpha2;
    this->at(3,3) = gamma_uu(2,2) - betaz*betaz*ialpha2;
  }
};

//------------------------------------------------------------------------------
// Fluid 4-velocity builder (Valencia -> u^a) with given W and v^i
// u^0 = W / alpha,  u^i = W*( v^i - beta^i / alpha )
//------------------------------------------------------------------------------
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
inline void u_from_W_v(CCTK_REAL alp, CCTK_REAL betax, CCTK_REAL betay, CCTK_REAL betaz,
                       CCTK_REAL W, CCTK_REAL vx, CCTK_REAL vy, CCTK_REAL vz,
                       CCTK_REAL* u0, CCTK_REAL* ux, CCTK_REAL* uy, CCTK_REAL* uz) {
  CCTK_REAL const ialp = CCTK_REAL(1)/alp;
  *u0 = W * ialp;
  *ux = W * (vx - betax*ialp);
  *uy = W * (vy - betay*ialp);
  *uz = W * (vz - betaz*ialp);
}

//------------------------------------------------------------------------------
// fluid_velocity_field_const: exposes u^a and v^i over the grid
//------------------------------------------------------------------------------
class fluid_velocity_field_const {
public:
  // Order: alp, betax, betay, betaz, W, velx, vely, velz
  CCTK_DEVICE CCTK_HOST
  fluid_velocity_field_const(CCTK_REAL const* alp,
                             CCTK_REAL const* betax,
                             CCTK_REAL const* betay,
                             CCTK_REAL const* betaz,
                             CCTK_REAL const* W,
                             CCTK_REAL const* velx,
                             CCTK_REAL const* vely,
                             CCTK_REAL const* velz) {
    m_data[0]=alp; m_data[1]=betax; m_data[2]=betay; m_data[3]=betaz;
    m_data[4]=W;   m_data[5]=velx;  m_data[6]=vely;  m_data[7]=velz;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get(ptrdiff_t ijk, generic<CCTK_REAL,4,1>* u) const {
    CCTK_REAL u0,ux,uy,uz;
    u_from_W_v(m_data[0][ijk], m_data[1][ijk], m_data[2][ijk], m_data[3][ijk],
               m_data[4][ijk], m_data[5][ijk], m_data[6][ijk], m_data[7][ijk],
               &u0,&ux,&uy,&uz);
    (*u)[0]=u0; (*u)[1]=ux; (*u)[2]=uy; (*u)[3]=uz;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get(ptrdiff_t ijk, generic<CCTK_REAL,3,1>* v) const {
    (*v)[0]=m_data[5][ijk]; (*v)[1]=m_data[6][ijk]; (*v)[2]=m_data[7][ijk];
  }

private:
  CCTK_REAL const* m_data[8];
};

//------------------------------------------------------------------------------
// slicing_geometry_const: geometry accessors and derived quantities
//------------------------------------------------------------------------------
class slicing_geometry_const {
public:
  // m_data layout:
  // [0]=alp, [1]=betax,[2]=betay,[3]=betaz,
  // [4..9]=gxx,gxy,gxz,gyy,gyz,gzz,
  // [10..15]=kxx,kxy,kxz,kyy,kyz,kzz, [16]=volform (sqrt(gamma))
  CCTK_DEVICE CCTK_HOST
  slicing_geometry_const(CCTK_REAL const* alp,
                         CCTK_REAL const* betax,
                         CCTK_REAL const* betay,
                         CCTK_REAL const* betaz,
                         CCTK_REAL const* gxx,
                         CCTK_REAL const* gxy,
                         CCTK_REAL const* gxz,
                         CCTK_REAL const* gyy,
                         CCTK_REAL const* gyz,
                         CCTK_REAL const* gzz,
                         CCTK_REAL const* kxx,
                         CCTK_REAL const* kxy,
                         CCTK_REAL const* kxz,
                         CCTK_REAL const* kyy,
                         CCTK_REAL const* kyz,
                         CCTK_REAL const* kzz,
                         CCTK_REAL const* volform) {
    m_data[0]=alp; m_data[1]=betax; m_data[2]=betay; m_data[3]=betaz;
    m_data[4]=gxx; m_data[5]=gxy; m_data[6]=gxz; m_data[7]=gyy; m_data[8]=gyz; m_data[9]=gzz;
    m_data[10]=kxx; m_data[11]=kxy; m_data[12]=kxz; m_data[13]=kyy; m_data[14]=kyz; m_data[15]=kzz;
    m_data[16]=volform;
  }

  // n^a = (1/alpha, -beta^i/alpha)
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get_normal(ptrdiff_t ijk, generic<CCTK_REAL,4,1>* n_u) const {
    CCTK_REAL const alp = m_data[0][ijk];
    CCTK_REAL const ialp = CCTK_REAL(1)/alp;
    (*n_u)[0] = ialp;
    (*n_u)[1] = -m_data[1][ijk]*ialp;
    (*n_u)[2] = -m_data[2][ijk]*ialp;
    (*n_u)[3] = -m_data[3][ijk]*ialp;
  }

  // n_a = (-alpha, 0, 0, 0)
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get_normal_form(ptrdiff_t ijk, generic<CCTK_REAL,4,1>* n_d) const {
    (*n_d)[0] = -m_data[0][ijk];
    (*n_d)[1] = CCTK_REAL(0);
    (*n_d)[2] = CCTK_REAL(0);
    (*n_d)[3] = CCTK_REAL(0);
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get_shift_vec(ptrdiff_t ijk, generic<CCTK_REAL,3,1>* beta_u) const {
    (*beta_u)[0] = m_data[1][ijk];
    (*beta_u)[1] = m_data[2][ijk];
    (*beta_u)[2] = m_data[3][ijk];
  }

  // 4-vector version (beta^0=0)
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get_shift_vec(ptrdiff_t ijk, generic<CCTK_REAL,4,1>* beta_u) const {
    (*beta_u)[0] = CCTK_REAL(0);
    (*beta_u)[1] = m_data[1][ijk];
    (*beta_u)[2] = m_data[2][ijk];
    (*beta_u)[3] = m_data[3][ijk];
  }

  // Direct component access
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  CCTK_REAL const* get_shift_comp(int a) const { return m_data[1 + a]; }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  CCTK_REAL const* get_space_metric_comp(int a, int b) const {
    return m_data[4 + symmetric2<CCTK_REAL,3,2>::multiindex(a,b)];
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  CCTK_REAL const* get_extr_curv_comp(int a, int b) const {
    return m_data[10 + symmetric2<CCTK_REAL,3,2>::multiindex(a,b)];
  }

  // gamma_ij at point
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get_metric(ptrdiff_t ijk, metric<3>* g) const {
    (*g)[0]=m_data[4][ijk]; (*g)[1]=m_data[5][ijk]; (*g)[2]=m_data[6][ijk];
    (*g)[3]=m_data[7][ijk]; (*g)[4]=m_data[8][ijk]; (*g)[5]=m_data[9][ijk];
  }

  // gamma^{ij} at point; if you know det (e.g., volform^2), use the overload
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get_inv_metric(ptrdiff_t ijk, inv_metric<3>* u) const {
    u->from_metric(m_data[4][ijk], m_data[5][ijk], m_data[6][ijk],
                   m_data[7][ijk], m_data[8][ijk], m_data[9][ijk]);
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get_inv_metric(ptrdiff_t ijk, inv_metric<3>* u, bool use_volform_det) const {
    if (use_volform_det) {
      CCTK_REAL const det = m_data[16][ijk]*m_data[16][ijk]; // volform = sqrt(det)
      u->from_metric(m_data[4][ijk], m_data[5][ijk], m_data[6][ijk],
                     m_data[7][ijk], m_data[8][ijk], m_data[9][ijk], det);
    } else {
      get_inv_metric(ijk,u);
    }
  }

  // Projector gamma^a_b = delta^a_b + n^a n_b
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get_space_proj(ptrdiff_t ijk, generic<CCTK_REAL,4,2>* gamma_ud) const {
    generic<CCTK_REAL,4,1> n_u, n_d;
    get_normal(ijk,&n_u);
    get_normal_form(ijk,&n_d);
    for (int a=0;a<4;++a)
    for (int b=0;b<4;++b)
      gamma_ud->at(a,b) = delta(a,b) + n_u(a)*n_d(b);
  }

  // g_ab and g^{ab}
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get_metric(ptrdiff_t ijk, metric<4>* g4_dn) const {
    g4_dn->from_adm(m_data[0][ijk], m_data[1][ijk], m_data[2][ijk], m_data[3][ijk],
                    m_data[4][ijk], m_data[5][ijk], m_data[6][ijk],
                    m_data[7][ijk], m_data[8][ijk], m_data[9][ijk]);
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get_inv_metric(ptrdiff_t ijk, inv_metric<4>* g4_uu) const {
    g4_uu->from_adm(m_data[0][ijk], m_data[1][ijk], m_data[2][ijk], m_data[3][ijk],
                    m_data[4][ijk], m_data[5][ijk], m_data[6][ijk],
                    m_data[7][ijk], m_data[8][ijk], m_data[9][ijk]);
  }

  // K_ij at point
  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE
  void get_extr_curv(ptrdiff_t ijk, symmetric2<CCTK_REAL,3,2>* k) const {
    (*k)[0]=m_data[10][ijk]; (*k)[1]=m_data[11][ijk]; (*k)[2]=m_data[12][ijk];
    (*k)[3]=m_data[13][ijk]; (*k)[4]=m_data[14][ijk]; (*k)[5]=m_data[15][ijk];
  }

private:
  CCTK_REAL const* m_data[17];
};

} // namespace tensor
} // namespace nuX_Utils

#endif // NUX_TENSOR_HXX

