#ifndef NUX_TENSOR_HXX
#define NUX_TENSOR_HXX

#include <algorithm>
#include <cstddef> // for ptrdiff_t

#include <cctk.h>

#include "nuX_metric.hxx"
#include "nuX_pow.hxx"
#include "nuX_valencia.hxx"

namespace nuX_Utils {
//! Classes for dealing with tensors and tensor fields
using namespace metric;

namespace tensor {

//! Kronecher delta symbol
CCTK_HOST CCTK_DEVICE inline CCTK_REAL delta(int a, int b) {
  return a == b ? 1.0 : 0.0;
}

///////////////////////////////////////////////////////////////////////////////
//! Generic tensor
///////////////////////////////////////////////////////////////////////////////
/*!
 *  \tparam T data type
 *  \tparam ndim number of dimensions
 *  \tparam rank tensor rank
 */
template <typename T, int ndim, int rank> class generic;

template <typename T, int ndim_> class generic<T, ndim_, 1> {
public:
  //! Data type
  typedef T data_t;
  //! Rank of the tensor
  static int const rank = 1;
  //! Number of dimensions
  static int const ndim = ndim_;
  //! Number of degrees of freedom
  static int const ndof = ndim;

  //! Computes the compressed index associated with the given indices
  /*!
   *  The data is stored in a row major format
   */
  static CCTK_HOST CCTK_DEVICE inline int multiindex(int const a) { return a; }

  //! Access the components of a tensor using a compressed index
  /*!
   *  The data is stored in a row major format
   */
  CCTK_HOST CCTK_DEVICE inline T &operator[](int const a) {
#ifdef NUX_UTILS_DEBUG
    assert(a >= 0 && a < ndof);
#endif
    return m_data[a];
  }
  //! Access the components of a tensor using a compressed index
  /*!
   *  The data is stored in a row major format
   */
  CCTK_HOST CCTK_DEVICE inline T const &operator[](int const a) const {
#ifdef NUX_UTILS_DEBUG
    assert(a >= 0 && a < ndof);
#endif
    return m_data[a];
  }

  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T &operator()(int const a) {
    return this->operator[](a);
  }
  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T const &operator()(int const a) const {
    return this->operator[](a);
  }

  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T &at(int const a) {
    return this->operator[](a);
  }
  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T const &at(int const a) const {
    return this->operator[](a);
  }
  CCTK_HOST CCTK_DEVICE inline T *data() { return &m_data[0]; }
  CCTK_HOST CCTK_DEVICE inline T const *data() const { return &m_data[0]; }

private:
  T m_data[ndim];
};

template <typename T, int ndim_> class generic<T, ndim_, 2> {
public:
  //! Data type
  typedef T data_t;
  //! Rank of the tensor
  static int const rank = 2;
  //! Number of dimensions
  static int const ndim = ndim_;
  //! Number of degrees of freedom
  static int const ndof = ndim * ndim;

  //! Computes the compressed index associated with the given indices
  /*!
   *  The data is stored in a row major format
   */
  static CCTK_HOST CCTK_DEVICE inline int multiindex(int const a, int const b) {
    return a * ndim + b;
  }

  //! Access the components of a tensor using a compressed index
  /*!
   *  The data is stored in a row major format
   */
  CCTK_HOST CCTK_DEVICE inline T &operator[](int const a) {
#ifdef NUX_UTILS_DEBUG
    assert(a >= 0 && a < ndof);
#endif
    return m_data[a];
  }
  //! Access the components of a tensor using a compressed index
  /*!
   *  The data is stored in a row major format
   */
  CCTK_HOST CCTK_DEVICE inline T const &operator[](int const a) const {
#ifdef NUX_UTILS_DEBUG
    assert(a >= 0 && a < ndof);
#endif
    return m_data[a];
  }

  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T &operator()(int const a, int const b) {
    return this->operator[](multiindex(a, b));
  }
  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T const &operator()(int const a,
                                                   int const b) const {
    return this->operator[](multiindex(a, b));
  }

  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T &at(int const a, int const b) {
    return this->operator[](multiindex(a, b));
  }
  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T const &at(int const a, int const b) const {
    return this->operator[](multiindex(a, b));
  }

private:
  T m_data[ndof];
};

template <typename T, int ndim_> class generic<T, ndim_, 3> {
public:
  //! Data type
  typedef T data_t;
  //! Rank of the tensor
  static int const rank = 3;
  //! Number of dimensions
  static int const ndim = ndim_;
  //! Number of degrees of freedom
  static int const ndof = ndim * ndim * ndim;

  //! Computes the compressed index associated with the given indices
  /*!
   *  The data is stored in a row major format
   */
  static CCTK_HOST CCTK_DEVICE inline int multiindex(int const a, int const b,
                                                     int const c) {
    return a * ndim * ndim + b * ndim + c;
  }

  //! Access the components of a tensor using a compressed index
  /*!
   *  The data is stored in a row major format
   */
  CCTK_HOST CCTK_DEVICE inline T &operator[](int const a) {
#ifdef NUX_UTILS_DEBUG
    assert(a >= 0 && a < ndof);
#endif
    return m_data[a];
  }
  //! Access the components of a tensor using a compressed index
  /*!
   *  The data is stored in a row major format
   */
  CCTK_HOST CCTK_DEVICE inline T const &operator[](int const a) const {
#ifdef NUX_UTILS_DEBUG
    assert(a >= 0 && a < ndof);
#endif
    return m_data[a];
  }

  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T &operator()(int const a, int const b,
                                             int const c) {
    return this->operator[](multiindex(a, b, c));
  }
  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T const &operator()(int const a, int const b,
                                                   int const c) const {
    return this->operator[](multiindex(a, b, c));
  }

  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T &at(int const a, int const b, int const c) {
    return this->operator[](multiindex(a, b, c));
  }
  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T const &at(int const a, int const b,
                                           int const c) const {
    return this->operator[](multiindex(a, b, c));
  }

private:
  T m_data[ndof];
};

///////////////////////////////////////////////////////////////////////////////
//! Symmetric tensor with respect to the last 2 indices
///////////////////////////////////////////////////////////////////////////////
/*!
 *  \tparam T data type
 *  \tparam ndim number of dimensions
 *  \tparam rank tensor rank
 */
template <typename T, int ndim, int rank> class symmetric2;

template <typename T, int ndim_> class symmetric2<T, ndim_, 2> {
public:
  //! Data type
  typedef T data_t;
  //! Rank of the tensor
  static int const rank = 2;
  //! Number of dimensions
  static int const ndim = ndim_;
  //! Number of degrees of freedom
  static int const ndof = (ndim * (ndim + 1)) / 2;

  //! Computes the compressed index associated with the given indices
  /*!
   *  The data is stored in a row major format
   *  Symmetric components g(i,j) with i > j are not stored so that, for
   *  instance a symmetric2<3, 2> tensor is stored as [Txx, Txy, Txz,
   *  Tyy, Tyz, Tzz].
   */
  static CCTK_HOST CCTK_DEVICE inline int multiindex(int const a, int const b) {
    if (b < a) {
      return multiindex(b, a);
    }
    int const offset = ndim * a - (a * (a - 1)) / 2;
    return offset + b - a;
  }

  //! Access the components of a tensor using a compressed index
  /*!
   *  The data is stored in a row major format.
   *  Symmetric components g(i,j) with i > j are not stored so that, for
   *  instance a symmetric2<3, 2> tensor is stored as [Txx, Txy, Txz,
   *  Tyy, Tyz, Tzz].
   */
  CCTK_HOST CCTK_DEVICE inline T &operator[](int const a) {
#ifdef NUX_UTILS_DEBUG
    assert(a >= 0 && a < ndof);
#endif
    return m_data[a];
  }
  //! Access the components of a tensor using a compressed index
  /*!
   *  The data is stored in a row major format.
   *  Symmetric components g(i,j) with i > j are not stored so that, for
   *  instance a symmetric2<3, 2> tensor is stored as [Txx, Txy, Txz,
   *  Tyy, Tyz, Tzz].
   */
  CCTK_HOST CCTK_DEVICE inline T const &operator[](int const a) const {
#ifdef NUX_UTILS_DEBUG
    assert(a >= 0 && a < ndof);
#endif
    return m_data[a];
  }

  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T &operator()(int const a, int const b) {
    return this->operator[](multiindex(a, b));
  }
  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T const &operator()(int const a,
                                                   int const b) const {
    return this->operator[](multiindex(a, b));
  }

  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T &at(int const a, int const b) {
    return this->operator[](multiindex(a, b));
  }
  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T const &at(int const a, int const b) const {
    return this->operator[](multiindex(a, b));
  }

private:
  T m_data[(ndim * (ndim + 1)) / 2];
};

template <typename T, int ndim_> class symmetric2<T, ndim_, 3> {
public:
  //! Data type
  typedef T data_t;
  //! Rank of the tensor
  static int const rank = 3;
  //! Number of dimensions
  static int const ndim = ndim_;
  //! Number of degrees of freedom
  static int const ndof = ndim * (ndim * (ndim + 1)) / 2;

  //! Computes the compressed index associated with the given indices
  /*!
   *  The data is stored in a row major format.
   *  Symmetric components g(i,j) with i > j are not stored so that, for
   *  instance a symmetric2<2, 3> tensor is stored as  [Txxx, Txxy, Txyy,
   *  Tyxx, Tyxy, Tyyy]
   */
  static CCTK_HOST CCTK_DEVICE inline int multiindex(int const a, int const b,
                                                     int const c) {
    if (c < b) {
      return multiindex(a, c, b);
    }
    int const offset = ndim * b - (b * (b - 1)) / 2;
    return a * (ndim * (ndim + 1)) / 2 + offset + c - b;
  }

  //! Access the components of a tensor using a compressed index
  /*!
   *  The data is stored in a row major format.
   *  Symmetric components g(i,j) with i > j are not stored so that, for
   *  instance a symmetric2<2, 3> tensor is stored as  [Txxx, Txxy, Txyy,
   *  Tyxx, Tyxy, Tyyy]
   */
  CCTK_HOST CCTK_DEVICE inline T &operator[](int const a) {
#ifdef NUX_UTILS_DEBUG
    assert(a >= 0 && a < ndof);
#endif
    return m_data[a];
  }
  //! Access the components of a tensor using a compressed index
  /*!
   *  The data is stored in a row major format.
   *  Symmetric components g(i,j) with i > j are not stored so that, for
   *  instance a symmetric2<2, 3> tensor is stored as  [Txxx, Txxy, Txyy,
   *  Tyxx, Tyxy, Tyyy]
   */
  CCTK_HOST CCTK_DEVICE inline T const &operator[](int const a) const {
#ifdef NUX_UTILS_DEBUG
    assert(a >= 0 && a < ndof);
#endif
    return m_data[a];
  }

  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T &operator()(int const a, int const b,
                                             int const c) {
    return this->operator[](multiindex(a, b, c));
  }
  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T const &operator()(int const a, int const b,
                                                   int const c) const {
    return this->operator[](multiindex(a, b, c));
  }

  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T &at(int const a, int const b, int const c) {
    return this->operator[](multiindex(a, b, c));
  }
  //! Access the components of a tensor using the natural indices
  CCTK_HOST CCTK_DEVICE inline T const &at(int const a, int const b,
                                           int const c) const {
    return this->operator[](multiindex(a, b, c));
  }

private:
  T m_data[ndof];
};

///////////////////////////////////////////////////////////////////////////////
//! Commonly used tensor operations (now device-callable)
///////////////////////////////////////////////////////////////////////////////
template <typename T, int N>
CCTK_HOST CCTK_DEVICE inline void contract(symmetric2<T, N, 2> const &mat,
                                           generic<T, N, 1> const &va,
                                           generic<T, N, 1> *vb) {
  for (int a = 0; a < N; ++a) {
    vb->operator[](a) = 0;
    for (int b = 0; b < N; ++b) {
      vb->operator[](a) += mat(a, b) * va(b);
    }
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE inline void contract(symmetric2<T, N, 2> const &mat,
                                           symmetric2<T, N, 2> const &A,
                                           generic<T, N, 2> *B) {
  for (int a = 0; a < N; ++a)
    for (int b = 0; b < N; ++b) {
      B->at(a, b) = 0.0;
      for (int c = 0; c < N; ++c) {
        B->at(a, b) += mat(a, c) * A(c, b);
      }
    }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE inline void contract2(symmetric2<T, N, 2> const &mat,
                                            symmetric2<T, N, 2> const &A,
                                            symmetric2<T, N, 2> *B) {
  for (int a = 0; a < N; ++a)
    for (int b = a; b < N; ++b) {
      B->at(a, b) = 0.0;
      for (int c = 0; c < N; ++c)
        for (int d = 0; d < N; ++d) {
          B->at(a, b) += mat(a, c) * mat(b, d) * A(c, d);
        }
    }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE inline T dot(generic<T, N, 1> const &va,
                                   generic<T, N, 1> const &vb) {
  T out(0);
  for (int a = 0; a < N; ++a) {
    out += va(a) * vb(a);
  }
  return out;
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE inline T dot(symmetric2<T, N, 2> const &A,
                                   symmetric2<T, N, 2> const &B) {
  T out(0);
  for (int a = 0; a < N; ++a)
    for (int b = 0; b < N; ++b) {
      out += A(a, b) * B(a, b);
    }
  return out;
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE inline T dot(symmetric2<T, N, 2> const &met,
                                   generic<T, N, 1> const &va,
                                   generic<T, N, 1> const &vb) {
  T out(0);
  for (int a = 0; a < N; ++a)
    for (int b = 0; b < N; ++b) {
      out += met(a, b) * va(a) * vb(b);
    }
  return out;
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE inline T dot(symmetric2<T, N, 2> const &met,
                                   symmetric2<T, N, 2> const &A,
                                   symmetric2<T, N, 2> const &B) {
  T out(0);
  for (int a = 0; a < N; ++a)
    for (int b = 0; b < N; ++b)
      for (int c = 0; c < N; ++c)
        for (int d = 0; d < N; ++d) {
          out += met(a, c) * met(b, d) * A(c, d) * B(a, b);
        }
  return out;
}

///////////////////////////////////////////////////////////////////////////////
//! Metric tensor
///////////////////////////////////////////////////////////////////////////////
template <int ndim> class metric;

template <int ndim> class inv_metric;

//! Spatial metric
template <> class metric<3> : public symmetric2<CCTK_REAL, 3, 2> {
public:
  //! Computes the metric determinant at the given point
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL det() const {
    return nuX_Utils::metric::spatial_det((*this)(0, 0), (*this)(0, 1),
                                          (*this)(0, 2), (*this)(1, 1),
                                          (*this)(1, 2), (*this)(2, 2));
  }
};

//! Inverse of the spatial metric
template <> class inv_metric<3> : public symmetric2<CCTK_REAL, 3, 2> {
public:
  //! Constructs the inverse metric from the metric
  CCTK_HOST CCTK_DEVICE inline void
  from_metric(CCTK_REAL const gxx, CCTK_REAL const gxy, CCTK_REAL const gxz,
              CCTK_REAL const gyy, CCTK_REAL const gyz, CCTK_REAL const gzz) {
    CCTK_REAL const det =
        nuX_Utils::metric::spatial_det(gxx, gxy, gxz, gyy, gyz, gzz);
    CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
    nuX_Utils::metric::spatial_inv(det, gxx, gxy, gxz, gyy, gyz, gzz, &uxx,
                                   &uxy, &uxz, &uyy, &uyz, &uzz);
    (*this)(0, 0) = uxx;
    (*this)(0, 1) = uxy;
    (*this)(0, 2) = uxz;
    (*this)(1, 1) = uyy;
    (*this)(1, 2) = uyz;
    (*this)(2, 2) = uzz;
  }
  //! Constructs the inverse metric from the metric and the spatial det
  CCTK_HOST CCTK_DEVICE inline void
  from_metric(CCTK_REAL const gxx, CCTK_REAL const gxy, CCTK_REAL const gxz,
              CCTK_REAL const gyy, CCTK_REAL const gyz, CCTK_REAL const gzz,
              CCTK_REAL const det) {
    CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
    nuX_Utils::metric::spatial_inv(det, gxx, gxy, gxz, gyy, gyz, gzz, &uxx,
                                   &uxy, &uxz, &uyy, &uyz, &uzz);
    (*this)(0, 0) = uxx;
    (*this)(0, 1) = uxy;
    (*this)(0, 2) = uxz;
    (*this)(1, 1) = uyy;
    (*this)(1, 2) = uyz;
    (*this)(2, 2) = uzz;
  }
  //! Constructs the inverse metric from the metric
  CCTK_HOST CCTK_DEVICE inline void from_metric(metric<3> const &g) {
    CCTK_REAL const det = g.det();
    CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
    nuX_Utils::metric::spatial_inv(det, g[0], g[1], g[2], g[3], g[4], g[5],
                                   &uxx, &uxy, &uxz, &uyy, &uyz, &uzz);
    (*this)(0, 0) = uxx;
    (*this)(0, 1) = uxy;
    (*this)(0, 2) = uxz;
    (*this)(1, 1) = uyy;
    (*this)(1, 2) = uyz;
    (*this)(2, 2) = uzz;
  }
  //! Constructs the inverse metric from the metric and the spatial det
  CCTK_HOST CCTK_DEVICE inline void from_metric_det(metric<3> const &g,
                                                    CCTK_REAL const det) {
    CCTK_REAL uxx, uxy, uxz, uyy, uyz, uzz;
    nuX_Utils::metric::spatial_inv(det, g[0], g[1], g[2], g[3], g[4], g[5],
                                   &uxx, &uxy, &uxz, &uyy, &uyz, &uzz);
    (*this)(0, 0) = uxx;
    (*this)(0, 1) = uxy;
    (*this)(0, 2) = uxz;
    (*this)(1, 1) = uyy;
    (*this)(1, 2) = uyz;
    (*this)(2, 2) = uzz;
  }
};

//! Spacetime metric
template <> class metric<4> : public symmetric2<CCTK_REAL, 4, 2> {
public:
  //! Construct the spacetime metric from the ADM quantities
  CCTK_HOST CCTK_DEVICE inline void
  from_adm(CCTK_REAL const alp, CCTK_REAL const betax, CCTK_REAL const betay,
           CCTK_REAL const betaz, CCTK_REAL const gxx, CCTK_REAL const gxy,
           CCTK_REAL const gxz, CCTK_REAL const gyy, CCTK_REAL const gyz,
           CCTK_REAL const gzz) {
    CCTK_REAL g[16];
    nuX_Utils::metric::spacetime(alp, betax, betay, betaz, gxx, gxy, gxz, gyy,
                                 gyz, gzz, &g[0]);
    for (int a = 0; a < 4; ++a)
      for (int b = a; b < 4; ++b) {
        this->at(a, b) = g[4 * a + b];
      }
  }
};

//! Spacetime inverse metric
template <> class inv_metric<4> : public symmetric2<CCTK_REAL, 4, 2> {
public:
  //! Construct the spacetime metric from the ADM quantities
  CCTK_HOST CCTK_DEVICE inline void
  from_adm(CCTK_REAL const alp, CCTK_REAL const betax, CCTK_REAL const betay,
           CCTK_REAL const betaz, CCTK_REAL const gxx, CCTK_REAL const gxy,
           CCTK_REAL const gxz, CCTK_REAL const gyy, CCTK_REAL const gyz,
           CCTK_REAL const gzz) {
    CCTK_REAL u[16];
    nuX_Utils::metric::spacetime_upper(alp, betax, betay, betaz, gxx, gxy, gxz,
                                       gyy, gyz, gzz, &u[0]);
    for (int a = 0; a < 4; ++a)
      for (int b = a; b < 4; ++b) {
        this->at(a, b) = u[4 * a + b];
      }
  }
};

///////////////////////////////////////////////////////////////////////////////
// Special tensor fields
///////////////////////////////////////////////////////////////////////////////
//! Fluid four velocity (as a vector)
class fluid_velocity_field_const {
public:
  //! Initialize the fluid 4 velocity
  fluid_velocity_field_const(CCTK_REAL const *alp, CCTK_REAL const *betax,
                             CCTK_REAL const *betay, CCTK_REAL const *betaz,
                             CCTK_REAL const *w_lorentz, CCTK_REAL const *velx,
                             CCTK_REAL const *vely, CCTK_REAL const *velz) {
    m_data[0] = alp;
    m_data[1] = betax;
    m_data[2] = betay;
    m_data[3] = betaz;
    m_data[4] = w_lorentz;
    m_data[5] = velx;
    m_data[6] = vely;
    m_data[7] = velz;
#ifdef NUX_UTILS_DEBUG
    for (int i = 0; i < 8; ++i) {
      assert(m_data[i]);
    }
#endif
  }

  //! Evaluate the fluid four velocity at a given location
  CCTK_HOST CCTK_DEVICE inline void get(ptrdiff_t const ijk,
                                        generic<CCTK_REAL, 4, 1> *u) const {
    nuX_Utils::valencia::uvel(m_data[0][ijk], m_data[1][ijk], m_data[2][ijk],
                              m_data[3][ijk], m_data[4][ijk], m_data[5][ijk],
                              m_data[6][ijk], m_data[7][ijk], u->data());
  }
  //! Evaluate the three velocity at a given location
  CCTK_HOST CCTK_DEVICE inline void get(ptrdiff_t const ijk,
                                        generic<CCTK_REAL, 3, 1> *v) const {
    (*v)(0) = m_data[5][ijk];
    (*v)(1) = m_data[6][ijk];
    (*v)(2) = m_data[7][ijk];
  }

private:
  CCTK_REAL const *m_data[8];
};

//! Class describing the geometry of the ADM slicing
class slicing_geometry_const {
public:
  //! Initialize the slicing geometry from the ADM quantities
  slicing_geometry_const(CCTK_REAL const *alp, CCTK_REAL const *betax,
                         CCTK_REAL const *betay, CCTK_REAL const *betaz,
                         CCTK_REAL const *gxx, CCTK_REAL const *gxy,
                         CCTK_REAL const *gxz, CCTK_REAL const *gyy,
                         CCTK_REAL const *gyz, CCTK_REAL const *gzz,
                         CCTK_REAL const *kxx, CCTK_REAL const *kxy,
                         CCTK_REAL const *kxz, CCTK_REAL const *kyy,
                         CCTK_REAL const *kyz, CCTK_REAL const *kzz,
                         CCTK_REAL const *volform) {
    m_data[0] = alp;
    m_data[1] = betax;
    m_data[2] = betay;
    m_data[3] = betaz;
    m_data[4] = gxx;
    m_data[5] = gxy;
    m_data[6] = gxz;
    m_data[7] = gyy;
    m_data[8] = gyz;
    m_data[9] = gzz;
    m_data[10] = kxx;
    m_data[11] = kxy;
    m_data[12] = kxz;
    m_data[13] = kyy;
    m_data[14] = kyz;
    m_data[15] = kzz;
    m_data[16] = volform;
#ifdef NUX_UTILS_DEBUG
    for (int i = 0; i < 16; ++i) {
      assert(m_data[i]);
    }
#endif
  }

  //! Get the normal vector to the spacelike hypersurface
  CCTK_HOST CCTK_DEVICE inline void
  get_normal(ptrdiff_t const ijk, generic<CCTK_REAL, 4, 1> *n) const {
    nuX_Utils::metric::normal(m_data[0][ijk], m_data[1][ijk], m_data[2][ijk],
                              m_data[3][ijk], n->data());
  }
  //! Get the normal one-form to the spacelike hypersurface
  CCTK_HOST CCTK_DEVICE inline void
  get_normal_form(ptrdiff_t const ijk, generic<CCTK_REAL, 4, 1> *n_d) const {
    (*n_d)(0) = -m_data[0][ijk];
    (*n_d)(1) = 0.0;
    (*n_d)(2) = 0.0;
    (*n_d)(3) = 0.0;
  }

  //! Get the shift vector
  CCTK_HOST CCTK_DEVICE inline void
  get_shift_vec(ptrdiff_t const ijk, generic<CCTK_REAL, 3, 1> *beta_u) const {
    (*beta_u)(0) = m_data[1][ijk];
    (*beta_u)(1) = m_data[2][ijk];
    (*beta_u)(2) = m_data[3][ijk];
  }
  //! Get the shift 4-vector β^μ = (0, β^i)
  CCTK_HOST CCTK_DEVICE inline void
  get_shift_vec(ptrdiff_t const ijk, generic<CCTK_REAL, 4, 1> *beta_u) const {
    (*beta_u)(0) = 0.0;
    (*beta_u)(1) = m_data[1][ijk];
    (*beta_u)(2) = m_data[2][ijk];
    (*beta_u)(3) = m_data[3][ijk];
  }

  //! Get a component of the shift vector on the entire grid
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL const *
  get_shift_comp(int const a) const {
    return m_data[1 + a];
  }
  //! Get a component of the spatial metric on the entire grid
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL const *
  get_space_metric_comp(int const a, int const b) const {
    return m_data[4 + metric<3>::multiindex(a, b)];
  }
  //! Get a component of the extrinsic curvature on the entire grid
  CCTK_HOST CCTK_DEVICE inline CCTK_REAL const *
  get_extr_curv_comp(int const a, int const b) const {
    return m_data[10 + symmetric2<CCTK_REAL, 3, 2>::multiindex(a, b)];
  }

  //! Get the spatial metric at a given location
  CCTK_HOST CCTK_DEVICE inline void get_metric(ptrdiff_t const ijk,
                                               metric<3> *g) const {
    (*g)(0, 0) = m_data[4][ijk];
    (*g)(0, 1) = m_data[5][ijk];
    (*g)(0, 2) = m_data[6][ijk];
    (*g)(1, 1) = m_data[7][ijk];
    (*g)(1, 2) = m_data[8][ijk];
    (*g)(2, 2) = m_data[9][ijk];
  }
  //! Get inverse spatial metric at a given location
  CCTK_HOST CCTK_DEVICE inline void get_inv_metric(ptrdiff_t const ijk,
                                                   inv_metric<3> *u) const {
    u->from_metric(m_data[4][ijk], m_data[5][ijk], m_data[6][ijk],
                   m_data[7][ijk], m_data[8][ijk], m_data[9][ijk],
                   m_data[16][ijk] * m_data[16][ijk]);
  }

  //! \brief Get the projector onto the spacelike hypersurface:
  //! \f$ \gamma^a_{\phantom{a}b} \f$
  CCTK_HOST CCTK_DEVICE inline void
  get_space_proj(ptrdiff_t const ijk,
                 generic<CCTK_REAL, 4, 2> *gamma_ud) const {
    generic<CCTK_REAL, 4, 1> n_d;
    generic<CCTK_REAL, 4, 1> n_u;
    get_normal(ijk, &n_u);
    get_normal_form(ijk, &n_d);

    for (int a = 0; a < 4; ++a)
      for (int b = 0; b < 4; ++b) {
        gamma_ud->at(a, b) = delta(a, b) + n_u(a) * n_d(b);
      }
  }

  //! Get the spacetime metric at a given location
  CCTK_HOST CCTK_DEVICE inline void get_metric(ptrdiff_t const ijk,
                                               metric<4> *g) const {
    g->from_adm(m_data[0][ijk], m_data[1][ijk], m_data[2][ijk], m_data[3][ijk],
                m_data[4][ijk], m_data[5][ijk], m_data[6][ijk], m_data[7][ijk],
                m_data[8][ijk], m_data[9][ijk]);
  }

  //! Get the spacetime inverse metric at a given location
  CCTK_HOST CCTK_DEVICE inline void get_inv_metric(ptrdiff_t const ijk,
                                                   inv_metric<4> *u) const {
    u->from_adm(m_data[0][ijk], m_data[1][ijk], m_data[2][ijk], m_data[3][ijk],
                m_data[4][ijk], m_data[5][ijk], m_data[6][ijk], m_data[7][ijk],
                m_data[8][ijk], m_data[9][ijk]);
  }

  //! Get the extrinsic curvature at a given location
  CCTK_HOST CCTK_DEVICE inline void
  get_extr_curv(ptrdiff_t const ijk, symmetric2<CCTK_REAL, 3, 2> *k) const {
    (*k)(0, 0) = m_data[10][ijk];
    (*k)(0, 1) = m_data[11][ijk];
    (*k)(0, 2) = m_data[12][ijk];
    (*k)(1, 1) = m_data[13][ijk];
    (*k)(1, 2) = m_data[14][ijk];
    (*k)(2, 2) = m_data[15][ijk];
  }

private:
  CCTK_REAL const *m_data[17];
};

} // namespace tensor
} // namespace nuX_Utils

#endif
