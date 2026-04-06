#ifndef NUX_ROOTS_HXX
#define NUX_ROOTS_HXX

#include <algorithm>
#include <cmath>
#include <limits>

#include <cctk.h>
#include <mat.hxx>
#include <vec.hxx>

namespace nuX_Utils {
namespace roots {

enum class status : int {
  success = 0,
  failure = -1,
  continue_iter = -2,
  edom = 1,
  erange = 2,
  efault = 3,
  einval = 4,
  efailed = 5,
  efactor = 6,
  esanity = 7,
  enomem = 8,
  ebadfunc = 9,
  erunaway = 10,
  emaxiter = 11,
  ezerodiv = 12,
  ebadtol = 13,
  etol = 14,
  eundrflw = 15,
  eovrflw = 16,
  eloss = 17,
  eround = 18,
  ebadlen = 19,
  enotsqr = 20,
  esing = 21,
  ediverge = 22,
  eunsup = 23,
  eunimpl = 24,
  ecache = 25,
  etable = 26,
  enoprog = 27,
  enoprogj = 28,
  etolf = 29,
  etolx = 30,
  etolg = 31,
  eof = 32
};

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
code(const status s) {
  return static_cast<int>(s);
}

template <typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool finite(
    const T x) {
  using std::isfinite;
  return isfinite(x);
}

template <typename F, typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
safe_func_call(const F &f, const T x, T *y) {
  *y = f(x);
  return finite(*y) ? code(status::success) : code(status::ebadfunc);
}

template <typename T> struct brent_solver {
  T root{};
  T x_lower{};
  T x_upper{};
  T a{};
  T b{};
  T c{};
  T d{};
  T e{};
  T fa{};
  T fb{};
  T fc{};
};

template <typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
solver_root(const brent_solver<T> *const s) {
  return s->root;
}

template <typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
solver_x_lower(const brent_solver<T> *const s) {
  return s->x_lower;
}

template <typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
solver_x_upper(const brent_solver<T> *const s) {
  return s->x_upper;
}

template <typename F, typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
solver_set(brent_solver<T> *const s, const F &f, const T x_lower,
           const T x_upper) {
  T f_lower{};
  T f_upper{};

  s->root = T(0.5) * (x_lower + x_upper);

  int ierr = safe_func_call(f, x_lower, &f_lower);
  if (ierr != code(status::success))
    return ierr;

  ierr = safe_func_call(f, x_upper, &f_upper);
  if (ierr != code(status::success))
    return ierr;

  s->x_lower = x_lower;
  s->x_upper = x_upper;
  s->a = x_lower;
  s->fa = f_lower;
  s->b = x_upper;
  s->fb = f_upper;
  s->c = x_upper;
  s->fc = f_upper;
  s->d = x_upper - x_lower;
  s->e = x_upper - x_lower;

  if ((f_lower < T(0) && f_upper < T(0)) ||
      (f_lower > T(0) && f_upper > T(0))) {
    return code(status::einval);
  }

  return code(status::success);
}

template <typename F, typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
solver_iterate(brent_solver<T> *const s, const F &f) {
  using std::abs;
  using std::min;

  auto swap1 = [](T &x, T &y) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    T z = x;
    x = y;
    y = z;
  };

  int ac_equal = 0;

  T a = s->a;
  T b = s->b;
  T c = s->c;
  T fa = s->fa;
  T fb = s->fb;
  T fc = s->fc;
  T d = s->d;
  T e = s->e;

  if ((fb < T(0) && fc < T(0)) || (fb > T(0) && fc > T(0))) {
    ac_equal = 1;
    c = a;
    fc = fa;
    d = b - a;
    e = b - a;
  }

  if (abs(fc) < abs(fb)) {
    ac_equal = 1;
    a = b;
    b = c;
    c = a;
    fa = fb;
    fb = fc;
    fc = fa;
  }

  const T tol = T(0.5) * std::numeric_limits<T>::epsilon() * abs(b);
  const T m = T(0.5) * (c - b);

  if (fb == T(0)) {
    s->root = b;
    s->x_lower = b;
    s->x_upper = b;
    return code(status::success);
  }

  if (abs(m) <= tol) {
    s->root = b;
    s->x_lower = min(b, c);
    s->x_upper = (b < c) ? c : b;
    return code(status::success);
  }

  if (abs(e) < tol || abs(fa) <= abs(fb)) {
    d = m;
    e = m;
  } else {
    T p{};
    T q{};
    T r{};
    const T sfa = fb / fa;

    if (ac_equal) {
      p = T(2) * m * sfa;
      q = T(1) - sfa;
    } else {
      q = fa / fc;
      r = fb / fc;
      p = sfa * (T(2) * m * q * (q - r) - (b - a) * (r - T(1)));
      q = (q - T(1)) * (r - T(1)) * (sfa - T(1));
    }

    if (p > T(0)) {
      q = -q;
    } else {
      p = -p;
    }

    if (T(2) * p <
        std::min(T(3) * m * q - abs(tol * q), abs(e * q))) {
      e = d;
      d = p / q;
    } else {
      d = m;
      e = m;
    }
  }

  a = b;
  fa = fb;

  if (abs(d) > tol) {
    b += d;
  } else {
    b += (m > T(0) ? +tol : -tol);
  }

  int ierr = safe_func_call(f, b, &fb);
  if (ierr != code(status::success))
    return ierr;

  s->a = a;
  s->b = b;
  s->c = c;
  s->d = d;
  s->e = e;
  s->fa = fa;
  s->fb = fb;
  s->fc = fc;
  s->root = b;

  if ((fb < T(0) && fc < T(0)) || (fb > T(0) && fc > T(0))) {
    c = a;
  }

  s->x_lower = min(b, c);
  s->x_upper = (b < c) ? c : b;

  return code(status::success);
}

template <typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
test_interval(const T x_lower, const T x_upper, const T epsabs,
              const T epsrel) {
  using std::abs;

  if (epsrel < T(0))
    return code(status::ebadtol);
  if (epsabs < T(0))
    return code(status::ebadtol);
  if (x_lower > x_upper)
    return code(status::einval);

  const T abs_lower = abs(x_lower);
  const T abs_upper = abs(x_upper);
  const T min_abs = ((x_lower > T(0) && x_upper > T(0)) ||
                     (x_lower < T(0) && x_upper < T(0)))
                        ? std::min(abs_lower, abs_upper)
                        : T(0);
  const T tolerance = epsabs + epsrel * min_abs;

  return (abs(x_upper - x_lower) < tolerance)
             ? code(status::success)
             : code(status::continue_iter);
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
test_delta(const Arith::vec<T, N> &dx, const Arith::vec<T, N> &x,
           const T epsabs, const T epsrel) {
  using std::abs;

  if (epsrel < T(0))
    return code(status::ebadtol);

  for (int i = 0; i < N; ++i) {
    const T tolerance = epsabs + epsrel * abs(x(i));
    if (!(abs(dx(i)) < tolerance || dx(i) == T(0))) {
      return code(status::continue_iter);
    }
  }

  return code(status::success);
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
test_residual(const Arith::vec<T, N> &f, const T epsabs) {
  using std::abs;

  if (epsabs < T(0))
    return code(status::ebadtol);

  T residual = T(0);
  for (int i = 0; i < N; ++i) {
    residual += abs(f(i));
  }

  return (residual < epsabs) ? code(status::success)
                             : code(status::continue_iter);
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
finite(const Arith::vec<T, N> &x) {
  for (int i = 0; i < N; ++i) {
    if (!finite(x(i))) {
      return false;
    }
  }
  return true;
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline bool
finite(const Arith::mat<T, N> &A) {
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      if (!finite(A(i, j))) {
        return false;
      }
    }
  }
  return true;
}

template <typename T, int N> struct hybridsj_solver {
  int iter{};
  int ncfail{};
  int ncsuc{};
  int nslow1{};
  int nslow2{};
  T fnorm{};
  T delta{};
  Arith::vec<T, N> x;
  Arith::vec<T, N> f;
  Arith::vec<T, N> dx;
  Arith::mat<T, N> J;
  Arith::mat<T, N> q;
  Arith::mat<T, N> r;
  Arith::vec<T, N> diag;
  Arith::vec<T, N> qtf;
  Arith::vec<T, N> newton;
  Arith::vec<T, N> gradient;
  Arith::vec<T, N> x_trial;
  Arith::vec<T, N> f_trial;
  Arith::vec<T, N> df;
  Arith::vec<T, N> qtdf;
  Arith::vec<T, N> rdx;
  Arith::vec<T, N> w;
  Arith::vec<T, N> v;
};

namespace detail {

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
set_zero(Arith::vec<T, N> *x) {
  for (int i = 0; i < N; ++i) {
    (*x)(i) = T(0);
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
set_zero(Arith::mat<T, N> *A) {
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      (*A)(i, j) = T(0);
    }
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
set_identity(Arith::mat<T, N> *A) {
  set_zero(A);
  for (int i = 0; i < N; ++i) {
    (*A)(i, i) = T(1);
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
enorm(const Arith::vec<T, N> &f) {
  T e2 = T(0);
  for (int i = 0; i < N; ++i) {
    e2 += f(i) * f(i);
  }
  return sqrt(e2);
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
scaled_enorm(const Arith::vec<T, N> &d, const Arith::vec<T, N> &f) {
  T e2 = T(0);
  for (int i = 0; i < N; ++i) {
    const T u = d(i) * f(i);
    e2 += u * u;
  }
  return sqrt(e2);
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
enorm_sum(const Arith::vec<T, N> &a, const Arith::vec<T, N> &b) {
  T e2 = T(0);
  for (int i = 0; i < N; ++i) {
    const T u = a(i) + b(i);
    e2 += u * u;
  }
  return sqrt(e2);
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
compute_df(const Arith::vec<T, N> &f_trial, const Arith::vec<T, N> &f,
           Arith::vec<T, N> *df) {
  for (int i = 0; i < N; ++i) {
    (*df)(i) = f_trial(i) - f(i);
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
compute_diag(const Arith::mat<T, N> &J, Arith::vec<T, N> *diag) {
  for (int j = 0; j < N; ++j) {
    T sum = T(0);
    for (int i = 0; i < N; ++i) {
      const T Jij = J(i, j);
      sum += Jij * Jij;
    }
    if (sum == T(0)) {
      sum = T(1);
    }
    (*diag)(j) = sqrt(sum);
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
update_diag(const Arith::mat<T, N> &J, Arith::vec<T, N> *diag) {
  for (int j = 0; j < N; ++j) {
    T sum = T(0);
    for (int i = 0; i < N; ++i) {
      const T Jij = J(i, j);
      sum += Jij * Jij;
    }
    if (sum == T(0)) {
      sum = T(1);
    }
    const T cnorm = sqrt(sum);
    if (cnorm > (*diag)(j)) {
      (*diag)(j) = cnorm;
    }
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
compute_delta(const Arith::vec<T, N> &diag, const Arith::vec<T, N> &x) {
  const T Dx = scaled_enorm(diag, x);
  const T factor = T(100);
  return (Dx > T(0)) ? factor * Dx : factor;
}

template <typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
compute_actual_reduction(const T fnorm, const T fnorm1) {
  if (fnorm1 < fnorm) {
    const T u = fnorm1 / fnorm;
    return T(1) - u * u;
  }
  return T(-1);
}

template <typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
compute_predicted_reduction(const T fnorm, const T fnorm1) {
  if (fnorm1 < fnorm) {
    const T u = fnorm1 / fnorm;
    return T(1) - u * u;
  }
  return T(0);
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
compute_qtf(const Arith::mat<T, N> &q, const Arith::vec<T, N> &f,
            Arith::vec<T, N> *qtf) {
  for (int j = 0; j < N; ++j) {
    T sum = T(0);
    for (int i = 0; i < N; ++i) {
      sum += q(i, j) * f(i);
    }
    (*qtf)(j) = sum;
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
compute_rdx(const Arith::mat<T, N> &r, const Arith::vec<T, N> &dx,
            Arith::vec<T, N> *rdx) {
  for (int i = 0; i < N; ++i) {
    T sum = T(0);
    for (int j = i; j < N; ++j) {
      sum += r(i, j) * dx(j);
    }
    (*rdx)(i) = sum;
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
compute_trial_step(const Arith::vec<T, N> &x, const Arith::vec<T, N> &dx,
                   Arith::vec<T, N> *x_trial) {
  for (int i = 0; i < N; ++i) {
    (*x_trial)(i) = x(i) + dx(i);
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
r_solve(const Arith::mat<T, N> &R, const Arith::vec<T, N> &b,
        Arith::vec<T, N> *x) {
  using std::abs;

  for (int i = N - 1; i >= 0; --i) {
    T sum = b(i);
    for (int j = i + 1; j < N; ++j) {
      sum -= R(i, j) * (*x)(j);
    }
    if (abs(R(i, i)) == T(0)) {
      return code(status::ezerodiv);
    }
    (*x)(i) = sum / R(i, i);
  }
  return code(status::success);
}

template <typename T>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
givens(const T a, const T b, T *c, T *s) {
  using std::abs;
  using std::sqrt;

  if (b == T(0)) {
    *c = T(1);
    *s = T(0);
  } else if (abs(b) > abs(a)) {
    const T t = -a / b;
    const T s1 = T(1) / sqrt(T(1) + t * t);
    *s = s1;
    *c = s1 * t;
  } else {
    const T t = -b / a;
    const T c1 = T(1) / sqrt(T(1) + t * t);
    *c = c1;
    *s = c1 * t;
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
givens_gv(Arith::vec<T, N> *v, const int i, const int j, const T c,
          const T s) {
  const T vi = (*v)(i);
  const T vj = (*v)(j);
  (*v)(i) = c * vi - s * vj;
  (*v)(j) = s * vi + c * vj;
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
apply_givens_qr(Arith::mat<T, N> *Q, Arith::mat<T, N> *R, const int i,
                const int j, const T c, const T s) {
  for (int k = 0; k < N; ++k) {
    const T qki = (*Q)(k, i);
    const T qkj = (*Q)(k, j);
    (*Q)(k, i) = qki * c - qkj * s;
    (*Q)(k, j) = qki * s + qkj * c;
  }

  for (int k = ((i < j) ? i : j); k < N; ++k) {
    const T rik = (*R)(i, k);
    const T rjk = (*R)(j, k);
    (*R)(i, k) = c * rik - s * rjk;
    (*R)(j, k) = s * rik + c * rjk;
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline T
householder_transform_column(Arith::mat<T, N> *A, const int col) {
  using std::abs;
  using std::hypot;

  if (col == N - 1) {
    return T(0);
  }

  T xnorm = T(0);
  for (int i = col + 1; i < N; ++i) {
    xnorm += (*A)(i, col) * (*A)(i, col);
  }
  xnorm = sqrt(xnorm);

  if (xnorm == T(0)) {
    return T(0);
  }

  const T alpha = (*A)(col, col);
  const T beta = -((alpha >= T(0)) ? T(1) : T(-1)) * hypot(alpha, xnorm);
  const T tau = (beta - alpha) / beta;

  const T s = alpha - beta;
  if (abs(s) > std::numeric_limits<T>::min()) {
    for (int i = col + 1; i < N; ++i) {
      (*A)(i, col) /= s;
    }
  } else {
    const T eps = std::numeric_limits<T>::epsilon();
    for (int i = col + 1; i < N; ++i) {
      (*A)(i, col) *= eps / s;
    }
    for (int i = col + 1; i < N; ++i) {
      (*A)(i, col) /= eps;
    }
  }

  (*A)(col, col) = beta;
  return tau;
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
householder_apply_left_from_column(const T tau, const Arith::mat<T, N> &src,
                                   const int col, Arith::mat<T, N> *A,
                                   const int row0, const int col0) {
  if (tau == T(0)) {
    return;
  }

  for (int j = col0; j < N; ++j) {
    T wj = (*A)(row0, j);
    for (int i = row0 + 1; i < N; ++i) {
      wj += (*A)(i, j) * src(i, col);
    }

    (*A)(row0, j) -= tau * wj;
    for (int i = row0 + 1; i < N; ++i) {
      (*A)(i, j) -= tau * src(i, col) * wj;
    }
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
qr_decomp_unpack(const Arith::mat<T, N> &A, Arith::mat<T, N> *Q,
                 Arith::mat<T, N> *R) {
  Arith::mat<T, N> qr = A;
  Arith::vec<T, N> tau;
  set_zero(&tau);

  for (int i = 0; i < N; ++i) {
    tau(i) = householder_transform_column(&qr, i);
    if (i + 1 < N) {
      householder_apply_left_from_column(tau(i), qr, i, &qr, i, i + 1);
    }
  }

  set_identity(Q);
  for (int i = N - 1; i >= 0; --i) {
    householder_apply_left_from_column(tau(i), qr, i, Q, i, i);
  }

  set_zero(R);
  for (int i = 0; i < N; ++i) {
    for (int j = i; j < N; ++j) {
      (*R)(i, j) = qr(i, j);
    }
  }

  return code(status::success);
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
compute_wv(const Arith::vec<T, N> &qtdf, const Arith::vec<T, N> &rdx,
           const Arith::vec<T, N> &dx, const Arith::vec<T, N> &diag,
           const T pnorm, Arith::vec<T, N> *w, Arith::vec<T, N> *v) {
  for (int i = 0; i < N; ++i) {
    (*w)(i) = (qtdf(i) - rdx(i)) / pnorm;
    (*v)(i) = diag(i) * diag(i) * dx(i) / pnorm;
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
qr_update(Arith::mat<T, N> *Q, Arith::mat<T, N> *R, Arith::vec<T, N> *w,
          const Arith::vec<T, N> &v) {
  for (int k = N - 1; k > 0; --k) {
    T c{};
    T s{};
    givens((*w)(k - 1), (*w)(k), &c, &s);
    givens_gv(w, k - 1, k, c, s);
    apply_givens_qr(Q, R, k - 1, k, c, s);
  }

  const T w0 = (*w)(0);
  for (int j = 0; j < N; ++j) {
    (*R)(0, j) += w0 * v(j);
  }

  for (int k = 1; k < N; ++k) {
    T c{};
    T s{};
    givens((*R)(k - 1, k - 1), (*R)(k, k - 1), &c, &s);
    apply_givens_qr(Q, R, k - 1, k, c, s);
    (*R)(k, k - 1) = T(0);
  }

  return code(status::success);
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
gradient_direction(const Arith::mat<T, N> &r, const Arith::vec<T, N> &qtf,
                   const Arith::vec<T, N> &diag, Arith::vec<T, N> *g) {
  for (int j = 0; j < N; ++j) {
    T sum = T(0);
    for (int i = 0; i < N; ++i) {
      sum += r(i, j) * qtf(i);
    }
    (*g)(j) = -sum / diag(j);
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
minimum_step(const T gnorm, const Arith::vec<T, N> &diag,
             Arith::vec<T, N> *g) {
  for (int i = 0; i < N; ++i) {
    (*g)(i) = ((*g)(i) / gnorm) / diag(i);
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
compute_Rg(const Arith::mat<T, N> &r, const Arith::vec<T, N> &gradient,
           Arith::vec<T, N> *Rg) {
  for (int i = 0; i < N; ++i) {
    T sum = T(0);
    for (int j = i; j < N; ++j) {
      sum += r(i, j) * gradient(j);
    }
    (*Rg)(i) = sum;
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
scaled_addition(const T alpha, const Arith::vec<T, N> &newton, const T beta,
                const Arith::vec<T, N> &gradient, Arith::vec<T, N> *p) {
  for (int i = 0; i < N; ++i) {
    (*p)(i) = alpha * newton(i) + beta * gradient(i);
  }
}

template <typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
dogleg(const Arith::mat<T, N> &r, const Arith::vec<T, N> &qtf,
       const Arith::vec<T, N> &diag, const T delta,
       Arith::vec<T, N> *newton, Arith::vec<T, N> *gradient,
       Arith::vec<T, N> *p) {
  int ierr = r_solve(r, qtf, newton);
  if (ierr != code(status::success)) {
    return ierr;
  }

  for (int i = 0; i < N; ++i) {
    (*newton)(i) = -(*newton)(i);
  }

  const T qnorm = scaled_enorm(diag, *newton);
  if (qnorm <= delta) {
    *p = *newton;
    return code(status::success);
  }

  gradient_direction(r, qtf, diag, gradient);

  const T gnorm = enorm(*gradient);
  if (gnorm == T(0)) {
    scaled_addition(delta / qnorm, *newton, T(0), *gradient, p);
    return code(status::success);
  }

  minimum_step(gnorm, diag, gradient);
  compute_Rg(r, *gradient, p);

  const T temp = enorm(*p);
  const T sgnorm = (gnorm / temp) / temp;

  if (sgnorm > delta) {
    scaled_addition(T(0), *newton, delta, *gradient, p);
    return code(status::success);
  }

  const T bnorm = enorm(qtf);
  const T bg = bnorm / gnorm;
  const T bq = bnorm / qnorm;
  const T dq = delta / qnorm;
  const T dq2 = dq * dq;
  const T sd = sgnorm / delta;
  const T sd2 = sd * sd;
  const T t1 = bg * bq * sd;
  const T u = t1 - dq;
  const T t2 = t1 - dq * sd2 + sqrt(u * u + (T(1) - dq2) * (T(1) - sd2));
  const T alpha = dq * (T(1) - sd2) / t2;
  const T beta = (T(1) - alpha) * sgnorm;

  scaled_addition(alpha, *newton, beta, *gradient, p);
  return code(status::success);
}

} // namespace detail

template <typename FDF, typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
hybridsj_set(hybridsj_solver<T, N> *const s, const FDF &fdf,
             const Arith::vec<T, N> &x0) {
  s->x = x0;

  int ierr = fdf(s->x, s->f, s->J);
  if (ierr != code(status::success) || !finite(s->f) || !finite(s->J)) {
    return code(status::ebadfunc);
  }

  s->iter = 1;
  s->fnorm = detail::enorm(s->f);
  s->ncfail = 0;
  s->ncsuc = 0;
  s->nslow1 = 0;
  s->nslow2 = 0;

  detail::set_zero(&s->dx);
  detail::compute_diag(s->J, &s->diag);
  s->delta = detail::compute_delta(s->diag, s->x);

  return detail::qr_decomp_unpack(s->J, &s->q, &s->r);
}

template <typename F, typename DF, typename T, int N>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline int
hybridsj_iterate(hybridsj_solver<T, N> *const s, const F &f, const DF &df) {
  using std::abs;

  const T fnorm = s->fnorm;
  const T p1 = T(0.1);
  const T p5 = T(0.5);
  const T p001 = T(0.001);
  const T p0001 = T(0.0001);

  detail::compute_qtf(s->q, s->f, &s->qtf);

  int ierr = detail::dogleg(s->r, s->qtf, s->diag, s->delta, &s->newton,
                            &s->gradient, &s->dx);
  if (ierr != code(status::success)) {
    return ierr;
  }

  detail::compute_trial_step(s->x, s->dx, &s->x_trial);

  const T pnorm = detail::scaled_enorm(s->diag, s->dx);
  if (s->iter == 1 && pnorm < s->delta) {
    s->delta = pnorm;
  }

  ierr = f(s->x_trial, s->f_trial);
  if (ierr != code(status::success) || !finite(s->f_trial)) {
    return code(status::ebadfunc);
  }

  detail::compute_df(s->f_trial, s->f, &s->df);

  const T fnorm1 = detail::enorm(s->f_trial);
  const T actred = detail::compute_actual_reduction(fnorm, fnorm1);

  detail::compute_rdx(s->r, s->dx, &s->rdx);

  const T fnorm1p = detail::enorm_sum(s->qtf, s->rdx);
  const T prered = detail::compute_predicted_reduction(fnorm, fnorm1p);
  const T ratio = (prered > T(0)) ? actred / prered : T(0);

  if (ratio < p1) {
    s->ncsuc = 0;
    s->ncfail++;
    s->delta *= p5;
  } else {
    s->ncfail = 0;
    s->ncsuc++;

    if (ratio >= p5 || s->ncsuc > 1) {
      s->delta = std::max(s->delta, pnorm / p5);
    }
    if (abs(ratio - T(1)) <= p1) {
      s->delta = pnorm / p5;
    }
  }

  if (ratio >= p0001) {
    s->x = s->x_trial;
    s->f = s->f_trial;
    s->fnorm = fnorm1;
    s->iter++;
  }

  s->nslow1++;
  if (actred >= p001) {
    s->nslow1 = 0;
  }
  if (actred >= p1) {
    s->nslow2 = 0;
  }

  if (s->ncfail == 2) {
    ierr = df(s->x, s->J);
    if (ierr != code(status::success) || !finite(s->J)) {
      return code(status::ebadfunc);
    }

    s->nslow2++;
    if (s->iter == 1) {
      detail::compute_diag(s->J, &s->diag);
      s->delta = detail::compute_delta(s->diag, s->x);
    } else {
      detail::update_diag(s->J, &s->diag);
    }

    return detail::qr_decomp_unpack(s->J, &s->q, &s->r);
  }

  detail::compute_qtf(s->q, s->df, &s->qtdf);
  detail::compute_wv(s->qtdf, s->rdx, s->dx, s->diag, pnorm, &s->w, &s->v);

  ierr = detail::qr_update(&s->q, &s->r, &s->w, s->v);
  if (ierr != code(status::success)) {
    return ierr;
  }

  if (s->nslow2 == 5) {
    return code(status::enoprogj);
  }
  if (s->nslow1 == 10) {
    return code(status::enoprog);
  }

  return code(status::success);
}

} // namespace roots
} // namespace nuX_Utils

#endif // NUX_ROOTS_HXX
