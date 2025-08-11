#ifndef NUX_POW_HXX
#define NUX_POW_HXX

#include <cctk.h> // for CCTK_HOST / CCTK_DEVICE if mapped; harmless otherwise

namespace nuX_Utils {

#define UTILS_POW(T)                                                           \
template<int N>                                                                \
CCTK_HOST CCTK_DEVICE inline T pow(T x) {                                      \
  return T(x * pow<N-1>(x));                                                   \
}                                                                              \
template<>                                                                     \
CCTK_HOST CCTK_DEVICE inline T pow<1>(T x) {                                   \
  return T(x);                                                                 \
}                                                                              \
template<>                                                                     \
CCTK_HOST CCTK_DEVICE inline T pow<0>(T) {                                     \
  return T(1);                                                                 \
}

UTILS_POW(int)
UTILS_POW(unsigned)
UTILS_POW(float)
UTILS_POW(double)

#undef UTILS_POW

} // namespace nuX_Utils

#endif

