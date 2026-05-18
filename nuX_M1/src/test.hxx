#ifndef NUX_M1_TEST_HXX
#define NUX_M1_TEST_HXX

#include <cmath>
#include <cctk.h>
#include <limits>
#include <random>

namespace nuX_M1_Tests {

constexpr int random_seed = 1847;

template <typename T>
bool isapprox(const T x, const T y, const T atol = T(0),
              const T rtol = std::sqrt(std::numeric_limits<T>::epsilon())) {
  using std::abs;
  using std::max;
  return abs(x - y) <= max(atol, rtol * max(abs(x), abs(y)));
}

inline void require(const bool condition, const char *message) {
  if (!condition) {
    CCTK_VERROR("nuX_M1 unit test failed: %s", message);
  }
}

inline void require_close(const CCTK_REAL actual, const CCTK_REAL expected,
                          const CCTK_REAL atol, const char *message) {
  if (!isapprox(actual, expected, atol)) {
    CCTK_VERROR("nuX_M1 unit test failed: %s; expected %.17g got %.17g",
                message, static_cast<double>(expected),
                static_cast<double>(actual));
  }
}

void test_closure(std::mt19937_64 &engine, int repetitions);
void test_packing_and_flux_helpers(std::mt19937_64 &engine, int repetitions);
void test_roots(std::mt19937_64 &engine, int repetitions);
void test_sources(std::mt19937_64 &engine, int repetitions);

} // namespace nuX_M1_Tests

#endif // NUX_M1_TEST_HXX
