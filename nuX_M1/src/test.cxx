#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "test.hxx"

extern "C" void nuX_M1_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (!unit_test) {
    CCTK_INFO("Skipping nuX_M1 unit tests");
    return;
  }

  std::mt19937_64 engine{nuX_M1_Tests::random_seed};
  const int repetitions = unit_test_repetitions;

  CCTK_VINFO("Running nuX_M1 unit tests with %d repetitions", repetitions);
  nuX_M1_Tests::test_packing_and_flux_helpers(engine, repetitions);
  nuX_M1_Tests::test_closure(engine, repetitions);
  nuX_M1_Tests::test_roots(engine, repetitions);
  nuX_M1_Tests::test_sources(engine, repetitions);
  CCTK_INFO("nuX_M1 unit tests passed");
}
