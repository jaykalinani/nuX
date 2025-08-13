#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "nuX_fakerates.hxx"

namespace nuX_FakeRates {

FakeRatesDef global_fakerates;

extern "C" void nuX_FakeRates_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_FakeRates_Setup;
  DECLARE_CCTK_PARAMETERS;

  global_fakerates.init();

}

} // namespace
