#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

namespace nuX_M1 {

using namespace Loop;

namespace {

template <typename CopyFn>
CCTK_HOST void copy_state_snapshot(CCTK_ARGUMENTS, CopyFn copy_fn) {
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const int groupspec = ngroups * nspecies;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < groupspec; ++ig) {
          copy_fn(layout_cc.linear(p.i, p.j, p.k, ig));
        }
      });
}

} // namespace

extern "C" void nuX_M1_StorePostUpdateState(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_StorePostUpdateState;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_StorePostUpdateState");
  }

  copy_state_snapshot(CCTK_PASS_CTOC, [=] CCTK_DEVICE(const int i4D) {
    rN_postupdate_snapshot[i4D] = rN[i4D];
    rE_postupdate_snapshot[i4D] = rE[i4D];
    rFx_postupdate_snapshot[i4D] = rFx[i4D];
    rFy_postupdate_snapshot[i4D] = rFy[i4D];
    rFz_postupdate_snapshot[i4D] = rFz[i4D];
  });
}

extern "C" void nuX_M1_StorePrePostBCState(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_StorePrePostBCState;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_StorePrePostBCState");
  }

  copy_state_snapshot(CCTK_PASS_CTOC, [=] CCTK_DEVICE(const int i4D) {
    rN_prepostbc_snapshot[i4D] = rN[i4D];
    rE_prepostbc_snapshot[i4D] = rE[i4D];
    rFx_prepostbc_snapshot[i4D] = rFx[i4D];
    rFy_prepostbc_snapshot[i4D] = rFy[i4D];
    rFz_prepostbc_snapshot[i4D] = rFz[i4D];
  });
}

extern "C" void nuX_M1_StorePreFluxState(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_StorePreFluxState;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_StorePreFluxState");
  }

  copy_state_snapshot(CCTK_PASS_CTOC, [=] CCTK_DEVICE(const int i4D) {
    rN_preflux_snapshot[i4D] = rN[i4D];
    rE_preflux_snapshot[i4D] = rE[i4D];
    rFx_preflux_snapshot[i4D] = rFx[i4D];
    rFy_preflux_snapshot[i4D] = rFy[i4D];
    rFz_preflux_snapshot[i4D] = rFz[i4D];
    rPxx_preflux_snapshot[i4D] = rPxx[i4D];
    rPxy_preflux_snapshot[i4D] = rPxy[i4D];
    rPxz_preflux_snapshot[i4D] = rPxz[i4D];
    rPyy_preflux_snapshot[i4D] = rPyy[i4D];
    rPyz_preflux_snapshot[i4D] = rPyz[i4D];
    rPzz_preflux_snapshot[i4D] = rPzz[i4D];
    chi_preflux_snapshot[i4D] = chi[i4D];
  });
}

extern "C" void nuX_M1_StoreTransportInputs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_M1_StoreTransportInputs;
  DECLARE_CCTK_PARAMETERS;

  if (verbose) {
    CCTK_INFO("nuX_M1_StoreTransportInputs");
  }

  copy_state_snapshot(CCTK_PASS_CTOC, [=] CCTK_DEVICE(const int i4D) {
    rN_transport_input[i4D] = rN[i4D];
    rE_transport_input[i4D] = rE[i4D];
    rFx_transport_input[i4D] = rFx[i4D];
    rFy_transport_input[i4D] = rFy[i4D];
    rFz_transport_input[i4D] = rFz[i4D];
  });
}

} // namespace nuX_M1
