#include <cstdio>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "../../../CarpetX/CarpetX/src/fillpatch.hxx"
#include "../../../CarpetX/CarpetX/src/schedule.hxx"
#include "../../../CarpetX/CarpetX/src/task_manager.hxx"

namespace nuX_M1 {

using namespace CarpetX;

extern "C" void nuX_M1_Sync(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_nuX_M1_Sync;
    DECLARE_CCTK_PARAMETERS;

    if(verbose) {
        CCTK_INFO("nuX_M1_Sync");
    }

    // Do nothing, let SYNC handle BCs.
}

void ApplyOuterBC(CCTK_ARGUMENTS, const std::vector<int> &groups) {
  task_manager tasks1;
  task_manager tasks2;

  for (const int gi : groups) {
    active_levels->loop_serially([&](auto &restrict leveldata) {
      auto &restrict groupdata = *leveldata.groupdata.at(gi);

      const int ntls = groupdata.mfab.size();
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;

      // Copy from adjacent boxes on same level and apply boundary conditions
      // Even though this introduces additional communication, it makes the code
      // compatible with symmetries.
      for (int tl = 0; tl < sync_tl; ++tl) {
        tasks1.submit_serially([&tasks2, &leveldata, &groupdata, tl]() {
          FillPatch_Sync(tasks2, groupdata, *groupdata.mfab.at(tl),
                         ghext->patchdata.at(leveldata.patch)
                             .amrcore->Geom(leveldata.level));
        });
      } // for tl
    });
  } // for gi

  tasks1.run_tasks_serially();
  synchronize();
  tasks2.run_tasks_serially();
  synchronize();

  assert(ghext->num_patches() == 1);
}

extern "C" void nuX_M1_ApplyOuterBC(CCTK_ARGUMENTS) {
  static const std::vector<int> groups = {
      CCTK_GroupIndex("nuX_M1::rN"),
      CCTK_GroupIndex("nuX_M1::rE"),
      CCTK_GroupIndex("nuX_M1::rF")};

  ApplyOuterBC(CCTK_PASS_CTOC, groups);
}

}
