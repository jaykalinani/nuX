#include <cmath>
#include <loop_device.hxx>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "nuX_tensor.hxx"
#include "setup_eos.hxx"

namespace nuX_Seeds {

using namespace Loop;
using namespace EOSX;
using namespace nuX_Utils;

CCTK_DEVICE CCTK_HOST inline int kerr_old_to_comp_spatial(const int old_axis,
                                                          const int axis_mode) {
  if (axis_mode == 0)
    return old_axis;

  switch (old_axis) {
  case 1:
    return 2;
  case 2:
    return 3;
  case 3:
    return 1;
  default:
    return old_axis;
  }
}

CCTK_DEVICE CCTK_HOST inline int kerr_comp_to_old_spatial(const int comp_axis,
                                                          const int axis_mode) {
  if (axis_mode == 0)
    return comp_axis;

  switch (comp_axis) {
  case 1:
    return 3;
  case 2:
    return 1;
  case 3:
    return 2;
  default:
    return comp_axis;
  }
}

CCTK_DEVICE CCTK_HOST inline CCTK_REAL
kerr_coord_in_old_axes(const PointDesc &p, const int old_axis,
                       const int axis_mode) {
  if (axis_mode == 0) {
    switch (old_axis) {
    case 1:
      return p.x;
    case 2:
      return p.y;
    case 3:
      return p.z;
    default:
      return 0.0;
    }
  }

  switch (old_axis) {
  case 1:
    return p.y;
  case 2:
    return p.z;
  case 3:
    return p.x;
  default:
    return 0.0;
  }
}

CCTK_DEVICE CCTK_HOST inline void
kerr_beam_state(const tensor::slicing_geometry_const &geom, const PointDesc &p,
                const int axis_mode, const CCTK_REAL kerr_beam_width,
                const CCTK_REAL kerr_beam_position, CCTK_REAL &E_new,
                CCTK_REAL &Fx_new, CCTK_REAL &Fy_new, CCTK_REAL &Fz_new,
                CCTK_REAL &N_new) {
  E_new = 0.0;
  Fx_new = 0.0;
  Fy_new = 0.0;
  Fz_new = 0.0;
  N_new = 0.0;

  const CCTK_REAL y_old = kerr_coord_in_old_axes(p, 2, axis_mode);
  const CCTK_REAL z_old = kerr_coord_in_old_axes(p, 3, axis_mode);

  const bool in_beam = (fabs(y_old) <= kerr_beam_width) &&
                       (z_old >= kerr_beam_position) &&
                       (z_old <= kerr_beam_position + kerr_beam_width);
  if (!in_beam)
    return;

  CCTK_REAL const eps_flux = 0.01;
  tensor::metric<4> g_comp;
  tensor::metric<4> g_old;
  tensor::generic<CCTK_REAL, 4, 1> beta_u_comp;
  tensor::generic<CCTK_REAL, 4, 1> beta_u_old;
  tensor::generic<CCTK_REAL, 4, 1> beta_d_old;
  tensor::generic<CCTK_REAL, 4, 1> F_u_old;
  tensor::generic<CCTK_REAL, 4, 1> F_d_old;

  geom.get_metric(p, &g_comp);
  geom.get_shift_vec(p, &beta_u_comp);

  g_old(0, 0) = g_comp(0, 0);
  for (int i_old = 1; i_old <= 3; ++i_old) {
    const int i_comp = kerr_old_to_comp_spatial(i_old, axis_mode);
    g_old(0, i_old) = g_comp(0, i_comp);
    g_old(i_old, 0) = g_comp(i_comp, 0);
    beta_u_old(i_old) = beta_u_comp(i_comp);
    for (int j_old = 1; j_old <= 3; ++j_old) {
      const int j_comp = kerr_old_to_comp_spatial(j_old, axis_mode);
      g_old(i_old, j_old) = g_comp(i_comp, j_comp);
    }
  }
  beta_u_old(0) = 0.0;
  tensor::contract(g_old, beta_u_old, &beta_d_old);

  CCTK_REAL const alp_ijk = geom.get_lapse(p);
  CCTK_REAL const g_xx = g_old(1, 1);
  CCTK_REAL const beta_x = beta_d_old(1);
  CCTK_REAL const beta2 = tensor::dot(beta_u_old, beta_d_old);

  CCTK_REAL const disc =
      beta_x * beta_x - beta2 + alp_ijk * alp_ijk * (1.0 - eps_flux);
  CCTK_REAL a = 0.0;
  if (g_xx > 0.0) {
    a = (-beta_x + sqrt(fmax(disc, 0.0))) / g_xx;
  }

  E_new = sqrt(fmax(metric::spatial_det(g_old(1, 1), g_old(1, 2), g_old(1, 3),
                                        g_old(2, 2), g_old(2, 3), g_old(3, 3)),
                    0.0));

  CCTK_REAL const inv_alp = (alp_ijk > 0.0) ? (1.0 / alp_ijk) : 0.0;
  F_u_old(0) = 0.0;
  F_u_old(1) = (a + beta_u_old(1)) * E_new * inv_alp;
  F_u_old(2) = beta_u_old(2) * E_new * inv_alp;
  F_u_old(3) = beta_u_old(3) * E_new * inv_alp;

  tensor::contract(g_old, F_u_old, &F_d_old);
  Fx_new = F_d_old(kerr_comp_to_old_spatial(1, axis_mode));
  Fy_new = F_d_old(kerr_comp_to_old_spatial(2, axis_mode));
  Fz_new = F_d_old(kerr_comp_to_old_spatial(3, axis_mode));
  N_new = 1.0;
}

extern "C" void nuX_Seeds_SetupHydroTest_kerrschild(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupHydroTest_kerrschild;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupHydroTest_kerrschild");

  auto eos_3p_ig = global_eos_3p_ig;
  if (not CCTK_EQUALS(evolution_eos, "IdealGas")) {
    CCTK_VERROR("Invalid evolution EOS type '%s'. Please, set "
                "EOSX::evolution_eos = \"IdealGas\" in your parameter file.",
                evolution_eos);
  }

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout3(cctkGH, {1, 0, 0});
  const GF3D2layout layout4(cctkGH, {0, 1, 0});
  const GF3D2layout layout5(cctkGH, {0, 0, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p) {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        rho[ijk] = static_rho;
        eps[ijk] = static_eps;
        velx[ijk] = static_velx;
        vely[ijk] = static_vely;
        velz[ijk] = static_velz;
        Ye[ijk] = static_ye;

        press[ijk] =
            eos_3p_ig->press_from_valid_rho_eps_ye(rho[ijk], eps[ijk], Ye[ijk]);
      });
  grid.loop_all_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout3.linear(p.i, p.j, p.k);
        Avec_x[ijk] = 0.;
      });

  grid.loop_all_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout4.linear(p.i, p.j, p.k);
        Avec_y[ijk] = 0.;
      });

  grid.loop_all_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout5.linear(p.i, p.j, p.k);
        Avec_z[ijk] = 0.;
      });
}

extern "C" void nuX_Seeds_SetupNeutTest_kerrschild(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_SetupNeutTest_kerrschild;
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
    CCTK_INFO("nuX_Seeds_SetupNeutTest_kerrschild");

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          int const i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = 0.0;
          rFx[i4D] = 0.0;
          rFy[i4D] = 0.0;
          rFz[i4D] = 0.0;
          rN[i4D] = 0.0;
        }
      });
}

extern "C" void nuX_Seeds_KerrInitialBeam(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_KerrInitialBeam;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  const bool rotate_yzx = CCTK_EQUALS(kerr_axis_permutation, "yzx");
  const int axis_mode = rotate_yzx ? 1 : 0;
  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);

  // THC seeds the first physical cell at the inflow face before the initial
  // closure. Runtime KerrInflow below maintains both ghost cells and the active
  // physical boundary rows that THC overwrites through its ghost-zone loops.
  const int inflow_axis = rotate_yzx ? 1 : 0;
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (p.BI[inflow_axis] != -1)
          return;

        CCTK_REAL E_new = 0.0;
        CCTK_REAL Fx_new = 0.0;
        CCTK_REAL Fy_new = 0.0;
        CCTK_REAL Fz_new = 0.0;
        CCTK_REAL N_new = 0.0;
        kerr_beam_state(geom, p, axis_mode, kerr_beam_width, kerr_beam_position,
                        E_new, Fx_new, Fy_new, Fz_new, N_new);

        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rE[i4D] = E_new;
          rFx[i4D] = Fx_new;
          rFy[i4D] = Fy_new;
          rFz[i4D] = Fz_new;
          rN[i4D] = N_new;
        }
      });
}

extern "C" void nuX_Seeds_KerrInflow(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_KerrInflow;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const GF3D2layout layout_vc(cctkGH, {0, 0, 0});
  const int nghy = cctk_nghostzones[1];
  const int nghz = cctk_nghostzones[2];
  const bool rotate_yzx = CCTK_EQUALS(kerr_axis_permutation, "yzx");
  const int axis_mode = rotate_yzx ? 1 : 0;
  tensor::slicing_geometry_const geom(layout_vc, layout_cc, alp, betax, betay,
                                      betaz, gxx, gxy, gxz, gyy, gyz, gzz, kxx,
                                      kxy, kxz, kyy, kyz, kzz);

  // Apply faces in THC order to get deterministic edge/corner precedence:
  // -X, +X, -Y, +Y, -Z, +Z.
  // CarpetX p.NI marks ghost/boundary points, while p.BI marks the physical
  // outermost interior point on true domain boundaries. The Kerr inflow needs
  // both: ghost cells for stencils, plus physical boundary rows on the faces
  // that THC overwrites through its cctk_nghostzones loops.
  //
  // This must be split into multiple passes (instead of one combined kernel):
  // Y-face copy BCs read neighboring cells, and doing that while writing all
  // faces in one kernel causes boundary read/write races on GPU.

  if (!rotate_yzx) {
    grid.loop_bnd_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[0] != -1)
            return;

          CCTK_REAL E_new = 0.0;
          CCTK_REAL Fx_new = 0.0;
          CCTK_REAL Fy_new = 0.0;
          CCTK_REAL Fz_new = 0.0;
          CCTK_REAL N_new = 0.0;
          kerr_beam_state(geom, p, axis_mode, kerr_beam_width,
                          kerr_beam_position, E_new, Fx_new, Fy_new, Fz_new,
                          N_new);

          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = E_new;
            rFx[i4D] = Fx_new;
            rFy[i4D] = Fy_new;
            rFz[i4D] = Fz_new;
            rN[i4D] = N_new;
          }
        });

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[0] != -1)
            return;

          CCTK_REAL E_new = 0.0;
          CCTK_REAL Fx_new = 0.0;
          CCTK_REAL Fy_new = 0.0;
          CCTK_REAL Fz_new = 0.0;
          CCTK_REAL N_new = 0.0;
          kerr_beam_state(geom, p, axis_mode, kerr_beam_width,
                          kerr_beam_position, E_new, Fx_new, Fy_new, Fz_new,
                          N_new);

          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = E_new;
            rFx[i4D] = Fx_new;
            rFy[i4D] = Fy_new;
            rFz[i4D] = Fz_new;
            rN[i4D] = N_new;
          }
        });

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.BI[0] != -1)
            return;

          CCTK_REAL E_new = 0.0;
          CCTK_REAL Fx_new = 0.0;
          CCTK_REAL Fy_new = 0.0;
          CCTK_REAL Fz_new = 0.0;
          CCTK_REAL N_new = 0.0;
          kerr_beam_state(geom, p, axis_mode, kerr_beam_width,
                          kerr_beam_position, E_new, Fx_new, Fy_new, Fz_new,
                          N_new);

          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = E_new;
            rFx[i4D] = Fx_new;
            rFy[i4D] = Fy_new;
            rFz[i4D] = Fz_new;
            rN[i4D] = N_new;
          }
        });

    // +X: vacuum
    grid.loop_bnd_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[0] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[0] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.BI[0] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    // -Y: copy from first interior Y plane.
    const int jsrc_lo = nghy;
    for (int d = 0; d < nghy; ++d) {
      const int jdst = d;
      grid.loop_all_device<1, 1, 1>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.NI[1] != -1 || p.j != jdst)
              return;
            for (int ig = 0; ig < ngroups * nspecies; ++ig) {
              const int i4Db = layout_cc.linear(p.i, jdst, p.k, ig);
              const int i4Di = layout_cc.linear(p.i, jsrc_lo, p.k, ig);
              rE[i4Db] = rE[i4Di];
              rFx[i4Db] = rFx[i4Di];
              rFy[i4Db] = rFy[i4Di];
              rFz[i4Db] = rFz[i4Di];
              rN[i4Db] = rN[i4Di];
            }
          });
    }

    // +Y: copy from last interior Y plane.
    const int jsrc_hi = cctk_lsh[1] - nghy - 1;
    for (int d = 0; d < nghy; ++d) {
      const int jdst = cctk_lsh[1] - nghy + d;
      grid.loop_all_device<1, 1, 1>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.NI[1] != 1 || p.j != jdst)
              return;
            for (int ig = 0; ig < ngroups * nspecies; ++ig) {
              const int i4Db = layout_cc.linear(p.i, jdst, p.k, ig);
              const int i4Di = layout_cc.linear(p.i, jsrc_hi, p.k, ig);
              rE[i4Db] = rE[i4Di];
              rFx[i4Db] = rFx[i4Di];
              rFy[i4Db] = rFy[i4Di];
              rFz[i4Db] = rFz[i4Di];
              rN[i4Db] = rN[i4Di];
            }
          });
    }

    // -Z: vacuum
    grid.loop_bnd_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[2] != -1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[2] != -1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.BI[2] != -1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    // +Z: vacuum
    grid.loop_bnd_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[2] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[2] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.BI[2] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });
  } else {
    // Rotated branch for x->y, y->z, z->x. Keep the same boundary semantics
    // as the baseline branch, but permute which computational face gets each
    // role.
    grid.loop_bnd_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[1] != -1)
            return;

          CCTK_REAL E_new = 0.0;
          CCTK_REAL Fx_new = 0.0;
          CCTK_REAL Fy_new = 0.0;
          CCTK_REAL Fz_new = 0.0;
          CCTK_REAL N_new = 0.0;
          kerr_beam_state(geom, p, axis_mode, kerr_beam_width,
                          kerr_beam_position, E_new, Fx_new, Fy_new, Fz_new,
                          N_new);

          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = E_new;
            rFx[i4D] = Fx_new;
            rFy[i4D] = Fy_new;
            rFz[i4D] = Fz_new;
            rN[i4D] = N_new;
          }
        });

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[1] != -1)
            return;

          CCTK_REAL E_new = 0.0;
          CCTK_REAL Fx_new = 0.0;
          CCTK_REAL Fy_new = 0.0;
          CCTK_REAL Fz_new = 0.0;
          CCTK_REAL N_new = 0.0;
          kerr_beam_state(geom, p, axis_mode, kerr_beam_width,
                          kerr_beam_position, E_new, Fx_new, Fy_new, Fz_new,
                          N_new);

          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = E_new;
            rFx[i4D] = Fx_new;
            rFy[i4D] = Fy_new;
            rFz[i4D] = Fz_new;
            rN[i4D] = N_new;
          }
        });

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.BI[1] != -1)
            return;

          CCTK_REAL E_new = 0.0;
          CCTK_REAL Fx_new = 0.0;
          CCTK_REAL Fy_new = 0.0;
          CCTK_REAL Fz_new = 0.0;
          CCTK_REAL N_new = 0.0;
          kerr_beam_state(geom, p, axis_mode, kerr_beam_width,
                          kerr_beam_position, E_new, Fx_new, Fy_new, Fz_new,
                          N_new);

          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = E_new;
            rFx[i4D] = Fx_new;
            rFy[i4D] = Fy_new;
            rFz[i4D] = Fz_new;
            rN[i4D] = N_new;
          }
        });

    // +Y: vacuum
    grid.loop_bnd_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[1] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[1] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.BI[1] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    // -Z: copy from first interior Z plane.
    const int ksrc_lo = nghz;
    for (int d = 0; d < nghz; ++d) {
      const int kdst = d;
      grid.loop_all_device<1, 1, 1>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.NI[2] != -1 || p.k != kdst)
              return;
            for (int ig = 0; ig < ngroups * nspecies; ++ig) {
              const int i4Db = layout_cc.linear(p.i, p.j, kdst, ig);
              const int i4Di = layout_cc.linear(p.i, p.j, ksrc_lo, ig);
              rE[i4Db] = rE[i4Di];
              rFx[i4Db] = rFx[i4Di];
              rFy[i4Db] = rFy[i4Di];
              rFz[i4Db] = rFz[i4Di];
              rN[i4Db] = rN[i4Di];
            }
          });
    }

    // +Z: copy from last interior Z plane.
    const int ksrc_hi = cctk_lsh[2] - nghz - 1;
    for (int d = 0; d < nghz; ++d) {
      const int kdst = cctk_lsh[2] - nghz + d;
      grid.loop_all_device<1, 1, 1>(
          grid.nghostzones,
          [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.NI[2] != 1 || p.k != kdst)
              return;
            for (int ig = 0; ig < ngroups * nspecies; ++ig) {
              const int i4Db = layout_cc.linear(p.i, p.j, kdst, ig);
              const int i4Di = layout_cc.linear(p.i, p.j, ksrc_hi, ig);
              rE[i4Db] = rE[i4Di];
              rFx[i4Db] = rFx[i4Di];
              rFy[i4Db] = rFy[i4Di];
              rFz[i4Db] = rFz[i4Di];
              rN[i4Db] = rN[i4Di];
            }
          });
    }

    // -X: vacuum
    grid.loop_bnd_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[0] != -1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[0] != -1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.BI[0] != -1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    // +X: vacuum
    grid.loop_bnd_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[0] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.NI[0] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.BI[0] != 1)
            return;
          for (int ig = 0; ig < ngroups * nspecies; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
            rN[i4D] = 0.0;
          }
        });
  }
}

extern "C" void nuX_Seeds_StorePostBCState(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_StorePostBCState;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        for (int ig = 0; ig < ngroups * nspecies; ++ig) {
          const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
          rN_bc_snapshot[i4D] = rN[i4D];
          rE_bc_snapshot[i4D] = rE[i4D];
          rFx_bc_snapshot[i4D] = rFx[i4D];
          rFy_bc_snapshot[i4D] = rFy[i4D];
          rFz_bc_snapshot[i4D] = rFz[i4D];
        }
      });
}

extern "C" void nuX_Seeds_KerrSchild_Mask(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_nuX_Seeds_KerrSchild_Mask;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  const GF3D2layout layout_cc(cctkGH, {1, 1, 1});
  const CCTK_REAL mask_r2 = kerr_mask_radius * kerr_mask_radius;

  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const int ijk = layout_cc.linear(p.i, p.j, p.k);
        const bool masked = (p.x * p.x + p.y * p.y + p.z * p.z) < mask_r2;
        nuX_m1_mask[ijk] = masked ? 1.0 : 0.0;

        if (masked) {
          for (int ig = 0; ig < nspecies * ngroups; ++ig) {
            const int i4D = layout_cc.linear(p.i, p.j, p.k, ig);
            rN[i4D] = 0.0;
            rE[i4D] = 0.0;
            rFx[i4D] = 0.0;
            rFy[i4D] = 0.0;
            rFz[i4D] = 0.0;
          }
        }
      });
}

} // namespace nuX_Seeds
