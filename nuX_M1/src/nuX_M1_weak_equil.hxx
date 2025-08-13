#include "setup_eos.hxx"

#define MAX_SPECIES 3

namespace nuX_M1 {

// Set constants
const CCTK_REAL nu_2DNR_eps_lim  = 1.e-7;
const int nu_2DNR_n_max    = 100;
const int nu_bis_n_cut_max = 8;

// Neutrino equilibrium physical constants
const CCTK_REAL hc_mevfm = 1.23984172e3;           // hc    [MeV fm] (not reduced)
const CCTK_REAL pi       = 3.14159265358979323846; // pi    [-]
const CCTK_REAL pi2      = pi*pi;                  // pi**2 [-]
const CCTK_REAL pi4      = pi2*pi2;                // pi**4 [-]

const CCTK_REAL nu_n_prefactor = 4.0/3.0*pi/(hc_mevfm*hc_mevfm*hc_mevfm); // 4/3 *pi/(hc)**3 [1/MeV^3/fm^3]
const CCTK_REAL nu_e_prefactor = 4.0*pi/(hc_mevfm*hc_mevfm*hc_mevfm);     // 4*pi/(hc)**3    [1/MeV^3 fm^3]

const CCTK_REAL nu_7pi4_60 = 7.0*pi4/60.0;  // 7*pi**4/60  [-]
const CCTK_REAL nu_7pi4_30 = 7.0*pi4/30.0;  // 7*pi**4/30  [-]
const CCTK_REAL nu_7pi4_15 = 7.0*pi4/15.0;  // 7*pi**4/15  [-]
const CCTK_REAL nu_14pi4_15 = 14.0*pi4/15.0; // 14*pi**4/15 [-]
  
/// Low level function for neutrino equilibrium, not intended for outside use
CCTK_DEVICE inline CCTK_REAL error_func_eq_weak(CCTK_REAL Yle, CCTK_REAL e_eq, CCTK_REAL y[2]) {
  CCTK_REAL err = abs(y[0]/Yle) + abs(y[1]/1.0);
  return err;
}

CCTK_DEVICE inline void inv_jacobi(CCTK_REAL det, CCTK_REAL J[2][2], CCTK_REAL invJ[2][2]) {
  CCTK_REAL inv_det = 1.0/det;
  invJ[0][0] =  J[1][1]*inv_det;
  invJ[1][1] =  J[0][0]*inv_det;
  invJ[0][1] = -J[0][1]*inv_det;
  invJ[1][0] = -J[1][0]*inv_det;
}

template<typename EOSType>
CCTK_DEVICE inline int eta_e_gradient(CCTK_REAL rho, CCTK_REAL T, CCTK_REAL *Y, CCTK_REAL eta, CCTK_REAL &deta_dT, CCTK_REAL &deta_dYe, CCTK_REAL &de_dT, CCTK_REAL &de_dYe, const EOSType* tabeos) {
  int ierr=1;

  const CCTK_REAL min_T = tabeos->rgtemp.min;
  const CCTK_REAL max_T = tabeos->rgtemp.max;
  const CCTK_REAL min_Y = tabeos->rgye.min;
  const CCTK_REAL max_Y = tabeos->rgye.max;

  const CCTK_REAL Ye_delta = 0.005;
  const CCTK_REAL T_delta = 0.01;

  CCTK_REAL Y1[MAX_SPECIES] = {0.0};
  CCTK_REAL Y2[MAX_SPECIES] = {0.0};

  Y1[0] = fmax(Y[0] - Ye_delta, min_Y);
  CCTK_REAL mu_l1 = tabeos->mu_lepton_from_valid_rho_temp_ye(rho, T, Y1[0]);
  CCTK_REAL e1 = tabeos->eps_from_valid_rho_temp_ye(rho, T, Y1[0]);
  
  Y2[0] = fmin(Y[0] + Ye_delta, max_Y);
  CCTK_REAL mu_l2 = tabeos->mu_lepton_from_valid_rho_temp_ye(rho, T, Y2[0]);
  CCTK_REAL e2 = tabeos->eps_from_valid_rho_temp_ye(rho, T, Y2[0]);

  CCTK_REAL dmu_l_dYe = (mu_l2-mu_l1)/(Y2[0] - Y1[0]);
  de_dYe         = (e2-e1)/(Y2[0] - Y1[0]);

  CCTK_REAL T1 = fmax(T - T_delta, min_T);
  mu_l1 = tabeos->mu_lepton_from_valid_rho_temp_ye(rho, T, Y1[0]);
  e1 = tabeos->eps_from_valid_rho_temp_ye(rho, T, Y1[0]);

  CCTK_REAL T2 = fmin(T + T_delta, max_T);
  mu_l2 = tabeos->mu_lepton_from_valid_rho_temp_ye(rho, T, Y2[0]);
  e2 = tabeos->eps_from_valid_rho_temp_ye(rho, T, Y2[0]);

  CCTK_REAL dmu_l_dT   = (mu_l2 - mu_l1)/(T2 - T1);
  de_dT          = (e2 - e1)/(T2 - T1);
  
  deta_dT  = (dmu_l_dT - eta )/T; // [1/MeV] TODO: Check
  deta_dYe = dmu_l_dYe/T;      // [-]

  if (isnan(deta_dT)||isnan(deta_dYe)||isnan(de_dT)||isnan(de_dYe)) {
    ierr = 1;
  } else {
    ierr = 0;
  }

  return ierr;
}

template<typename EOSType>
CCTK_DEVICE inline int jacobi_eq_weak(CCTK_REAL rho, CCTK_REAL n, CCTK_REAL e_eq, CCTK_REAL Yle, CCTK_REAL x[2], CCTK_REAL J[2][2], const EOSType* tabeos) {
  int ierr = 0;

  CCTK_REAL T = x[0];
  CCTK_REAL Y[MAX_SPECIES] = {0.0};
  Y[0] = x[1];

  if (isnan(T)) {
    ierr = 1;
    return ierr;
  }

  CCTK_REAL mu_l = tabeos->mu_lepton_from_valid_rho_temp_ye(rho, T, Y[0]);
  CCTK_REAL eta = mu_l/T;
  CCTK_REAL eta2 = eta*eta;

  if (isnan(eta)) {
    ierr = 1;
    return ierr;
  }

  CCTK_REAL detadt,detadye,dedt,dedye;
  ierr = eta_e_gradient(rho,T,Y,eta,detadt,detadye,dedt,dedye,tabeos);
  if (ierr != 0){
    return ierr;
  }

  CCTK_REAL T2 = T*T;
  CCTK_REAL T3 = T2*T;
  // CCTK_REAL T4 = T3*T;

  J[0][0] = nu_n_prefactor/n*T2*(3.e0*eta*(pi2+eta2)+T*(pi2+3.e0*eta2)*detadt);
  J[0][1] = 1.e0+nu_n_prefactor/n*T3*(pi2+3.e0*eta2)*detadye;

  J[1][0] = (dedt+nu_e_prefactor*T3*(nu_7pi4_15+nu_14pi4_15+2.e0*eta2*(pi2+0.5*eta2)+eta*T*(pi2+eta2)*detadt))/e_eq;
  J[1][1] = (dedye+nu_e_prefactor*T3*eta*(pi2+eta2)*detadye)/e_eq;

  return ierr;
}

template<typename EOSType>
CCTK_DEVICE inline void func_eq_weak(CCTK_REAL rho, CCTK_REAL n, CCTK_REAL e_eq, CCTK_REAL Yle, CCTK_REAL x[2], CCTK_REAL y[2], const EOSType* tabeos) {
  CCTK_REAL T = x[0];

  CCTK_REAL Y[MAX_SPECIES] = {0.0};
  Y[0] = x[1];

  CCTK_REAL mu_l = tabeos->mu_lepton_from_valid_rho_temp_ye(rho, T, Y[0]);
  CCTK_REAL e = tabeos->eps_from_valid_rho_temp_ye(rho, T, Y[0]);
  CCTK_REAL eta = mu_l/T;
  CCTK_REAL eta2 = eta*eta;

  CCTK_REAL t3 = T*T*T;
  CCTK_REAL t4 = t3*T;
  y[0] = Y[0] + nu_n_prefactor*t3*eta*(pi2 + eta2)/n - Yle;
  y[1] = (e+nu_e_prefactor*t4*((nu_7pi4_60+0.5*eta2*(pi2+0.5*eta2))+nu_7pi4_30))/e_eq - 1.0;

  return;
}

template<typename EOSType>
CCTK_DEVICE inline int trapped_equilibrium_2DNR(CCTK_REAL rho, CCTK_REAL n, CCTK_REAL e, CCTK_REAL Yle, CCTK_REAL x0[2], CCTK_REAL x1[2], const EOSType* tabeos) {
  int ierr = 1;

  const CCTK_REAL min_T = tabeos->rgtemp.min;
  const CCTK_REAL max_T = tabeos->rgtemp.max;
  const CCTK_REAL min_Y = tabeos->rgye.min;
  const CCTK_REAL max_Y = tabeos->rgye.max;

  // initialize the solution
  x1[0] = x0[0];
  x1[1] = x0[1];
  bool KKT = false;

  //compute the initial residuals
  CCTK_REAL y[2] = {0.0};
  func_eq_weak(rho,n,e,Yle,x1,y, tabeos);

  // compute the error from the residuals
  CCTK_REAL err = error_func_eq_weak(Yle,e,y);

  // initialize the iteration variables
  int n_iter = 0;
  CCTK_REAL J[2][2] = {0.0};
  CCTK_REAL invJ[2][2] = {0.0};
  CCTK_REAL dx1[2] = {0.0};
  CCTK_REAL dxa[2] = {0.0};
  CCTK_REAL norm[2] = {0.0};
  CCTK_REAL x1_tmp[2] = {0.0};

  // loop until a low enough residual is found or until  a too
  // large number of steps has been performed
  while (err>nu_2DNR_eps_lim && n_iter<=nu_2DNR_n_max && !KKT) {
    // compute the Jacobian
    ierr = jacobi_eq_weak(rho,n,e,Yle,x1,J,tabeos);
    if (ierr != 0) {
      return ierr;
    }

    // compute and check the determinant of the Jacobian
    CCTK_REAL det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    if (det==0.0) {
      ierr = 1;
      return ierr;
    }

    // invert the Jacobian
    inv_jacobi(det,J,invJ);

    // compute the next step
    dx1[0] = - (invJ[0][0]*y[0] + invJ[0][1]*y[1]);
    dx1[1] = - (invJ[1][0]*y[0] + invJ[1][1]*y[1]);

    // check if we are the boundary of the table
    if (x1[0] == min_T) {
      norm[0] = -1.0;
    } else if (x1[0] == max_T) {
      norm[0] = 1.0;
    } else { 
      norm[0] = 0.0;
    }

    if (x1[1] == min_Y) {
      norm[1] = -1.0;
    } else if (x1[1] == max_Y) {
      norm[1] = 1.0;
    } else {
      norm[1] = 0.0;
    }

    // Take the part of the gradient that is active (pointing within the eos domain)
    CCTK_REAL scal = norm[0]*norm[0] + norm[1]*norm[1];
    if (scal <= 0.5) { // this can only happen if norm = (0, 0)
      scal = 1.0;
    }
    dxa[0] = dx1[0] - (dx1[0]*norm[0] + dx1[1]*norm[1])*norm[0]/scal;
    dxa[1] = dx1[1] - (dx1[0]*norm[0] + dx1[1]*norm[1])*norm[1]/scal;

    if ((dxa[0]*dxa[0] + dxa[1]*dxa[1]) < (nu_2DNR_eps_lim*nu_2DNR_eps_lim * (dx1[0]*dx1[0] + dx1[1]*dx1[1]))) {
      KKT = true;
      ierr = 2;
      return ierr;
    }

    int n_cut = 0;
    CCTK_REAL fac_cut = 1.0;
    CCTK_REAL err_old = err;

    while (n_cut <= nu_bis_n_cut_max && err >= err_old) {
      // the variation of x1 is divided by an powers of 2 if the
      // error is not decreasing along the gradient direction
      
      x1_tmp[0] = x1[0] + (dx1[0]*fac_cut);
      x1_tmp[1] = x1[1] + (dx1[1]*fac_cut);

      // check if the next step calculation had problems
      if (isnan(x1_tmp[0])) {
        ierr = 1;
        return ierr;
      }

      // tabBoundsFlag = enforceTableBounds(rho, x1_tmp[0], x1_tmp[1]);
      x1_tmp[0] = fmin(fmax(x1_tmp[0],min_T),max_T);
      x1_tmp[1] = fmin(fmax(x1_tmp[1],min_Y),max_Y);

      // assign the new point
      x1[0] = x1_tmp[0];
      x1[1] = x1_tmp[1];

      // compute the residuals for the new point
      func_eq_weak(rho,n,e,Yle,x1,y,tabeos);

      // compute the error
      err = error_func_eq_weak(Yle,e,y);

      // update the bisection cut along the gradient
      n_cut += 1;
      fac_cut *= 0.5;
    }

    // update the iteration
    n_iter += 1;
  }
    
  if (n_iter <= nu_2DNR_n_max) {
    ierr = 0;
  } else {
    ierr = 1;
  }
  
  return ierr;
}

/// Master Function: Calculate hot (neutrino trapped) beta equilibrium T_eq and Y_eq given rho, e, and Yl
template<typename EOSType>
CCTK_DEVICE inline int BetaEquilibriumTrapped(CCTK_REAL rho, CCTK_REAL n, CCTK_REAL e, CCTK_REAL Yl, CCTK_REAL &T_eq, CCTK_REAL &Y_eq, CCTK_REAL T_guess, CCTK_REAL Y_guess, const EOSType* tabeos) {
  const int n_at = 16;
  CCTK_REAL vec_guess[n_at][2] = { 
    {1.00e0, 1.00e0},
    {0.90e0, 1.25e0},
    {0.90e0, 1.10e0},
    {0.90e0, 1.00e0},
    {0.90e0, 0.90e0},
    {0.90e0, 0.75e0},
    {0.75e0, 1.25e0},
    {0.75e0, 1.10e0},
    {0.75e0, 1.00e0},
    {0.75e0, 0.90e0},
    {0.75e0, 0.75e0},
    {0.50e0, 1.25e0},
    {0.50e0, 1.10e0},
    {0.50e0, 1.00e0},
    {0.50e0, 0.90e0},
    {0.50e0, 0.75e0},
  };

  // ierr = 0    Equilibrium found
  // ierr = 1    Equilibrium not found
  int ierr = 1;
  int na = 0; // counter for the number of attempts

  CCTK_REAL x0[2], x1[2]; // T,Ye guess and T,Ye result

  while (ierr!=0 && na<n_at) {
    x0[0] = vec_guess[na][0] * T_guess;
    x0[1] = vec_guess[na][1] * Y_guess;

    ierr = trapped_equilibrium_2DNR(rho, n, e, Yl, x0, x1, tabeos);

    na += 1;
  }

  if (ierr==0){ // Success
    T_eq = x1[0];
    Y_eq = x1[1];
  } else {      // Failure
    T_eq = T_guess;       // Set results to guesses
    Y_eq = Y_guess;
  }

  return ierr;
}

}; // namespace
