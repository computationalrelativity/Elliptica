/*
// Alireza Rashti
// November 2020
*/


/* collection of functions to deal with tuning of the physical system */

#include "sys_tune.h"

/* tuning P_ADMs to be zero. */
int sys_tune_ADM_momenta(Physics_T *const phys)
{
  FUNC_TIC
  
  const char *const par = sysGets("P_ADM_control_method");
  char *adjust[3] = {0};
  double adm_p[3] = {0};
  double adm_j[3] = {0};
  
  printf(Pretty0"method        = %s\n",par);
  printf(Pretty0"update weight = %e\n",
        sysGetd("P_ADM_control_update_weight"));
  
  parse_adjust_parameter(par,adjust);
  
  void (*P_ADM_control_0)(Physics_T *const phys) =
                              get_func_P_ADM_adjustment(adjust[0]);
  
  void (*P_ADM_control_1)(Physics_T *const phys) =
                              get_func_P_ADM_adjustment(adjust[1]);
                              
  void (*P_ADM_control_2)(Physics_T *const phys) =
                              get_func_P_ADM_adjustment(adjust[2]);
  
  
  observe(phys,"ADM(P)",Gets("Observe_ADM_P"),adm_p);
  observe(phys,"ADM(J)",Gets("Observe_ADM_J"),adm_j);
  
  printf(Pretty0"Current P_ADM = (%e,%e,%e)\n",adm_p[0],adm_p[1],adm_p[2]);
  printf(Pretty0"Current J_ADM = (%e,%e,%e)\n",adm_j[0],adm_j[1],adm_j[2]);
  
  sysSetd("Px_ADM",adm_p[0]);
  sysSetd("Py_ADM",adm_p[1]);
  sysSetd("Pz_ADM",adm_p[2]);
  
  sysSetd("Jx_ADM",adm_j[0]);
  sysSetd("Jy_ADM",adm_j[1]);
  sysSetd("Jz_ADM",adm_j[2]);
  
  if (P_ADM_control_0)
    P_ADM_control_0(phys);
  
  if (P_ADM_control_1)
    P_ADM_control_1(phys);
    
  if (P_ADM_control_2)
    P_ADM_control_2(phys);
    
  Free(adjust[0]);
  Free(adjust[1]);
  Free(adjust[2]);
  
  FUNC_TOC
  
  return EXIT_SUCCESS;
}

/* getting adjustment str, returns the relevant function. */
static fAdjustment_t *get_func_P_ADM_adjustment(const char *const adjust)
{
  fAdjustment_t *f = 0;
  
  if (!adjust)
  {
    f = 0;
  }
  else if (strcmp_i(adjust,"none"))
  {
    f = 0;
  }
  else if (strcmp_i(adjust,"x_CM"))
  {
    f = Py_ADM_is0_by_x_CM;
  }
  else if (strcmp_i(adjust,"y_CM"))
  {
    f = Px_ADM_is0_by_y_CM;
  }
  else
    Error0(NO_OPTION);
  
  return f;
}

/* parsing adjust parameter value and fill components consequently */
static void parse_adjust_parameter(const char *const par,char *adjust[3])
{
  if (!strstr_i(par,"adjust(") && !strstr_i(par,"none"))
    Errors("Syntax error for '%s'.\n",par);
  
  /* if it is none */  
  if (strcmp_i(par,"none"))
  {
    adjust[0] = 0;
    adjust[1] = 0;
    adjust[2] = 0;
    
    return;
  }
  
  /* parse if not none */
  char *str = dup_s(par);
  char *save;
  char *main_str = sub_s(str,'(',')',&save);
  char *tok  = tok_s(main_str,',',&save);
  
  adjust[0] = dup_s(tok);
  tok = tok_s(0,',',&save);
  adjust[1] = dup_s(tok);
  tok = tok_s(0,',',&save);
  adjust[2] = dup_s(tok);
  
  free(str);
}

/* find y_CM by demanding Px_ADM = 0 */
static void Px_ADM_is0_by_y_CM(Physics_T *const phys)
{
  double dy_CM = 0,y_CM_new;
  const double px    = sysGetd("Px_ADM");
  const double py    = sysGetd("Py_ADM");
  const double W     = sysGetd("P_ADM_control_update_weight");
  const double dP    = sysGetd("P_ADM_control_tolerance");
  const double P_thr = sysGetd("P_ADM_control_threshold");
  const double Omega = sysGetd("angular_velocity");
  const double a     = infall_rate(phys);
  const double y_CM0 = sysGetd("y_CM");
  const double admM  = sysGetd("ADM_mass");
  const double P_err = (a == 0.) ? fabs(px) : sqrt(px*px+py*py);
  
  printf(Pretty0"adjusting Px_ADM by y_CM ...\n");
  
  /* changing center of mass */
  dy_CM    = (a*py-Omega*px)/(admM*(a*a+Omega*Omega));
  y_CM_new = y_CM0+dy_CM*W;
  
  /* having found new x_CM now update */
  if (GRT(P_err,dP) && LSS(P_err,P_thr))
  {
    sysSetd("y_CM",y_CM_new);
  }
}
  
/* find x_CM by demanding Py_ADM = 0 */
static void Py_ADM_is0_by_x_CM(Physics_T *const phys)
{
  double  dx_CM = 0,x_CM_new;
  const double px    = sysGetd("Px_ADM");
  const double py    = sysGetd("Py_ADM");
  const double W     = sysGetd("P_ADM_control_update_weight");
  const double dP    = sysGetd("P_ADM_control_tolerance");
  const double P_thr = sysGetd("P_ADM_control_threshold");
  const double Omega = sysGetd("angular_velocity");
  const double a     = infall_rate(phys);
  const double x_CM0 = sysGetd("x_CM");
  const double admM  = sysGetd("ADM_mass");
  const double P_err = (a == 0.) ? fabs(py) : sqrt(px*px+py*py);
  
  printf(Pretty0"adjusting Py_ADM by x_CM ...\n");
  
  /* changing center of mass */
  dx_CM    = (a*px+Omega*py)/(admM*(a*a+Omega*Omega));
  x_CM_new = x_CM0+dx_CM*W;
  
  /* having found new x_CM now update */
  if (GRT(P_err,dP) && LSS(P_err,P_thr))
  {
    sysSetd("x_CM",x_CM_new);
  }
}

static double infall_rate(Physics_T *const phys)
{
  double Vr = 0., D = 0.;
  
  sprintf(phys->par,"%s_%s",phys->ssys,"infall_velocity");
  if (PgetdEZ(phys->par) != DBL_MAX) Vr = PgetdEZ(phys->par);
  sprintf(phys->par,"%s_%s",phys->ssys,"separation");
  if (PgetdEZ(phys->par) != DBL_MAX) D = PgetdEZ(phys->par);
  
  return (D > 0.) ? Vr/D : 0.;
}

/* parse target_Jz: returns 0 if "none", else 1 and fills value */
static int parse_target_Jz(Physics_T *const phys, double *const Jt)
{
  const char *const par = sysGets("target_Jz");
  char *end = 0;

  if (!par || strcmp_i(par, "none"))
    return 0;

  *Jt = strtod(par, &end);

  if (end == par || *end != '\0' || !isfinite(*Jt))
    Errors("Syntax error for target_Jz = '%s'.\n", par);

  return 1;
}

/* Tune angular_velocity such that J_z_ADM approaches target_Jz.
// Intended for fixed-Omega runs (force balance = none).
// Call after sys_tune_ADM_momenta. */
int sys_tune_Omega_by_target_Jz(Physics_T *const phys)
{
  FUNC_TIC

  static int    have_prev = 0;
  static double Om_prev = 0., J_prev = 0.;

  double Jt = 0.;
  const char *const par = sysGets("target_Jz");
  const double W = sysGetd("Jz_control_update_weight");

  printf(Pretty0
         "Jz controller: raw target_Jz = '%s', update weight = %.8e\n",
         par ? par : "(null)", W);

  if (!parse_target_Jz(phys, &Jt))
  {
    printf(Pretty0 "Jz controller: inactive (target_Jz = none).\n");
    have_prev = 0;
  }
  else
  {
    /* force balance also tunes Omega, so the two would fight */
    if (phys->sys == BHNS &&
        !Pcmps("NS_force_balance_equation", "none"))
      Error0("target_Jz requires NS_force_balance_equation = none.\n");

    if (Jt == 0.)
      Error0("target_Jz is zero: relative error is undefined.\n");

    const double J   = sysGetd("Jz_ADM");
    const double Om0 = sysGetd("angular_velocity");
    const double tol = sysGetd("Jz_control_tolerance");
    const double thr = sysGetd("Jz_control_threshold");
    const double max_step = 0.05;
    const double err = (J-Jt)/fabs(Jt);

    printf(Pretty0
           "Jz controller: target = %.16e, current Jz = %.16e\n", Jt, J);
    printf(Pretty0
           "Jz controller: Omega = %.16e, rel. error = %.8e\n", Om0, err);
    printf(Pretty0
           "Jz controller: weight = %.8e, tolerance = %.8e, "
           "threshold = %.8e\n", W, tol, thr);

    if (Om0 == 0. || J == 0.)
      Error0("Omega or Jz_ADM is zero: cannot rescale.\n");

    if (GRT(fabs(err), tol) && LSS(fabs(err), thr))
    {
      const double s0 = J/Om0;
      double slope = s0;

      if (have_prev && fabs(Om0-Om_prev) > 1E-3*fabs(Om0))
      {
        const double s = (J-J_prev)/(Om0-Om_prev);

        /* accept secant only if within a factor of 2 of s0 */
        if (isfinite(s) && s/s0 > 0.5 && s/s0 < 2.)
          slope = s;
      }

      Om_prev   = Om0;
      J_prev    = J;
      have_prev = 1;

      double dOm = -W*(J-Jt)/slope;

      if (fabs(dOm) > max_step*fabs(Om0))
        dOm = (dOm > 0. ? 1. : -1.)*max_step*fabs(Om0);

      printf(Pretty0
             "Jz controller: updating Omega: dOmega = %.16e, "
             "new Omega = %.16e\n", dOm, Om0+dOm);

      sysSetd("angular_velocity", Om0+dOm);
    }
    else
    {
      printf(Pretty0
             "Jz controller: NO UPDATE; abs(rel. error) = %.8e, "
             "tolerance = %.8e, threshold = %.8e\n",
             fabs(err), tol, thr);
    }
  }

  FUNC_TOC
  return EXIT_SUCCESS;
}

/* parse target_ADM_mass: returns 0 if "none", else 1 and fills value */
static int parse_target_ADM_mass(Physics_T *const phys, double *const Mt)
{
  const char *const par = sysGets("target_ADM_mass");
  char *end = 0;

  if (!par || strcmp_i(par, "none"))
    return 0;

  *Mt = strtod(par, &end);

  if (end == par || *end != '\0' || !isfinite(*Mt))
    Errors("Syntax error for target_ADM_mass = '%s'.\n", par);

  return 1;
}

/* Tune infall_velocity such that ADM_mass approaches target_ADM_mass.
// M depends on v_r^2 so sensitivity is dM/dv ~ mu*v (Newtonian).
// Call after sys_tune_Omega_by_target_Jz. */
int sys_tune_infall_by_target_ADM_mass(Physics_T *const phys)
{
  FUNC_TIC

  static int    have_prev = 0;
  static double V_prev = 0., M_prev = 0.;

  double Mt = 0.;

  if (!parse_target_ADM_mass(phys, &Mt))
  {
    printf(Pretty0 "ADM mass controller: inactive (target = none).\n");
    have_prev = 0;
  }
  else
  {
    if (Mt == 0.)
      Error0("target_ADM_mass is zero: relative error is undefined.\n");

    const double M   = sysGetd("ADM_mass");
    const double V0  = sysGetd("infall_velocity");
    const double W   = sysGetd("ADM_mass_control_update_weight");
    const double tol = sysGetd("ADM_mass_control_tolerance");
    const double thr = sysGetd("ADM_mass_control_threshold");
    const double max_step = 0.05;
    const double err = (M-Mt)/fabs(Mt);

    /* Newtonian reduced mass from the masses in the parfile */
    const double m1 = Pgetd("BH_irreducible_mass");
    const double m2 = Pgetd("NS_baryonic_mass");
    const double mu = m1*m2/(m1+m2);

    printf(Pretty0
           "ADM mass controller: target = %.16e, current M = %.16e\n",
           Mt, M);
    printf(Pretty0
           "ADM mass controller: infall velocity = %.16e, "
           "rel. error = %.8e\n", V0, err);
    printf(Pretty0
           "ADM mass controller: weight = %.8e, tolerance = %.8e, "
           "threshold = %.8e\n", W, tol, thr);

    if (V0 == 0.)
    {
      printf(Pretty0
             "ADM mass controller: SKIPPED; infall_velocity = 0, "
             "so dM/dv vanishes and the target cannot be reached "
             "through it.\n");
      have_prev = 0;
    }
    else if (GRT(fabs(err), tol) && LSS(fabs(err), thr))
    {
      const double s0 = mu*V0;   /* dM/dv from M ~ M0 + mu*v^2/2 */
      double slope = s0;

      if (have_prev && fabs(V0-V_prev) > 1E-3*fabs(V0))
      {
        const double s = (M-M_prev)/(V0-V_prev);

        /* accept secant only if same sign and within a factor of 4 */
        if (isfinite(s) && s/s0 > 0.25 && s/s0 < 4.)
          slope = s;
      }

      V_prev    = V0;
      M_prev    = M;
      have_prev = 1;

      double dV = -W*(M-Mt)/slope;

      if (fabs(dV) > max_step*fabs(V0))
        dV = (dV > 0. ? 1. : -1.)*max_step*fabs(V0);

      printf(Pretty0
             "ADM mass controller: updating infall velocity: "
             "dV = %.16e, new V = %.16e\n", dV, V0+dV);

      sysSetd("infall_velocity", V0+dV);
    }
    else
    {
      printf(Pretty0
             "ADM mass controller: NO UPDATE; abs(rel. error) = %.8e, "
             "tolerance = %.8e, threshold = %.8e\n",
             fabs(err), tol, thr);
    }
  }

  FUNC_TOC
  return EXIT_SUCCESS;
}


