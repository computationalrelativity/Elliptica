/*
// Alireza Rashti
// August 2023
*/

#include "bhbh_main.h"


/* initial data for BH-BH binary system */
int BH_BH_Binary_Initial_Data(void *vp)
{
  /* if this is a generic ID reader call */
  if (strcmp_i(PgetsEZ("IDR_BHBH_export_id"),"generic"))
    bhbh_export_id_generic(vp);
  
  /* otherwise construct initial data */
  else
    construct_initial_data(vp);
  
  return EXIT_SUCCESS;
}

/* main algorithm for construction of ID */
static void construct_initial_data(void *vp)
{
  FUNC_TIC
  
  Physics_T *new_phys = 0,*old_phys = 0;
  int Stop = 0;/* if 1 it stops the main loop */
  Uint iter = 0;
  
  /* setting the default parameters */
  set_default_parameters();
  
  /* main iteration loop */
  while(!Stop)
  {
    printf("{ Outermost iteration %u ...\n",iter);
    
    Stop = update_iteration_params(iter,P_,P_"%s_%ux%ux%u");
    
    new_phys = bhbh_initialize_new_physics(old_phys);
    
    write_checkpoint(new_phys,Pgets(P_"my_directory"));
    
    bhbh_solve_equation(new_phys);
    
    bhbh_analyze(new_phys,Pgeti(P_"resolution_iteration"));
    
    free_physics(old_phys);

    fill_blackhole(new_phys);
    
    old_phys = new_phys;
    
    printf("} Outermost iteration %u ==> Done.\n",iter);
    
    iter++;
  }
  
  free_physics(new_phys);
  UNUSED(vp);
  FUNC_TOC
}

/* default parameters for this project */
static void set_default_parameters(void)
{
  /* general parameters: */
  // =================== */
  
  /* which fields to be saved at checkpoint file */
  Pset_default("checkpoint_save","enthalpy,phi,psi,alphaPsi,"
                                 "B0_U0,B0_U1,B0_U2");
  
  /* how to computer derivative */
  Pset_default("Derivative_method","Spectral");
  
  /* how to computer interpolation */
  Pset_default("Interpolation_method","Spectral");
  
  /* how to compute coefficients of expansion:
  // option:
  // =======
  // RFT : regular Fourier transformation. */
  Pset_default("Fourier_transformation_method","RFT");
  
  /* how to computer functional derivatives used in Jacobian */
  Pset_default("dF/du_for_Newton_method","Spectral");
  
  /* optimize ccs reader with 10 splits; if gives error make it smaller */
  Pset_default("matrix_ccs_reader_split","10");
  /* drop matrix entries if less than this. */
  Pset_default("matrix_ccs_drop_below","1.0E-12");
  /* use long UMFPACK */
  Pset_default("solve_UMFPACK_size","1");
  /* set the portion of cores to be used in schur algorithm.
  // the main use of this param is when the memory is scarce and
  // not all the threads should be used simultaneously otherwise 
  // it is killed by the OS. values are in [0,1]. e.g., 0.5 means
  // use 50% of the max number of available threads. */
  Pset_default("solve_ddm_schur_thread_cap","1");
  
  
  /* BH-BH parameters:
  // ================== */
 
  /* how far are BH-BH */
  Pset_default(P_"separation","50.");
  
  /* how fast BH-BH angular velocity */
  Pset_default(P_"angular_velocity","0.");
  
  /* how fast BH-BH infall velocity */
  Pset_default(P_"infall_velocity","0.");
 
  /* how to start off:
  // options:
  // "parameter_file" : it reads parameter file and initialize.
  // "checkpoint_file": it uses checkpoint file and initialize. */
  Pset_default(P_"start_off","parameter_file");
  
  /* total number of iterations that have been executed */
  Pset_default(P_"iteration","0");
  /* total number of iterations executed for this resolution */
  Pset_default(P_"resolution_iteration","0");
  /* stop the main loop if it is 1. */
  Pset_default(P_"STOP","0");
  
  /* system center of mass. */
  Pset_default(P_"x_CM","0."); 
  Pset_default(P_"y_CM","0."); 
  Pset_default(P_"z_CM","0."); 
  
  /* boost velocity for BH-BH */
  Pset_default(P_"boost_Vx","0."); 
  Pset_default(P_"boost_Vy","0."); 
  Pset_default(P_"boost_Vz","0."); 
  
  /* masses */
  Pset_default(P_"ADM_mass","1.");
  Pset_default(P_"Komar_mass","1.");
  
  /* what to print for properties of BHBH, add and separate with comma */
  Pset_default(P_"BHBH_properties",
    "separation,x_CM,y_CM,z_CM,ADM_mass,Komar_mass,mass_ratio,"
    "angular_velocity,infall_velocity,"
    "Px_ADM,Py_ADM,Pz_ADM,"
    "Jx_ADM,Jy_ADM,Jz_ADM,"
    "binding_energy,Virial_error");
    
  /* how to tune P_ADM:
  // options:[none,adjust(?_CM)], cf. */
  Pset_default(P_"P_ADM_control_method","adjust(x_CM,y_CM)");
  Pset_default(P_"P_ADM_control_update_weight","0.");
  Pset_default(P_"P_ADM_control_tolerance","1E-5");
  Pset_default(P_"P_ADM_control_threshold","1E-1");
  
  /* observer method */
  Pset_default(P_"Observe_ADM_P","S+V,constraint");
  Pset_default(P_"Observe_ADM_J","S+V,constraint");
  Pset_default(P_"Observe_ADM_M","S+V,conformal");
  Pset_default(P_"Observe_Komar_M","S_inf,default");
  
  /* equation */
  /* kind of dF/du to compute before calling elliptic solve.
  // if set "none" it won't prepare in advance. 
  // this is used for optimization.
  // NOTE: to optimize further it's better to start from second order. */
  //Pset_default(P_"dF/du_prepare","J_D0D0,J_D0D1,J_D0D2,"
  //                               "J_D1D1,J_D1D2,J_D2D2,"
  //                               "J_D0,J_D1,J_D2");
  /* NOTE: a fast and analytic method implemented thus no further 
  // dF/du preparation is requried. */
  Pset_default(P_"dF/du_prepare","none");
  

  /* BH1 paramters:
  // ============= */
  
  /* BH irreducible mass */
  Pset_default("BH1_irreducible_mass","1.");
  
  /* geometrical center of BH.
  // NOTE: geometrical center can be different from patch->c. */ 
  Pset_default("BH1_center_x","0."); 
  Pset_default("BH1_center_y","0."); 
  Pset_default("BH1_center_z","0."); 
  
  Pset_default("BH1_x_CM","0."); 
  Pset_default("BH1_y_CM","0."); 
  Pset_default("BH1_z_CM","0."); 
  
  /* box length at the center of BH */
  Pset_default("grid_BH1_central_box_length","auto");
  
  /* boost velocity for BH, 
  // if no boost required like in conformal flat metric, 
  // set them to value "off". */
  Pset_default("BH1_boost_Vx","0."); 
  Pset_default("BH1_boost_Vy","0."); 
  Pset_default("BH1_boost_Vz","0."); 
  
  /* spin vector to adjust spin for BH */
  Pset_default("BH1_Omega_x","0."); 
  Pset_default("BH1_Omega_y","0."); 
  Pset_default("BH1_Omega_z","0."); 
  
  /* tune spin */
  Pset_default("BH1_spin_update_weight","0."); 
  Pset_default("BH1_spin_tolerance","1E-3"); 
  
  /* spin */
  Pset_default("BH1_chi_x","0."); 
  Pset_default("BH1_chi_y","0."); 
  Pset_default("BH1_chi_z","0."); 
  Pset_default("BH1_spin_x","0."); 
  Pset_default("BH1_spin_y","0."); 
  Pset_default("BH1_spin_z","0."); 
  
  /* ADM momentum */
  Pset_default("BH1_Px_ADM","0."); 
  Pset_default("BH1_Py_ADM","0."); 
  Pset_default("BH1_Pz_ADM","0."); 
  
  /* ADM angular momentum */
  Pset_default("BH1_Jx_ADM","0."); 
  Pset_default("BH1_Jy_ADM","0."); 
  Pset_default("BH1_Jz_ADM","0."); 
  
  /* what to print for properties of BH, add and separate with comma */
  Pset_default(P_"BH1_properties",
   "center_x,center_y,center_z,x_CM,y_CM,z_CM,max_radius,min_radius,"
   "Komar_mass,irreducible_mass,irreducible_mass_current,"
   "Christodoulou_mass,Christodoulou_mass_current,"
   "Omega_x,Omega_y,Omega_z,chi_x,chi_y,chi_z,"
   "chi_x_current,chi_y_current,chi_z_current,"
   "spin_x,spin_y,spin_z,"
   "Px_ADM,Py_ADM,Pz_ADM,"
   "Jx_ADM,Jy_ADM,Jz_ADM,"
   "boost_Vx,boost_Vy,boost_Vz");
 
  /* the very first BH approximation
  // options:
  // ========
  // o. CloseKerrSchild (see BH physics)
  // o. IsoSchild (see BH physics) */
  Pset_default("BH1_start_off","CloseKerrSchild"); 
  
  /* roll off function to stich free data.
  // options:
  // ========
  // o. exp(-lambda*(r/rmax)^p):r<rmax.
  // o. exp(-lambda*(r/rmax)^p). */
  Pset_default("BH1_RollOff_function","exp(-lambda*(r/rmax)^p)");
  
  /* lambda in "BH1_RollOff_function".
  // options:
  // ========
  // o. constant_1.
  // o. |(r-rmin)/(rmax-r)|. */  
  Pset_default("BH1_RollOff_lambda","constant_1");
  
  /* rmax in "BH1_RollOff_function" */
  Pset_default("BH1_RollOff_rmax","auto");
  
  /* p in "BH1_RollOff_function" */
  Pset_default("BH1_RollOff_power","4.");
  
  /* max l in Ylm expansion */
  Pset_default("BH1_surface_Ylm_max_l","1"); 
  
  /* BH surface type */
  Pset_default("BH1_surface_type","perfect_s2"); 
  
  /* did_BH_surface_change? used for interpolation purposes and efficiency */
  Pset_default("BH1_did_BH_surface_change?","1");
  
  /* was BH filled? used to tack whether BH is filled or not */
  Pset_default("BH1_was_BH_filled?","0");
  
  /* equation related: */
  /* XCTS means: alpha,beta,psi */
  Pset_default("BH1_Eq_inner_BC_fields","XCTS");
  
  /* set alpha on AH to be W*KerrSchild value. */
  Pset_default("BH1_Eq_inner_BC_alpha","w*KerrSchild");
  
  /* set bete^i on AH to be alpha*s^i + Omega x r.  */
  Pset_default("BH1_Eq_inner_BC_beta","alpha+Omega*r");
  
  /* observe method pertinet to BH */
  Pset_default("BH1_Observe_ADM_M","S_obj,default");
  Pset_default("BH1_Observe_Komar_M","S_obj,default");
  Pset_default("BH1_Observe_Irreducible_M","S_obj,default");
  Pset_default("BH1_Observe_ADM_P","S_obj,default");
  Pset_default("BH1_Observe_ADM_J","S_obj,default");
  Pset_default("BH1_Observe_CM","S_obj,default");
  Pset_default("BH1_Observe_spin","S_obj,Campanelli");
  
  /* how to tune BH Radius */
  Pset_default("BH1_tune_BH_radius_criteria","fix_irreducible_mass");
  Pset_default("BH1_mass_tolerance","1E-5");
  Pset_default("BH1_radius_update_weight","0.");
  
  /* BH filler parameters:
  // { */
  /* how to fill: 
  // options:
  // ChebTn_Ylm_perfect_s2:
  //    fill PERFECT S2 surface excised BH with data 
  //    demanding C2 continuity across horizon. extrapolant is:
  //    f(r,th,ph) = C_{ilm}*ChebyshevT(i,r)*Ylm(th,ph). 
  //    thie method is faster than ChebTn_general_s2
  // ChebTn_general_s2: 
  //    fill a GENERAL S2 surface excised BH with data 
  //    demanding C2 continuity across horizon. extrapolant is:
  //    f(r) = C_{i}ChebyshevT(i,r) along radius. 
  // None:
  //    No filling. 
  // expmr_C0_perfect_s2:
  //    using f(r) = f(r0)e^-r (this is preferred during solve).
  // r_expmr_C1_perfect_s2: 
  //    using f(r) = (a+b*r)e^-r.
  */
  Pset_default("BH1_filler_method","expmr_C0_perfect_s2");
  /* max l for Y_{lm} expansion in ChebTn_Ylm_perfect_s2 */
  Pset_default("BH1_filler_Ylm_expansion_lmax","15");
  /* it won't fill inside the BH if r < Rmin_cutoff.
  // note: this only kicks in for ChebTn_general_s2 method.
  // the ChebTn_Ylm_perfect_s2 is so fast which we don't bother. */
  Pset_default("BH1_filler_Rmin_cutoff","-1.");
  /* which fields to be filled */
  Pset_default("BH1_filler_fields","alphaPsi,psi,beta_U0,beta_U1,beta_U2,"
                                  "adm_Kij_D0D0,adm_Kij_D0D1,adm_Kij_D0D2,"
                                  "adm_Kij_D1D1,adm_Kij_D1D2,adm_Kij_D2D2,"
                                  "gConf_D0D0,gConf_D0D1,gConf_D0D2,"
                                  "gConf_D1D1,gConf_D1D2,gConf_D2D2");
  /* value of extrapolated fields at the center of BH r = 0 */
  Pset_default("BH1_filler_r0_trK","0.");
  Pset_default("BH1_filler_r0_alpha","0.2");
  Pset_default("BH1_filler_r0_psi","2.");
  Pset_default("BH1_filler_r0_alphaPsi","0.4");
  Pset_default("BH1_filler_r0_beta_U0","0.");
  Pset_default("BH1_filler_r0_beta_U1","0.");
  Pset_default("BH1_filler_r0_beta_U2","0.");
  Pset_default("BH1_filler_r0_B0_U0","0.");
  Pset_default("BH1_filler_r0_B0_U1","0.");
  Pset_default("BH1_filler_r0_B0_U2","0.");
  Pset_default("BH1_filler_r0_gConf_D0D0","4.");
  Pset_default("BH1_filler_r0_gConf_D0D1","0.");
  Pset_default("BH1_filler_r0_gConf_D0D2","0.");
  Pset_default("BH1_filler_r0_gConf_D1D1","4.");
  Pset_default("BH1_filler_r0_gConf_D1D2","0.");
  Pset_default("BH1_filler_r0_gConf_D2D2","4.");
  Pset_default("BH1_filler_r0_adm_Kij_D0D0","0.");
  Pset_default("BH1_filler_r0_adm_Kij_D0D1","0.");
  Pset_default("BH1_filler_r0_adm_Kij_D0D2","0.");
  Pset_default("BH1_filler_r0_adm_Kij_D1D1","0.");
  Pset_default("BH1_filler_r0_adm_Kij_D1D2","0.");
  Pset_default("BH1_filler_r0_adm_Kij_D2D2","0.");
  /* verbos? [1 = yes; 0 = no] */
  Pset_default("BH1_filler_verbose","0");
  
  /* test BH filler: interpolate along this line:
  // x = init_x + t*slop_x
  // y = init_y + t*slop_y
  // z = init_z + t*slop_z. */
  Pset_default("BH1_filler_test_print_1d","no");
  Pset_default("BH1_filler_test_print_1d_slop_x","0.");
  Pset_default("BH1_filler_test_print_1d_slop_y","0.");
  Pset_default("BH1_filler_test_print_1d_slop_z","0.");
  Pset_default("BH1_filler_test_print_1d_init_x","0.");
  Pset_default("BH1_filler_test_print_1d_init_y","0.");
  Pset_default("BH1_filler_test_print_1d_init_z","0.");
  Pset_default("BH1_filler_test_print_1d_length","20.");
  Pset_default("BH1_filler_test_print_1d_points","1");
  /* } */
  
  /* modify parameters if needed */
  
  /* auto roll off param. 
  // note: too small roll off might not get you a converge solution
  // the bigger the better however, the metric should become flat at star
  // and around it, so not too big! */
  if (Pcmps("BH1_RollOff_rmax","auto"))
  {
    double rmax = Pgetd(P_"separation")/2.;
    assert(rmax > 0);
    Psetd("BH1_RollOff_rmax",rmax);
  }

  /* BH2 paramters:
  // ============= */
  
  /* BH irreducible mass */
  Pset_default("BH2_irreducible_mass","1.");
  
  /* geometrical center of BH.
  // NOTE: geometrical center can be different from patch->c. */ 
  Pset_default("BH2_center_x","0."); 
  Pset_default("BH2_center_y","0."); 
  Pset_default("BH2_center_z","0."); 
  
  Pset_default("BH2_x_CM","0."); 
  Pset_default("BH2_y_CM","0."); 
  Pset_default("BH2_z_CM","0."); 
  
  /* box length at the center of BH */
  Pset_default("grid_BH2_central_box_length","auto");
  
  /* boost velocity for BH, 
  // if no boost required like in conformal flat metric, 
  // set them to value "off". */
  Pset_default("BH2_boost_Vx","0."); 
  Pset_default("BH2_boost_Vy","0."); 
  Pset_default("BH2_boost_Vz","0."); 
  
  /* spin vector to adjust spin for BH */
  Pset_default("BH2_Omega_x","0."); 
  Pset_default("BH2_Omega_y","0."); 
  Pset_default("BH2_Omega_z","0."); 
  
  /* tune spin */
  Pset_default("BH2_spin_update_weight","0."); 
  Pset_default("BH2_spin_tolerance","1E-3"); 
  
  /* spin */
  Pset_default("BH2_chi_x","0."); 
  Pset_default("BH2_chi_y","0."); 
  Pset_default("BH2_chi_z","0."); 
  Pset_default("BH2_spin_x","0."); 
  Pset_default("BH2_spin_y","0."); 
  Pset_default("BH2_spin_z","0."); 
  
  /* ADM momentum */
  Pset_default("BH2_Px_ADM","0."); 
  Pset_default("BH2_Py_ADM","0."); 
  Pset_default("BH2_Pz_ADM","0."); 
  
  /* ADM angular momentum */
  Pset_default("BH2_Jx_ADM","0."); 
  Pset_default("BH2_Jy_ADM","0."); 
  Pset_default("BH2_Jz_ADM","0."); 
  
  /* what to print for properties of BH, add and separate with comma */
  Pset_default(P_"BH2_properties",
   "center_x,center_y,center_z,x_CM,y_CM,z_CM,max_radius,min_radius,"
   "Komar_mass,irreducible_mass,irreducible_mass_current,"
   "Christodoulou_mass,Christodoulou_mass_current,"
   "Omega_x,Omega_y,Omega_z,chi_x,chi_y,chi_z,"
   "chi_x_current,chi_y_current,chi_z_current,"
   "spin_x,spin_y,spin_z,"
   "Px_ADM,Py_ADM,Pz_ADM,"
   "Jx_ADM,Jy_ADM,Jz_ADM,"
   "boost_Vx,boost_Vy,boost_Vz");
 
  /* the very first BH approximation
  // options:
  // ========
  // o. CloseKerrSchild (see BH physics)
  // o. IsoSchild (see BH physics) */
  Pset_default("BH2_start_off","CloseKerrSchild"); 
  
  /* roll off function to stich free data.
  // options:
  // ========
  // o. exp(-lambda*(r/rmax)^p):r<rmax.
  // o. exp(-lambda*(r/rmax)^p). */
  Pset_default("BH2_RollOff_function","exp(-lambda*(r/rmax)^p)");
  
  /* lambda in "BH2_RollOff_function".
  // options:
  // ========
  // o. constant_1.
  // o. |(r-rmin)/(rmax-r)|. */  
  Pset_default("BH2_RollOff_lambda","constant_1");
  
  /* rmax in "BH2_RollOff_function" */
  Pset_default("BH2_RollOff_rmax","auto");
  
  /* p in "BH2_RollOff_function" */
  Pset_default("BH2_RollOff_power","4.");
  
  /* max l in Ylm expansion */
  Pset_default("BH2_surface_Ylm_max_l","1"); 
  
  /* BH surface type */
  Pset_default("BH2_surface_type","perfect_s2"); 
  
  /* did_BH_surface_change? used for interpolation purposes and efficiency */
  Pset_default("BH2_did_BH_surface_change?","1");
  
  /* was BH filled? used to tack whether BH is filled or not */
  Pset_default("BH2_was_BH_filled?","0");
  
  /* equation related: */
  /* XCTS means: alpha,beta,psi */
  Pset_default("BH2_Eq_inner_BC_fields","XCTS");
  
  /* set alpha on AH to be W*KerrSchild value. */
  Pset_default("BH2_Eq_inner_BC_alpha","w*KerrSchild");
  
  /* set bete^i on AH to be alpha*s^i + Omega x r.  */
  Pset_default("BH2_Eq_inner_BC_beta","alpha+Omega*r");
  
  /* observe method pertinet to BH */
  Pset_default("BH2_Observe_ADM_M","S_obj,default");
  Pset_default("BH2_Observe_Komar_M","S_obj,default");
  Pset_default("BH2_Observe_Irreducible_M","S_obj,default");
  Pset_default("BH2_Observe_ADM_P","S_obj,default");
  Pset_default("BH2_Observe_ADM_J","S_obj,default");
  Pset_default("BH2_Observe_CM","S_obj,default");
  Pset_default("BH2_Observe_spin","S_obj,Campanelli");
  
  /* how to tune BH Radius */
  Pset_default("BH2_tune_BH_radius_criteria","fix_irreducible_mass");
  Pset_default("BH2_mass_tolerance","1E-5");
  Pset_default("BH2_radius_update_weight","0.");
  
  /* BH filler parameters:
  // { */
  /* how to fill: 
  // options:
  // ChebTn_Ylm_perfect_s2:
  //    fill PERFECT S2 surface excised BH with data 
  //    demanding C2 continuity across horizon. extrapolant is:
  //    f(r,th,ph) = C_{ilm}*ChebyshevT(i,r)*Ylm(th,ph). 
  //    thie method is faster than ChebTn_general_s2
  // ChebTn_general_s2: 
  //    fill a GENERAL S2 surface excised BH with data 
  //    demanding C2 continuity across horizon. extrapolant is:
  //    f(r) = C_{i}ChebyshevT(i,r) along radius. 
  // None:
  //    No filling. 
  // expmr_C0_perfect_s2:
  //    using f(r) = f(r0)e^-r (this is preferred during solve).
  // r_expmr_C1_perfect_s2: 
  //    using f(r) = (a+b*r)e^-r.
  */
  Pset_default("BH2_filler_method","expmr_C0_perfect_s2");
  /* max l for Y_{lm} expansion in ChebTn_Ylm_perfect_s2 */
  Pset_default("BH2_filler_Ylm_expansion_lmax","15");
  /* it won't fill inside the BH if r < Rmin_cutoff.
  // note: this only kicks in for ChebTn_general_s2 method.
  // the ChebTn_Ylm_perfect_s2 is so fast which we don't bother. */
  Pset_default("BH2_filler_Rmin_cutoff","-1.");
  /* which fields to be filled */
  Pset_default("BH2_filler_fields","alphaPsi,psi,beta_U0,beta_U1,beta_U2,"
                                  "adm_Kij_D0D0,adm_Kij_D0D1,adm_Kij_D0D2,"
                                  "adm_Kij_D1D1,adm_Kij_D1D2,adm_Kij_D2D2,"
                                  "gConf_D0D0,gConf_D0D1,gConf_D0D2,"
                                  "gConf_D1D1,gConf_D1D2,gConf_D2D2");
  /* value of extrapolated fields at the center of BH r = 0 */
  Pset_default("BH2_filler_r0_trK","0.");
  Pset_default("BH2_filler_r0_alpha","0.2");
  Pset_default("BH2_filler_r0_psi","2.");
  Pset_default("BH2_filler_r0_alphaPsi","0.4");
  Pset_default("BH2_filler_r0_beta_U0","0.");
  Pset_default("BH2_filler_r0_beta_U1","0.");
  Pset_default("BH2_filler_r0_beta_U2","0.");
  Pset_default("BH2_filler_r0_B0_U0","0.");
  Pset_default("BH2_filler_r0_B0_U1","0.");
  Pset_default("BH2_filler_r0_B0_U2","0.");
  Pset_default("BH2_filler_r0_gConf_D0D0","4.");
  Pset_default("BH2_filler_r0_gConf_D0D1","0.");
  Pset_default("BH2_filler_r0_gConf_D0D2","0.");
  Pset_default("BH2_filler_r0_gConf_D1D1","4.");
  Pset_default("BH2_filler_r0_gConf_D1D2","0.");
  Pset_default("BH2_filler_r0_gConf_D2D2","4.");
  Pset_default("BH2_filler_r0_adm_Kij_D0D0","0.");
  Pset_default("BH2_filler_r0_adm_Kij_D0D1","0.");
  Pset_default("BH2_filler_r0_adm_Kij_D0D2","0.");
  Pset_default("BH2_filler_r0_adm_Kij_D1D1","0.");
  Pset_default("BH2_filler_r0_adm_Kij_D1D2","0.");
  Pset_default("BH2_filler_r0_adm_Kij_D2D2","0.");
  /* verbos? [1 = yes; 0 = no] */
  Pset_default("BH2_filler_verbose","0");
  
  /* test BH filler: interpolate along this line:
  // x = init_x + t*slop_x
  // y = init_y + t*slop_y
  // z = init_z + t*slop_z. */
  Pset_default("BH2_filler_test_print_1d","no");
  Pset_default("BH2_filler_test_print_1d_slop_x","0.");
  Pset_default("BH2_filler_test_print_1d_slop_y","0.");
  Pset_default("BH2_filler_test_print_1d_slop_z","0.");
  Pset_default("BH2_filler_test_print_1d_init_x","0.");
  Pset_default("BH2_filler_test_print_1d_init_y","0.");
  Pset_default("BH2_filler_test_print_1d_init_z","0.");
  Pset_default("BH2_filler_test_print_1d_length","20.");
  Pset_default("BH2_filler_test_print_1d_points","1");
  /* } */
  
  /* modify parameters if needed */
  
  /* auto roll off param. 
  // note: too small roll off might not get you a converge solution
  // the bigger the better however, the metric should become flat at star
  // and around it, so not too big! */
  if (Pcmps("BH2_RollOff_rmax","auto"))
  {
    double rmax = Pgetd(P_"separation")/2.;
    assert(rmax > 0);
    Psetd("BH2_RollOff_rmax",rmax);
  }

  /* use Kepler's law to set angular velocity */
  if (Pcmps(P_"angular_velocity","auto"))
  {
    const double m1 = Pgetd("BH1_irreducible_mass");
    const double m2 = Pgetd("BH2_irreducible_mass");
    
    const double r  = Pgetd(P_"separation");
    const double O  = sqrt((m1+m2)/pow(r,3.));
    
    Psetd(P_"angular_velocity",O);
  }
  
  /* set BH1 and BH2 centers */
  if (strstr(Pgets("grid_set_BH1"),"right"))
  {
    const double S = Pgetd(P_"separation");
    
    /* BH1 center in +y */
    Psetd("BH1_center_x",0.);
    Psetd("BH1_center_y",S/2.);
    Psetd("BH1_center_z",0.);
    
    /* BH2 center in -y */
    Psetd("BH2_center_x",0.);
    Psetd("BH2_center_y",-S/2.);
    Psetd("BH2_center_z",0.);
  }
  else
  {
    const double S = Pgetd(P_"separation");

    /* BH1 center in -y */
    Psetd("BH1_center_x",0.);
    Psetd("BH1_center_y",-S/2.);
    Psetd("BH1_center_z",0.);
    
    /* BH2 center in +y */
    Psetd("BH2_center_x",0.);
    Psetd("BH2_center_y",S/2.);
    Psetd("BH2_center_z",0.);
    
  }
  
  /* set center of mass */
  {
    const double irr_mass1 = Pgetd("BH1_irreducible_mass");
    const double irr_mass2 = Pgetd("BH2_irreducible_mass");
    const double y_CM = (irr_mass1*Pgetd("BH1_center_y") +
                         irr_mass2*Pgetd("BH2_center_y"))/(irr_mass1+irr_mass2);
   
    Psetd(P_"x_CM",0.);
    Psetd(P_"y_CM",y_CM);
    Psetd(P_"z_CM",0.);
  }
  
  /* boost speed for BH1 in x direction */
  if (!Pcmps("BH1_boost_Vx","off"))
  {
    const double Omega = Pgetd(P_"angular_velocity");
    const double y_CM  = Pgetd(P_"y_CM");
    const double BH_center_y = Pgetd("BH1_center_y");
    
    Psetd("BH1_boost_Vx",-Omega*(BH_center_y-y_CM));
  }
  
  /* boost speed for BH2 in x direction */
  if (!Pcmps("BH2_boost_Vx","off"))
  {
    const double Omega = Pgetd(P_"angular_velocity");
    const double y_CM  = Pgetd(P_"y_CM");
    const double BH_center_y = Pgetd("BH2_center_y");
    
    Psetd("BH2_boost_Vx",-Omega*(BH_center_y-y_CM));
  }
  
}

/* if nessesary fill the bh for interpolation purposes */
static void fill_blackhole(Physics_T *const phys)
{
  if (!phys) return;
 
  FUNC_TIC
  
  Physics_T *const bh1 = init_physics(phys,BH1);
  
  /* BH1 is not filled yet so: */
  Pseti("BH1_was_BH_filled?",0);
  
  if (Pgeti("BH1_did_BH_surface_change?"))
  {
    /* fill the hole if we have excised region */
    if (Pcmpss("grid_set_BH1","excised"))
    {
      Psets("BH1_filler_fields","alphaPsi,psi,B0_U0,B0_U1,B0_U2");
      physics(bh1,BH_FILL);
      Pseti("BH1_was_BH_filled?",1);
    }
  }
  free_physics(bh1);
  
  Physics_T *const bh2 = init_physics(phys,BH2);
  
  /* BH2 is not filled yet so: */
  Pseti("BH2_was_BH_filled?",0);
  
  if (Pgeti("BH2_did_BH_surface_change?"))
  {
    /* fill the hole if we have excised region */
    if (Pcmpss("grid_set_BH2","excised"))
    {
      Psets("BH2_filler_fields","alphaPsi,psi,B0_U0,B0_U1,B0_U2");
      physics(bh2,BH_FILL);
      Pseti("BH2_was_BH_filled?",1);
    }
  }
  free_physics(bh2);
  
  FUNC_TOC
}
