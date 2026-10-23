/*
// Alireza Rashti
// August 2023
*/

/* analyzing initial data such as mass, momentum, constraints etc.  */

#include "bhbh_analyze.h"

/* analyzing physics properties, constraints etc */
void bhbh_analyze(Physics_T *const phys,const int iteration)
{
  if (!phys) return;

  FUNC_TIC
  
  const char *const properties_file_name = P_"properties.txt";
  FILE *file = 0;
  char str[MAX_STR_LEN];
   
  /* compute properties and constraints */ 
  physics(phys,ADM_COMPUTE_CONSTRAINTS);
  
  /* compute various properties */
  compute_properties(phys);
  
  /* open properties file in "my_directory" and save */
  sprintf(str,"%s/%s",Pgets(P_"my_directory"),properties_file_name);
  file = Fopen(str,"w");
  bhbh_print_physical_system_properties(phys,file,iteration,0);
  Fclose(file);

  /* open properties file in "Diagnostics" and save */
  sprintf(str,"%s/%s",Pgets(P_"Diagnostics"),properties_file_name);
  file = Fopen(str,"a");
  bhbh_print_physical_system_properties(phys,file,iteration,0);
  Fclose(file);
  
  /* prints */
  print_fields_0D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_1D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_2D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  print_fields_3D(phys->grid,iteration,Pgets(P_"Diagnostics"));
  
  FUNC_TOC
}

/* print physical system properties such as mass, spin etc in the given
// file, if pr_screen is 1, it also prints in stdout */
void bhbh_print_physical_system_properties(Physics_T *const phys,
                                          FILE *const file,
                                          const int iteration,
                                          const int pr_screen)
{
  Physics_T *const bh1 = init_physics(phys,BH1);
  Physics_T *const bh2 = init_physics(phys,BH2);

  if (pr_screen)
  {
    printf(Pretty0"iteration = %d:\n",iteration);
  }
  fprintf(file,"%s\n",LINE_STR);
  fprintf(file,"# iteration = %d\n",iteration);
  fprintf(file,"\n");
  
  bh_print_properties(bh1,Pgets(P_"BH1_properties"),file,pr_screen);
  bh_print_properties(bh2,Pgets(P_"BH2_properties"),file,pr_screen);
  sys_print_properties(phys,Pgets(P_"BHBH_properties"),file,pr_screen);
  
  free_physics(bh1);
  free_physics(bh2);
}

/* compute variety of properties.
// NOTE: order of parameter calculations matter. 
// NOTE: if there is a confusion between target params and current
//       params, "current" suffix added to the latter. */
static void compute_properties(Physics_T *const phys/* bhbh */)
{
  Physics_T *const bh1 = init_physics(phys,BH1);
  Physics_T *const bh2 = init_physics(phys,BH2);
  
  const double x_CM = Pgetd(P_"x_CM");
  const double y_CM = Pgetd(P_"y_CM");
  const double z_CM = Pgetd(P_"z_CM");
  double im[2] = {0.};
  double p[3]  = {0.};
  double j[3]  = {0.};
  double s[3]  = {0.};
  double cm[3] = {0.};
  double m     = 0.;
  double irr_mass = 0.;
  double net_spin = 0.;
  
  /* BH1: */
  observe(bh1,"Komar(M)",Pgets("BH1_Observe_Komar_M"),&m);
  Psetd("BH1_Komar_mass",m);
  
  observe(bh1,"Irreducible(M)",Pgets("BH1_Observe_irreducible_M"),im);
  Psetd("BH1_irreducible_mass_current",im[0]);
  Psetd("BH1_AH_area",im[1]);
  
  observe(bh1,"CM",Pgets("BH1_Observe_CM"),cm);
  Psetd("BH1_x_CM",cm[0]+x_CM);
  Psetd("BH1_y_CM",cm[1]+y_CM);
  Psetd("BH1_z_CM",cm[2]+z_CM);

  observe(bh1,"ADM(P)",Pgets("BH1_Observe_ADM_P"),p);
  Psetd("BH1_Px_ADM",p[0]);
  Psetd("BH1_Py_ADM",p[1]);
  Psetd("BH1_Pz_ADM",p[2]);

  observe(bh1,"ADM(J)",Pgets("BH1_Observe_ADM_J"),j);
  Psetd("BH1_Jx_ADM",j[0]);
  Psetd("BH1_Jy_ADM",j[1]);
  Psetd("BH1_Jz_ADM",j[2]);
  
  observe(bh1,"spin",Pgets("BH1_Observe_spin"),s);
  Psetd("BH1_Spin_x",s[0]);
  Psetd("BH1_Spin_y",s[1]);
  Psetd("BH1_Spin_z",s[2]);
  
  /* calculate BH current Christodoulou mass.
  // NOTE: s[?] depend on BH spin on top */
  irr_mass = Pgetd("BH1_irreducible_mass_current");
  net_spin = sqrt(Pow2(s[0])+Pow2(s[1])+Pow2(s[2]));
  m = sqrt(Pow2(irr_mass)+Pow2(net_spin)/(4*Pow2(irr_mass)));
  Psetd("BH1_Christodoulou_mass_current",m);
  
  Psetd("BH1_chi_x_current",s[0]/Pow2(m));
  Psetd("BH1_chi_y_current",s[1]/Pow2(m));
  Psetd("BH1_chi_z_current",s[2]/Pow2(m));

  /* BH2: */
  observe(bh2,"Komar(M)",Pgets("BH2_Observe_Komar_M"),&m);
  Psetd("BH2_Komar_mass",m);
  
  observe(bh2,"Irreducible(M)",Pgets("BH2_Observe_irreducible_M"),im);
  Psetd("BH2_irreducible_mass_current",im[0]);
  Psetd("BH2_AH_area",im[1]);
  
  observe(bh2,"CM",Pgets("BH2_Observe_CM"),cm);
  Psetd("BH2_x_CM",cm[0]+x_CM);
  Psetd("BH2_y_CM",cm[1]+y_CM);
  Psetd("BH2_z_CM",cm[2]+z_CM);

  observe(bh2,"ADM(P)",Pgets("BH2_Observe_ADM_P"),p);
  Psetd("BH2_Px_ADM",p[0]);
  Psetd("BH2_Py_ADM",p[1]);
  Psetd("BH2_Pz_ADM",p[2]);

  observe(bh2,"ADM(J)",Pgets("BH2_Observe_ADM_J"),j);
  Psetd("BH2_Jx_ADM",j[0]);
  Psetd("BH2_Jy_ADM",j[1]);
  Psetd("BH2_Jz_ADM",j[2]);
  
  observe(bh2,"spin",Pgets("BH2_Observe_spin"),s);
  Psetd("BH2_Spin_x",s[0]);
  Psetd("BH2_Spin_y",s[1]);
  Psetd("BH2_Spin_z",s[2]);
  
  /* calculate BH current Christodoulou mass.
  // NOTE: s[?] depend on BH spin on top */
  irr_mass = Pgetd("BH2_irreducible_mass_current");
  net_spin = sqrt(Pow2(s[0])+Pow2(s[1])+Pow2(s[2]));
  m = sqrt(Pow2(irr_mass)+Pow2(net_spin)/(4*Pow2(irr_mass)));
  Psetd("BH2_Christodoulou_mass_current",m);
  
  Psetd("BH2_chi_x_current",s[0]/Pow2(m));
  Psetd("BH2_chi_y_current",s[1]/Pow2(m));
  Psetd("BH2_chi_z_current",s[2]/Pow2(m));
  
  /* BHBH: */
  observe(phys,"ADM(M)",Pgets(P_"Observe_ADM_M"),&m);
  Psetd(P_"adm_mass",m);
  
  observe(phys,"Komar(M)",Pgets(P_"Observe_Komar_M"),&m);
  Psetd(P_"Komar_mass",m);
  
  observe(phys,"ADM(P)",Pgets(P_"Observe_ADM_P"),p);
  Psetd(P_"Px_ADM",p[0]);
  Psetd(P_"Py_ADM",p[1]);
  Psetd(P_"Pz_ADM",p[2]);

  observe(phys,"ADM(J)",Pgets(P_"Observe_ADM_J"),j);
  Psetd(P_"Jx_ADM",j[0]);
  Psetd(P_"Jy_ADM",j[1]);
  Psetd(P_"Jz_ADM",j[2]);

  double q = 
   Pgetd("BH1_Christodoulou_mass_current")/Pgetd("BH2_Christodoulou_mass_current");
  Psetd(P_"mass_ratio",q);

  /* binding energy */
  double bin_e = Pgetd(P_"adm_mass")-
    (Pgetd("BH1_Christodoulou_mass_current")+Pgetd("BH2_Christodoulou_mass_current"));
  Psetd(P_"binding_energy",bin_e);

  /* virial error */
  double v_e = fabs(1.-Pgetd(P_"adm_mass")/Pgetd(P_"komar_mass"));
  Psetd(P_"virial_error",v_e);
 
  if (0)/* too much inaccurate */
  {
    double m1    = Pgetd("BH1_adm_mass");
    double m2    = Pgetd("BH2_adm_mass");
    double nu    = m1*m2/Pow2(m1+m2);
    double omega = Pgetd(P_"angular_velocity");
    double m_tot = Pgetd(P_"ADM_mass");
    double N_orb = pow(m_tot*omega,-5./3.)/(32.*nu)/(2.*M_PI);
    /* initially some vars might be off or negative like total adm_mass */
    N_orb = (isfinite(N_orb) && N_orb > 0. ? N_orb : 0.);
    Psetd(P_"number_of_orbits_1PN",N_orb);
  }
  
  free_physics(bh1);
  free_physics(bh2);
}
