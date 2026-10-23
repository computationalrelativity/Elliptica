/*
// Alireza Rashti
// June 2019
*/

#include "eos_test.h"

/* write a row of all thermodynamic variables at eos->h, then end */
static void write_tabular_row(EoS_T *const eos,FILE *const file,
                              const char *const end)
{
  const double h = eos->h;
  fprintf(file,"%+0.15e ",h);
  eos->h = h; fprintf(file,"%+0.15e ",eos->pressure(eos));
  eos->h = h; fprintf(file,"%+0.15e ",eos->rest_mass_density(eos));
  eos->h = h; fprintf(file,"%+0.15e ",eos->energy_density(eos));
  eos->h = h; fprintf(file,"%+0.15e ",eos->specific_internal_energy(eos));
  eos->h = h; fprintf(file,"%+0.15e ",eos->drho0_dh(eos));
  eos->h = h; fprintf(file,"%+0.15e ",eos->de_dh(eos));
  eos->h = h; fprintf(file,"%+0.15e ",eos->export_rest_mass_density(eos));
  eos->h = h; fprintf(file,"%+0.15e",eos->export_specific_internal_energy(eos));
  fprintf(file,"%s",end);
  eos->h = h;
}

void test_EoS(Physics_T *const phys)
{
  Grid_T *const grid = phys->grid;
  EoS_T *eos = init_EoS(phys);
  const char *const path_par = Pgets("top_directory");
  char *path,file_name[1000];
  FILE *file = 0;
  /* number of samples, optional parameter */
  const int N_par = PgetiEZ("Test_EoS_num_points");
  Uint N = (N_par != INT_MAX && N_par > 2) ? (Uint)N_par : 1000;
  const double h_max = PgetdEZ(Ftype("NS_EoS_enthalpy_ceiling")) == DBL_MAX ?
                       3.00 : Getd(P_"enthalpy_ceiling");
  const double h_min = PgetdEZ(Ftype("NS_EoS_enthalpy_floor")) == DBL_MAX ?
                       1.00 : Getd(P_"enthalpy_floor");
  double s = (h_max-h_min)/(N-1);
  Uint i;
  path = make_directory(path_par,"EoS_Tests");
  
  /* values */
  /* if it is pwp */
  if (strstr_i(eos->type,"piecewise_polytropic") ||
      strstr_i(eos->type,"pwp"))
  {
    sprintf(file_name,"%s/%s_pwp.txt",path,eos->description);
    file = Fopen(file_name,"w+");
    
    fprintf(file,"[1]:piece  [2]:Kappa  [3]:rho  [4]:gamma  [5]:a  [6]:h-1\n");
    for (i = 0; i < eos->N; ++i)  
    {
      fprintf(file,"%u %0.15f %0.15f %0.15f %0.15f %0.15f\n",
        i,eos->K[i],eos->rho0_th[i],eos->gamma[i],eos->a[i],eos->h_th[i]-1);
    }
    Fclose(file);
  }
    
  /* continuity */
  sprintf(file_name,"%s/%s.txt",path,"pressure");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  pressure\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->pressure(eos));
  }
  Fclose(file);

  sprintf(file_name,"%s/%s.txt",path,"rest_mass_density");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  rest_mass_density\n");  
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->rest_mass_density(eos));
  }
  Fclose(file);
  
  sprintf(file_name,"%s/%s.txt",path,"energy_density");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  energy_density\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->energy_density(eos));
  }
  Fclose(file);
  
  sprintf(file_name,"%s/%s.txt",path,"specific_internal_energy");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  specific_internal_energy\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->specific_internal_energy(eos));
  }
  Fclose(file);
  
  sprintf(file_name,"%s/%s.txt",path,"drho0_dh");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  drho0_dh\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->drho0_dh(eos));
  }
  Fclose(file);
      
  sprintf(file_name,"%s/%s.txt",path,"de_dh");
  file = Fopen(file_name,"w+");
  fprintf(file,"# enthalpy  de_dh\n");
  for (i = 0; i < N; ++i)
  {
    eos->h = h_min+s*i;
    fprintf(file,"%0.15e %0.15e\n",eos->h,eos->de_dh(eos));
  }
  Fclose(file);
 
  if (strcmp_i(eos->type, "tabular") || 
      strcmp_i(eos->type, "tab")     || 
      strcmp_i(eos->type, "table"))
  {
    sprintf(file_name,"%s/%s_tabular.txt",path,eos->description);
    file = Fopen(file_name,"w");
    // header
    fprintf(file,"# [1]:enthalpy(h) [2]:pressure(p) [3]:rest_mass_density(rho0) "
                 "[4]:energy_density(e) [5]:specific_internal_energy(e0) "
                 "[6]:drho0/dh [7]:de/dh [8]:export_rest_mass_density "
                 "[9]:export_specific_internal_energy\n");
                 
    for (i = 0; i < N; ++i)
    {
      eos->h = h_min+s*i;
      write_tabular_row(eos,file,"\n");
    }
    Fclose(file);
    
    /* the same, log spaced in h-1 to resolve the crust */
    sprintf(file_name,"%s/%s_tabular_logh.txt",path,eos->description);
    file = Fopen(file_name,"w");
    fprintf(file,"# same columns as %s_tabular.txt, "
                 "h-1 log spaced in [1e-9,h_max-1]\n",eos->description);
    const double l0 = log(1e-9), l1 = log(h_max-1.);
    for (i = 0; i < N; ++i)
    {
      eos->h = 1.+exp(l0+(l1-l0)*i/(N-1.));
      write_tabular_row(eos,file,"\n");
    }
    Fclose(file);
    
    /* the table rows vs the interpolant at the rows and the midpoints */
    sprintf(file_name,"%s/%s_tabular_rows.txt",path,eos->description);
    file = Fopen(file_name,"w");
    fprintf(file,"# [1]:row [2]:h [3]:p_table [4]:e_table "
                 "[5]:rho0_table(n m_b) | "
                 "[6-14]: columns of %s_tabular.txt at h of the row | "
                 "[15-23]: the same at the midpoint (h_i+h_{i+1})/2\n",
                 eos->description);
    const Uint Nrow = eos->spline->sample_size;
    for (i = 0; i < Nrow; ++i)
    {
      fprintf(file,"%u %+0.15e %+0.15e %+0.15e %+0.15e ",i,
        eos->spline->h_sample[i],eos->spline->p_sample[i],
        eos->spline->e_sample[i],eos->spline->rho0_sample[i]);
      eos->h = eos->spline->h_sample[i];
      write_tabular_row(eos,file," ");
      if (i+1 < Nrow)
      {
        eos->h = 0.5*(eos->spline->h_sample[i]+eos->spline->h_sample[i+1]);
        write_tabular_row(eos,file,"\n");
      }
      else
      {
        for (Uint c = 0; c < 9; ++c)
          fprintf(file,"%s",c < 8 ? "nan " : "nan\n");
      }
    }
    Fclose(file);
  }
 
  free_EoS(eos);
  free(path);
  UNUSED(grid);
}
