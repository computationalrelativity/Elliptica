#ifndef physics_EoS_LIB_H
#define physics_EoS_LIB_H
#include "elliptica_system_lib.h"

#define EOS_MAX_STR (400)

/* forward declaration */
struct PHYSICS_T;

/* struct for equation of states */
typedef struct EquationOfState_T
{
 struct PHYSICS_T *phys;
 char description[EOS_MAX_STR];
 char type[EOS_MAX_STR];
 char unit[EOS_MAX_STR];
 
 /* some flags for extra checks:
 // the default value is 0, so initially they are inactive, 
 // but when someone puts them 1, they become active. */
 Uint enthalpy_fatal: 1;/* checks for enthalpy values */
 double (*pressure)(struct EquationOfState_T *const eos);
 double (*energy_density)(struct EquationOfState_T *const eos);/* total energy density */
 double (*rest_mass_density)(struct EquationOfState_T *const eos);
 double (*specific_internal_energy)(struct EquationOfState_T *const eos);
 double (*de_dh)(struct EquationOfState_T *const eos);/* d(energy_density)/dh */
 double (*drho0_dh)(struct EquationOfState_T *const eos);/* d(rest_mass_density)/dh */
 /* rest mass density and specific internal energy written to the ID for
 // the evolution code, i.e., n*m_b of the table. they differ from
 // rest_mass_density and specific_internal_energy for a table, whose
 // enthalpy is defined by the first law (see eos_tabular.c). */
 double (*export_rest_mass_density)(struct EquationOfState_T *const eos);
 double (*export_specific_internal_energy)(struct EquationOfState_T *const eos);

 ///////////////
 // polytrope //
 ///////////////
 double *K;/* polytropic constant */
 double *rho0_th;/* thresholds of rest mass density */
 double *h_th;/* enthalpy thresholds */
 double *n;/* polytropic index n = 1/(gamma-1) */
 double *gamma;/* polytropic index */
 double *a;/* constant coefficient to ensure continuity */
 double h;/* specific enthalpy, h = (total_energy_density + pressure)/rest_mass_density */
 Uint N;/* number of intervals i.e number of pieces */
 
 //////////////
 // discrete //
 //////////////
 /* spline interpolation when eos format is discrete like table */
 struct
 {
   double h_floor;/* set h to h_floor if h < h_floor */
   double h_ceil;/* set h to h_ceil if h > h_ceil */
   Uint sample_size;/* the number of points for a spline fit */
   double *h_sample;/* enthalpy sample values */
   double *p_sample;/* pressure sample values */
   double *e_sample;/* energy_density sample values */
   double *rho0_sample;/* rest_mass_density sample values,
                       // for a table this is n*m_b of the table
                       // (used for the export). */
   void *interp_p;/* spline interpolation struct for p */
   void *interp_e;/* spline interpolation struct for e */
   void *interp_rho0;/* spline interpolation struct for rho0 */
   struct EOS_TAB_MONO_T *mono;/* interpolation_method = monotone_first_law */
 }spline[1];
 
}EoS_T;
#undef EOS_MAX_STR

EoS_T *init_EoS(struct PHYSICS_T *const phys);
void free_EoS(EoS_T *eos);
void test_EoS(struct PHYSICS_T *const phys);


#endif


