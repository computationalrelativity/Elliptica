#include "eos_header.h"
#include "maths_spectral_methods_lib.h"
#include "maths_equation_solvings_lib.h"
//#include "maths_calculus_lib.h"

// max number of eos table's row
#define EOS_MAX_NUM_ROW_TABLE (4000)

// h within this relative distance below 1 is the surface (h = 1)
#define EOS_TAB_H_FLOOR_TOL (1e-13)

// max number of eos table's column char
#define EOS_MAX_NUM_COL_TABLE (1000)

void eos_tab_read_table(EoS_T* const eos);
void eos_tab_set_monotone_first_law(EoS_T* const eos);
void eos_tab_free_monotone(void *const mono);

// interpolant for interpolation_method = monotone_first_law
typedef struct EOS_TAB_MONO_T
{
  Uint N;// number of rows
  double *x;// ln(n*m_b) at rows
  double *lp, *le;// ln p, ln e at rows
  double *dp, *de;// d(ln p)/dx, d(ln e)/dx at rows
  double *lh;// ln h at rows
  double *g;// d(ln h)/dx at rows
  double h_max;// h is clamped to this
  // values at the last h (the thermo. functions are called in a row)
  Uint c_valid;
  double c_h, c_h_clamp, c_p, c_e, c_rho0, c_de_dh, c_x_rho0;
}EoS_Tab_Mono_T;
static void mono_eval(const EoS_Tab_Mono_T *const m, const Uint i,
                      const double x, double *const lp, double *const le,
                      double *const dlp, double *const dle);
static double mono_dlnh_dx(const EoS_Tab_Mono_T *const m, const Uint i,
                           const double x);
static double mono_lnh_increment(const EoS_Tab_Mono_T *const m, const Uint i,
                                 const double x);
static void mono_fc_slopes(const double *const x, const double *const y,
                           double *const d, const Uint N);
static void mono_solve(EoS_T *const eos, double h);
static double p_mono(EoS_T *const eos);
static double e_mono(EoS_T *const eos);
static double rho0_mono(EoS_T *const eos);
static double e0_mono(EoS_T *const eos);
static double de_dh_mono(EoS_T *const eos);
static double drho0_dh_mono(EoS_T *const eos);
static double export_rho0_mono(EoS_T *const eos);
static double export_e0_mono(EoS_T *const eos);
