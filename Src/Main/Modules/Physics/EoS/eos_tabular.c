/*
// Andrew Noe & Alireza Rashti
// October 2023
*/

#include "eos_tabular.h"

// Implements thermodynamic functions for a tabular equation of state.

// 1) read thermo. vars. from a given table
// 2) convert to geometrised units
// 3) populate pressure, rest mass density (n*m_b) and total energy density.
//    the enthalpy is set in eos_tab_set_monotone_first_law.
void eos_tab_read_table(EoS_T* const eos)
{
  enum tab_format {
  UTab_format, // undef
  Tab_Format1, // "line,number_density,total_energy_density,pressure"
  Tab_Format2, // "rest_mass_density,specific_internal_energy,pressure"
  };

  Physics_T *const phys = eos->phys;
  enum tab_format format = UTab_format;
  FILE* table_file = Fopen(Gets(P_"table_path"),"r");
  double *p_tab    = alloc_double(EOS_MAX_NUM_ROW_TABLE);
  double *e_tab    = alloc_double(EOS_MAX_NUM_ROW_TABLE);
  double *rho0_tab = alloc_double(EOS_MAX_NUM_ROW_TABLE);
  double *h_tab    = 0;
  Uint num_tab_row = 0;
  double rho0_pnt = DBL_MAX;
  double e_pnt = DBL_MAX;
  double p_pnt = DBL_MAX;
  char line[EOS_MAX_NUM_COL_TABLE] = {'\0'};

  // open and read eos table
  if (strcmp_i(Gets(P_"table_format"),
       "line,number_density,total_energy_density,pressure"))
  {
    format = Tab_Format1;
    num_tab_row = 0;
    while (fgets(line, sizeof(line), table_file))
    {
      Uint l;
      // skip comment lines
      if (line[0] == '#')
      {
        continue;
      }
      // populate thermo. vars.
      // NOTE: it only reads the matching lines.
      if (sscanf(line,"%u %lf %lf %lf",&l, &rho0_pnt, &e_pnt, &p_pnt) == 4)
      {
        rho0_tab[num_tab_row] = rho0_pnt;
        e_tab[num_tab_row]    = e_pnt;
        p_tab[num_tab_row]    = p_pnt;
        num_tab_row++;
        if (num_tab_row >= EOS_MAX_NUM_ROW_TABLE)
        {
          Error0("number of table's row exceeds the max value. "
                 "please increase the macro EOS_MAX_NUM_ROW_TABLE.");
        }
      }
    }
  }
  else if (strcmp_i(Gets(P_"table_format"),
                  "rest_mass_density,specific_internal_energy,pressure"))
  {
    format = Tab_Format2;
    num_tab_row = 0;
    while (fgets(line, sizeof(line), table_file))
    {
      // skip comment lines
      if (line[0] == '#')
      {
        continue;
      }
      // populate thermo. vars.
      if (sscanf(line,"%lf %lf %lf",&rho0_pnt, &e_pnt, &p_pnt) == 3)
      {
        rho0_tab[num_tab_row] = rho0_pnt;
        e_tab[num_tab_row]    = e_pnt;
        p_tab[num_tab_row]    = p_pnt;
        num_tab_row++;
        if (num_tab_row >= EOS_MAX_NUM_ROW_TABLE)
        {
          Error0("number of table's row exceeds the max value. "
                 "please increase the macro EOS_MAX_NUM_ROW_TABLE.");
        }
      }
    }
  }
  else
  {
    Error0(NO_OPTION);
  }
  Fclose(table_file);

  if (num_tab_row == 0)
  {
    Error0("nothing read from the table.");
  }
  // adjust to the actual size
  p_tab = realloc(p_tab,num_tab_row*sizeof(*p_tab)); IsNull(p_tab);
  e_tab = realloc(e_tab,num_tab_row*sizeof(*e_tab)); IsNull(e_tab);
  rho0_tab = realloc(rho0_tab,num_tab_row*sizeof(*rho0_tab)); IsNull(rho0_tab);
  h_tab    = alloc_double(num_tab_row);

  // unit conversion
  if (strcmp_i(eos->unit,"compose") && format == Tab_Format1)
  {
    /* compose format in columns
    *          [line number] [number density] [(total) energy density] [pressure].
    * in units               [1/fm^3]         [g/cm^3]                 [dyn/cm^2].
    * Calculates rest-mass density (total) energy density,
    * and specific enthalpy, and converts to geometrized units.
    *
    * Physical constants from 2018 CODATA values,
    * https://physics.nist.gov/cuu/Constants/Table/allascii.txt
    *
    * G_const = 6.67430E-11;          // Gravitational constant in m^3/(kg s^2)
    * c_const = 299792458;            // Speed of light in m/s
    * M_const = 1.98841E30;           // Solar mass in kg
    * mn_const = 1.67492749804E-27;   // Neutron mass in kg
    * compose data is converted to SI units then geometrized.
    * Pressure conversion factor: G^3 * M_solar^2 / (10 * c^8)
    * Energy density conversion factor: G^3  * M_solar^2 / (10 * c^6)
    * Rest-mass conversion factor: 10^45 * neutron mass * G^3 * m_solar ^2 / (c^6)
    * (Rest-mass is calculated by multiplying the baryon density by the baryon mass.)
    */
    const double p_FACTOR    = 1.80162095578956E-39;
    const double rho0_FACTOR = 0.002712069678583313;
    const double e_FACTOR    = 1.619216164136643E-18;

    for (Uint i = 0; i < num_tab_row; i++)
    {
      p_tab[i]    *= p_FACTOR;
      rho0_tab[i] *= rho0_FACTOR;
      e_tab[i]    *= e_FACTOR;
    }
  }
  else if (strcmp_i(eos->unit,"geo") && format == Tab_Format2)
  {
    for (Uint i = 0; i < num_tab_row; i++)
    {
      e_tab[i] = e_tab[i]*rho0_tab[i]+rho0_tab[i];
    }
  }
  else
  {
    Error0(NO_OPTION);
  }

  eos->spline->h_sample    = h_tab;
  eos->spline->p_sample    = p_tab;
  eos->spline->e_sample    = e_tab;
  eos->spline->rho0_sample = rho0_tab;
  eos->spline->sample_size = num_tab_row;

  h_tab = 0;
  p_tab = 0;
  rho0_tab = 0;
}

// interpolation of tables (monotone_first_law):
// the enthalpy follows from the first law, d(ln h) = dp/(e+p), with h = 1
// at the first row (surface). this holds for any (even finite
// temperature) slice of a table, whereas h = (e+p)/(n*m_b) satisfies it
// only for T = 0; thus rho0 = (e+p)/h and n*m_b is only used for the export.
// everything is a function of the table's own variable x = ln(n*m_b).
// ln p(x) and ln e(x) are monotone cubic (Fritsch-Carlson) Hermite
// interpolants, so p, e are reproduced at all rows and are increasing,
// and ln h(x) = int p/(e+p) d(ln p) (Gauss-Legendre per interval), so
// dp = (e+p) d(ln h) = rho0 dh holds between the rows as well.
// for a given h, x is found by Newton's method. then
//   rho0 = (e+p)/h, de/dh = e (e+p) d(ln e)/dx / (h p d(ln p)/dx),
//   drho0/dh = (de/dh)/h, c_s^2 = dp/de = p d(ln p)/dx / (e d(ln e)/dx) > 0.
// p, e, rho0 are C^1 in h and drho0/dh is C^0.

// 8-point Gauss-Legendre on [0,1]
static const double mono_gl_t[8] = {
  0.0198550717512318841, 0.1016667612931866302,
  0.2372337950418355070, 0.4082826787521750975,
  0.5917173212478249025, 0.7627662049581644930,
  0.8983332387068133698, 0.9801449282487681159 };
static const double mono_gl_w[8] = {
  0.0506142681451881296, 0.1111905172266872352,
  0.1568533229389436437, 0.1813418916891809915,
  0.1813418916891809915, 0.1568533229389436437,
  0.1111905172266872352, 0.0506142681451881296 };

// ln p, ln e and their x derivatives at x in the interval [x_i,x_{i+1}]
static void mono_eval(const EoS_Tab_Mono_T *const m, const Uint i,
                      const double x, double *const lp, double *const le,
                      double *const dlp, double *const dle)
{
  const double dx  = m->x[i+1]-m->x[i];
  const double t   = (x-m->x[i])/dx;
  const double h00 = 2*t*t*t - 3*t*t + 1, h10 = t*t*t - 2*t*t + t;
  const double h01 = -2*t*t*t + 3*t*t,    h11 = t*t*t - t*t;
  const double g00 = 6*t*t - 6*t,         g10 = 3*t*t - 4*t + 1;
  const double g01 = -6*t*t + 6*t,        g11 = 3*t*t - 2*t;

  *lp  = h00*m->lp[i] + h10*dx*m->dp[i] + h01*m->lp[i+1] + h11*dx*m->dp[i+1];
  *le  = h00*m->le[i] + h10*dx*m->de[i] + h01*m->le[i+1] + h11*dx*m->de[i+1];
  *dlp = (g00*m->lp[i] + g01*m->lp[i+1])/dx + g10*m->dp[i] + g11*m->dp[i+1];
  *dle = (g00*m->le[i] + g01*m->le[i+1])/dx + g10*m->de[i] + g11*m->de[i+1];
}

// d(ln h)/dx = p/(e+p) d(ln p)/dx
static double mono_dlnh_dx(const EoS_Tab_Mono_T *const m, const Uint i,
                           const double x)
{
  double lp, le, dlp, dle;
  mono_eval(m,i,x,&lp,&le,&dlp,&dle);
  return dlp / (1. + exp(le - lp));
}

// ln h(x) - ln h(x_i) for x in the interval [x_i,x_{i+1}]
static double mono_lnh_increment(const EoS_Tab_Mono_T *const m, const Uint i,
                                 const double x)
{
  const double a = m->x[i];
  double sum = 0.;
  for (Uint k = 0; k < 8; k++)
  {
    sum += mono_gl_w[k]*mono_dlnh_dx(m,i,a + mono_gl_t[k]*(x-a));
  }
  return sum*(x-a);
}

// Fritsch-Carlson slopes of y(x)
static void mono_fc_slopes(const double *const x, const double *const y,
                           double *const d, const Uint N)
{
  for (Uint i = 0; i < N; i++)
  {
    if (i == 0 || i == N-1)
    {
      const Uint j = (i == 0 ? 0 : N-2);
      d[i] = (y[j+1]-y[j])/(x[j+1]-x[j]);
      continue;
    }
    const double s0 = (y[i]-y[i-1])/(x[i]-x[i-1]);
    const double s1 = (y[i+1]-y[i])/(x[i+1]-x[i]);
    if (s0*s1 <= 0.)
    {
      d[i] = 0.;
    }
    else
    {
      const double w0 = 2.*(x[i+1]-x[i]) + (x[i]-x[i-1]);
      const double w1 = (x[i+1]-x[i]) + 2.*(x[i]-x[i-1]);
      d[i] = (w0+w1)/(w0/s0 + w1/s1);
    }
  }
}

// find x with ln h(x) = lnh, and set the cached thermo. vars. at h
static void mono_solve(EoS_T *const eos, double h)
{
  EoS_Tab_Mono_T *const m = eos->spline->mono;

  if (m->c_valid && h == m->c_h)
  {
    return;
  }
  m->c_valid = 1;
  m->c_h     = h;

  // no matter below the surface
  if (h < 1.-EOS_TAB_H_FLOOR_TOL)
  {
    m->c_p = m->c_e = m->c_rho0 = m->c_de_dh = 0.;
    m->c_x_rho0 = 0.;
    return;
  }
  h = (h < 1.) ? 1. : h;
  h = (h > m->h_max) ? m->h_max : h;
  const double lnh = log(h);

  // interval: lh[i] <= lnh <= lh[i+1]
  Uint lo = 0, hi = m->N-1;
  while (hi - lo > 1)
  {
    const Uint mid = (lo+hi)/2;
    if (m->lh[mid] <= lnh) lo = mid;
    else                   hi = mid;
  }
  const Uint i = lo;

  // safeguarded Newton, F(x) = ln h(x) - lnh is increasing.
  // F is smooth in the interval, so after a Newton step below 1e-8 of the
  // interval the error of the new x is ~1e-16 of it: stop there.
  double a = m->x[i], b = m->x[i+1];
  const double dx_tol = 1e-8*(b-a);
  const double ab_tol = 1e-14*(b-a) + 4.*DBL_EPSILON*fmax(fabs(a),fabs(b));// roundoff
  // initial guess: cubic Hermite of x(ln h) with dx/d(ln h) = 1/g at rows
  const double dlh = m->lh[i+1]-m->lh[i];
  const double u   = (lnh-m->lh[i])/dlh;
  double x = (2*u*u*u - 3*u*u + 1)*a + (-2*u*u*u + 3*u*u)*b +
             (u*u*u - 2*u*u + u)*dlh/m->g[i] + (u*u*u - u*u)*dlh/m->g[i+1];
  if (!(x >= a && x <= b))
  {
    x = a + (b-a)*u;
  }
  Uint it;
  for (it = 0; it < 100; it++)
  {
    const double F = m->lh[i] + mono_lnh_increment(m,i,x) - lnh;
    if (F == 0.)
    {
      break;
    }
    if (F > 0.) b = x;
    else        a = x;
    double xn = x - F/mono_dlnh_dx(m,i,x);
    Flag_T newton = YES;
    if (!(xn > a && xn < b))
    {
      xn = 0.5*(a+b);
      newton = NO;
    }
    const double step = fabs(xn-x);
    x = xn;
    if ((newton == YES && step <= dx_tol) || b-a <= ab_tol)
    {
      break;
    }
  }
  if (it == 100)
  {
    Error0("monotone_first_law: no convergence for x(h).");
  }

  double lp, le, dlp, dle;
  mono_eval(m,i,x,&lp,&le,&dlp,&dle);
  const double p = exp(lp), e = exp(le);
  m->c_p       = p;
  m->c_e       = e;
  m->c_rho0    = (e+p)/h;
  m->c_de_dh   = e*(e+p)*dle / (h*p*dlp);
  m->c_x_rho0  = exp(x);
  m->c_h_clamp = h;
}

static double p_mono(EoS_T *const eos)
{
  mono_solve(eos,eos->h);
  return eos->spline->mono->c_p;
}

static double e_mono(EoS_T *const eos)
{
  mono_solve(eos,eos->h);
  return eos->spline->mono->c_e;
}

static double rho0_mono(EoS_T *const eos)
{
  mono_solve(eos,eos->h);
  return eos->spline->mono->c_rho0;
}

static double e0_mono(EoS_T *const eos)
{
  mono_solve(eos,eos->h);
  const EoS_Tab_Mono_T *const m = eos->spline->mono;
  return (m->c_rho0 > 0.) ? m->c_e/m->c_rho0 - 1. : 0.;
}

static double de_dh_mono(EoS_T *const eos)
{
  mono_solve(eos,eos->h);
  return eos->spline->mono->c_de_dh;
}

static double drho0_dh_mono(EoS_T *const eos)
{
  mono_solve(eos,eos->h);
  const EoS_Tab_Mono_T *const m = eos->spline->mono;
  return (m->c_rho0 > 0.) ? m->c_de_dh/m->c_h_clamp : 0.;
}

// n*m_b of the table, i.e., exp(x)
static double export_rho0_mono(EoS_T *const eos)
{
  mono_solve(eos,eos->h);
  return eos->spline->mono->c_x_rho0;
}

// e/(n*m_b) - 1 with n*m_b of the table
static double export_e0_mono(EoS_T *const eos)
{
  mono_solve(eos,eos->h);
  const EoS_Tab_Mono_T *const m = eos->spline->mono;
  return (m->c_x_rho0 > 0.) ? m->c_e/m->c_x_rho0 - 1. : 0.;
}

// set the interpolant of the table (monotone_first_law).
// the optional parameter enthalpy_ceiling caps h; by default h is capped
// at the last row.
void eos_tab_set_monotone_first_law(EoS_T* const eos)
{
  Physics_T *const phys = eos->phys;
  const Uint N = eos->spline->sample_size;
  const double *const rho0 = eos->spline->rho0_sample;// n*m_b
  const double *const p    = eos->spline->p_sample;
  const double *const e    = eos->spline->e_sample;
  double *const h          = eos->spline->h_sample;
  EoS_Tab_Mono_T *m = calloc(1,sizeof(*m)); IsNull(m);

  if (N < 3)
  {
    Error0("the table needs at least 3 rows.");
  }

  m->N  = N;
  m->x  = alloc_double(N);
  m->lp = alloc_double(N);
  m->le = alloc_double(N);
  m->dp = alloc_double(N);
  m->de = alloc_double(N);
  m->lh = alloc_double(N);
  m->g  = alloc_double(N);
  for (Uint i = 0; i < N; i++)
  {
    if (!(rho0[i] > 0. && e[i] > 0. && p[i] > 0.))
    {
      Error0("the table needs n, e, p > 0 in all rows.");
    }
    m->x[i]  = log(rho0[i]);
    m->lp[i] = log(p[i]);
    m->le[i] = log(e[i]);
    if (i > 0 && !(m->x[i] > m->x[i-1] && m->lp[i] > m->lp[i-1] &&
                   m->le[i] > m->le[i-1]))
    {
      Error0("n, e and p must be strictly increasing in the table.");
    }
  }
  mono_fc_slopes(m->x,m->lp,m->dp,N);
  mono_fc_slopes(m->x,m->le,m->de,N);

  // ln h at the rows, h = 1 at the first row
  m->lh[0] = 0.;
  for (Uint i = 0; i < N-1; i++)
  {
    m->lh[i+1] = m->lh[i] + mono_lnh_increment(m,i,m->x[i+1]);
  }
  for (Uint i = 0; i < N; i++)
  {
    m->g[i] = mono_dlnh_dx(m,(i < N-1 ? i : N-2),m->x[i]);
  }
  for (Uint i = 0; i < N; i++)
  {
    h[i] = exp(m->lh[i]);
  }

  eos->spline->mono = m;
  sprintf(phys->par,"%s_%s",phys->stype,P_"enthalpy_ceiling");
  eos->spline->h_ceil  = (PgetdEZ(phys->par) != DBL_MAX) ?
                         PgetdEZ(phys->par) : h[N-1];
  m->h_max = (eos->spline->h_ceil < h[N-1]) ? eos->spline->h_ceil : h[N-1];

  eos->pressure                 = p_mono;
  eos->energy_density           = e_mono;
  eos->rest_mass_density        = rho0_mono;
  eos->specific_internal_energy = e0_mono;
  eos->de_dh                    = de_dh_mono;
  eos->drho0_dh                 = drho0_dh_mono;
  eos->export_rest_mass_density        = export_rho0_mono;
  eos->export_specific_internal_energy = export_e0_mono;

  // diagnostics, once per run: sound speed and the table's own
  // consistency, (e+p)/h vs n*m_b and p(e) vs a first law integration.
  static int printed = 0;
  if (!printed)
  {
    printed = 1;
    const Uint S = 16;// samples per interval
    double cs2_min = DBL_MAX, cs2_max = 0.;
    double h_cs2_min = 0., h_cs2_max = 0., h_acausal = 0.;
    double r_min = DBL_MAX, r_max = 0.;
    for (Uint i = 0; i < N; i++)
    {
      const double r = (p[i] + e[i]) / (h[i]*rho0[i]);
      r_min = (r < r_min) ? r : r_min;
      r_max = (r > r_max) ? r : r_max;
    }
    for (Uint i = 0; i < N-1; i++)
    {
      for (Uint k = 0; k <= S; k++)
      {
        const double x = m->x[i] + (m->x[i+1]-m->x[i])*k/S;
        double lp, le, dlp, dle;
        mono_eval(m,i,x,&lp,&le,&dlp,&dle);
        const double cs2 = exp(lp-le)*dlp/dle;
        if (cs2 < cs2_min) { cs2_min = cs2; h_cs2_min = h[i]; }
        if (cs2 > cs2_max) { cs2_max = cs2; h_cs2_max = h[i]; }
        if (cs2 >= 1. && h_acausal == 0.) h_acausal = h[i];
      }
    }
    printf(Pretty0"EoS interpolation monotone_first_law: "
                  "ln p, ln e monotone cubic in ln(n m_b),\n");
    printf(Pretty0"dp = rho0 dh between rows, table rows reproduced.\n");
    printf(Pretty0"h-1 at rows 1 and 2  = %e, %e\n",h[1]-1.,h[2]-1.);
    printf(Pretty0"h at the last row    = %0.15f\n",h[N-1]);
    printf(Pretty0"(e+p)/(h n m_b): first row = %0.10f, last row = %0.10f\n",
      (p[0]+e[0])/(h[0]*rho0[0]),(p[N-1]+e[N-1])/(h[N-1]*rho0[N-1]));
    printf(Pretty0"(e+p)/(h n m_b): min = %0.10f, max = %0.10f\n",r_min,r_max);
    printf(Pretty0"c_s^2 min = %e (h = %0.10f), max = %e (h = %0.10f)\n",
           cs2_min,h_cs2_min,cs2_max,h_cs2_max);
    if (cs2_min <= 0.)
    {
      printf(Pretty0"WARNING: c_s^2 <= 0 in the table's range.\n");
    }
    if (h_acausal > 0.)
    {
      printf(Pretty0"WARNING: c_s^2 >= 1 from h = %0.10f on.\n",h_acausal);
    }
  }
}

// free the interpolant of the table
void eos_tab_free_monotone(void *const mono)
{
  EoS_Tab_Mono_T *const m = mono;
  if (!m)
    return;
  Free(m->x);
  Free(m->lp);
  Free(m->le);
  Free(m->dp);
  Free(m->de);
  Free(m->lh);
  Free(m->g);
  free(m);
}
