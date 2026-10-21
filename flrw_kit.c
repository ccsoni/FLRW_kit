/* FLRW_kit: implementation
 *
 * Background
 *   R(a) = a^4 E(a)^2 = Og + rho4_nu(a) + Ocb a + Ok a^2 + OL a^4,
 *   rho4_nu(a) = a^4 rho_nu(a) / rho_c0 = Onu1 sum_i F(y_i),   y_i = m_i a / (k_B T_nu0),
 *   where Onu1 = Og(T_cmb) (7/8) (4/11)^{4/3} N_eff/3 is one massless species and F(y) is the
 *   Fermi-Dirac energy density in units of its massless value (flrw_nu_table.h; the pressure
 *   P(y) is tabulated in the same way, P(0) = 1/3).
 *   With u = sqrt(a), which keeps the integrands smooth at a -> 0 with or without radiation,
 *     t(a)   = int_0^sqrt(a) 2u^3 du / sqrt(R(u^2)),   tau(a) = int_0^sqrt(a) 2u du / sqrt(R(u^2)),
 *     DC(a)  = int_sqrt(a)^1 (same integrand as tau),
 *   evaluated by adaptive Gauss-Kronrod (7-15 points).
 *
 * Linear growth (x = ln a, ' = d/dx, Newtonian sub-horizon equations)
 *   d_cb'' + (2 + dlnE/dx) d_cb' = (3/2) [Ocb(a) d_cb + sum_i On_i(a) d_i]
 *   d_i''  + (2 + dlnE/dx) d_i'  = (3/2) [Ocb(a) d_cb + sum_j On_j(a) d_j] - J_i d_i
 *     On_i(a) = (rho_i - 3 p_i)/rho_crit(a): matter-like (clustering) part of species i,
 *     J_i = c_s,i^2 k^2 / (a E)^2, c_s,i^2 = c_nr^2 / (1 + 3 c_nr^2), c_nr = sqrt(25 zeta(5) /
 *     (3 zeta(3))) k_B T_nu0 / (m_i a) (fluid approximation, cf. Shoji & Komatsu 2010).
 *   Small-scale limit (cb): d_i = 0. Large-scale limit (m): d_i = d_cb (J_i = 0).
 *   When J_i > J_OFF the species is treated quasi-statically (d_i'' and d_i' neglected), and it
 *   becomes dynamical again when J_i < J_ON; this avoids integrating fast acoustic oscillations.
 *   Initial conditions at a_i << a_eq: Meszaros growing mode d = a + (2/3) a_eq, d' = a, with
 *   a_eq = (Og + Onu,rel)/Ocb; dynamical neutrinos start with d_i = d_cb.
 *
 * Spherical collapse of cdm + baryons with smooth photons, neutrinos and OL
 *   With y = R / (a r_L) (r_L: Lagrangian radius, 1 + d = y^-3),
 *     y'' + (2 + dlnE/dx) y' = -(1/2) Ocb(a) (y^-2 - y),
 *   integrated for w = y - 1 to keep the precision of the initial contrast ~1e-9,
 *   started on the linear growing mode at a_i. The initial amplitude is found such that y -> 0
 *   at a_c, and delta_c is that amplitude evolved with the linear equation (cb limit).
 *   Turnaround is where (a y)' = 0. The virial radius y_v R_ta follows from energy conservation
 *   between turnaround and collapse and the virial theorem at collapse, with the potential energy
 *   of the smooth components ~ (rho_s + 3 p_s) R^2 evaluated at each epoch:
 *     4 eta_c y_v^3 - 2 (1 + eta_ta) y_v + 1 = 0,
 *     eta(a) = -S(a) R_ta^3 / (2 Ocb),   S = sum_smooth (rho_s + 3 p_s)/rho_c0  (= -2 OL for OL)
 *   (for OL alone eta_c = eta_ta = OL R_ta^3 / Ocb, the classical result).
 */
#include "flrw_kit.h"

#include <float.h>
#include <math.h>
#include <string.h>

#include "flrw_nu_table.h"

#define INT_EPSREL 1e-13
#define INT_MAXSEG 400
#define ROOT_MAXIT 200
#define A_MAX 1e12 /* upper bound of a when inverting t(a), tau(a) */
#define PI 3.14159265358979323846

/* physical constants */
#define OMEGA_GAMMA_H2_REF 2.472975333e-5 /* Og h^2 for T_cmb = 2.7255 K (CODATA 2018) */
#define T_CMB_REF 2.7255
#define K_B_EV 8.617333262e-5              /* Boltzmann constant [eV/K] */
#define ZETA3 1.2020569031595942
#define ZETA5 1.0369277551433699
#define ZETA7 1.0083492773819228
#define NU_CS_COEF 2.68115000876075        /* sqrt(25 zeta(5) / (3 zeta(3))) */

/* neutrino mass splittings [eV^2] (NuFIT 5.x); can be overridden at compile time */
#ifndef FLRW_DM21_SQ
#define FLRW_DM21_SQ 7.42e-5
#endif
#ifndef FLRW_DM31_SQ_NORMAL
#define FLRW_DM31_SQ_NORMAL 2.515e-3   /* m3^2 - m1^2, normal ordering */
#endif
#ifndef FLRW_DM32_SQ_INVERTED
#define FLRW_DM32_SQ_INVERTED 2.498e-3 /* m2^2 - m3^2, inverted ordering */
#endif

/* growth and collapse ODEs */
#define A_INIT 1e-9     /* initial scale factor (times a if a < 1) */
#define ODE_RTOL 1e-11
#define ODE_NMAX 8
#define ODE_MAXSTEP 1000000
#define J_ON 100.0      /* neutrino fluid: quasi-static -> dynamical below this J */
#define J_OFF 400.0     /* neutrino fluid: dynamical -> quasi-static above this J */
#define Y_COLLAPSE 1e-4 /* R / (a r_L) at which collapse is extrapolated */

typedef struct {
  double ocb, ol, ok, og; /* Ocb, OL, Ok, Og today */
  double onu1;            /* Omega of one massless neutrino species today */
  int nm;                 /* number of massive species */
  double nless;           /* number of massless species */
  double y0[3];           /* m_i / (k_B T_nu0) of the massive species */
  double cs0[3];          /* non-relativistic sound speed at a = 1 [c] (c_nr = cs0/a) */
  double onu, om;         /* Onu (all species) and Om (Ocb + matter-like Onu) today */
  double aeq;             /* (Og + Onu with all species relativistic) / Ocb */
} Bg;

typedef double (*IntFn)(double x, const void* par);

/* ------------------------------------------------------------------ */
/* Adaptive Gauss-Kronrod 7-15 point quadrature                       */
/* ------------------------------------------------------------------ */
static const double xgk[8] = {
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
static const double wgk[8] = {
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
static const double wg[4] = {
    0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

static double gk15(IntFn f, const void* par, double a, double b, double* err) {
  const double c = 0.5 * (a + b), h = 0.5 * (b - a);
  const double fc = f(c, par);
  double rk = fc * wgk[7], rg = fc * wg[3];
  for (int j = 0; j < 7; ++j) {
    const double dx = h * xgk[j];
    const double fs = f(c - dx, par) + f(c + dx, par);
    rk += wgk[j] * fs;
    if (j & 1) rg += wg[j / 2] * fs;
  }
  *err = fabs((rk - rg) * h);
  return rk * h;
}

/* Globally adaptive: repeatedly bisect the subinterval with the largest error until the total
 * error falls below the tolerance or the round-off limit is reached
 * (same strategy as QUADPACK QAG) */
static double integrate(IntFn f, const void* par, double a, double b) {
  if (a == b) return 0.0;
  double lo[INT_MAXSEG], hi[INT_MAXSEG], val[INT_MAXSEG], err[INT_MAXSEG];
  int n = 1;
  lo[0] = a;
  hi[0] = b;
  val[0] = gk15(f, par, a, b, &err[0]);
  double sum = val[0], esum = err[0];
  while (isfinite(sum) && n < INT_MAXSEG) {
    if (esum <= fmax(INT_EPSREL * fabs(sum), 1e3 * DBL_MIN)) break;
    int k = 0;
    for (int i = 1; i < n; ++i)
      if (err[i] > err[k]) k = i;
    if (err[k] <= 50.0 * DBL_EPSILON * fabs(val[k])) break; /* round-off limit */
    const double m = 0.5 * (lo[k] + hi[k]);
    if (!(m > fmin(lo[k], hi[k]) && m < fmax(lo[k], hi[k]))) break;
    double el, er;
    const double vl = gk15(f, par, lo[k], m, &el);
    const double vr = gk15(f, par, m, hi[k], &er);
    sum += vl + vr - val[k];
    esum += el + er - err[k];
    lo[n] = m, hi[n] = hi[k], val[n] = vr, err[n] = er;
    hi[k] = m, val[k] = vl, err[k] = el;
    ++n;
  }
  /* re-sum at the end to avoid accumulated round-off in the running sum */
  sum = 0.0;
  for (int i = 0; i < n; ++i) sum += val[i];
  return sum;
}

/* ------------------------------------------------------------------ */
/* Brent's method (f(a), f(b) of opposite signs; xtol: absolute tolerance)  */
/* ------------------------------------------------------------------ */
typedef double (*RootFn)(double x, const void* par);

static double brent(RootFn f, const void* par, double a, double b, double fa, double fb,
                    double xtol) {
  if (!(fa * fb <= 0.0)) return NAN;
  double c = b, fc = fb, d = 0.0, e = 0.0;
  for (int it = 0; it < ROOT_MAXIT; ++it) {
    if ((fb > 0.0) == (fc > 0.0)) {
      c = a;
      fc = fa;
      e = d = b - a;
    }
    if (fabs(fc) < fabs(fb)) {
      a = b, b = c, c = a;
      fa = fb, fb = fc, fc = fa;
    }
    const double tol = 2.0 * DBL_EPSILON * fabs(b) + 0.5 * xtol;
    const double xm = 0.5 * (c - b);
    if (fabs(xm) <= tol || fb == 0.0) return b;
    if (fabs(e) >= tol && fabs(fa) > fabs(fb)) {
      double p, q, r;
      const double s = fb / fa;
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) q = -q;
      p = fabs(p);
      if (2.0 * p < fmin(3.0 * xm * q - fabs(tol * q), fabs(e * q))) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    b += (fabs(d) > tol) ? d : copysign(tol, xm);
    fb = f(b, par);
    if (!isfinite(fb)) return NAN;
  }
  return b;
}


/* ------------------------------------------------------------------ */
/* Neutrino energy density and pressure                               */
/* ------------------------------------------------------------------ */
/* quintic Hermite interpolation from values and first/second derivatives (v[0..2]) */
static inline double herm5(const double* v0, const double* v1, double h, double s) {
  const double s2 = s * s, s3 = s2 * s, s4 = s3 * s, s5 = s4 * s;
  const double H0 = 1.0 - 10.0 * s3 + 15.0 * s4 - 6.0 * s5;
  const double H1 = s - 6.0 * s3 + 8.0 * s4 - 3.0 * s5;
  const double H2 = 0.5 * (s2 - 3.0 * s3 + 3.0 * s4 - s5);
  const double G0 = 10.0 * s3 - 15.0 * s4 + 6.0 * s5;
  const double G1 = -4.0 * s3 + 7.0 * s4 - 3.0 * s5;
  const double G2 = 0.5 * (s3 - 2.0 * s4 + s5);
  return v0[0] * H0 + v1[0] * G0 + h * (v0[1] * H1 + v1[1] * G1) +
         h * h * (v0[2] * H2 + v1[2] * G2);
}

/* energy density F(y) and pressure P(y) of one species in units of its massless energy density,
 * y = m / (k_B T_nu): small-y series, table, or large-y expansion */
static void nu_FP(double y, double* F, double* P) {
  if (y < 1e-3) {
    const double x = y * y / (PI * PI);
    *F = 1.0 + 5.0 / 7.0 * x;
    *P = 1.0 / 3.0 - 5.0 / 21.0 * x;
  } else if (y > 1e4) {
    const double c = 120.0 / (7.0 * PI * PI * PI * PI);
    const double iy = 1.0 / y, iy3 = iy * iy * iy;
    *F = c * (1.5 * ZETA3 * y + 11.25 * ZETA5 * iy - 88.59375 * ZETA7 * iy3);
    *P = c * (7.5 * ZETA5 * iy - 118.125 * ZETA7 * iy3);
  } else {
    const double u = (log(y) - NU_TAB_T0) / NU_TAB_H;
    int i = (int)u;
    if (i > NU_TAB_N - 2) i = NU_TAB_N - 2;
    if (i < 0) i = 0;
    const double s = u - i;
    *F = exp(herm5(nu_tab[i], nu_tab[i + 1], NU_TAB_H, s));
    *P = exp(herm5(nu_tab[i] + 3, nu_tab[i + 1] + 3, NU_TAB_H, s));
  }
}

/* a^4 rho_nu/rho_c0, a^4 p_nu/rho_c0 and the matter-like part a^4 (rho_nu - 3 p_nu)/rho_c0;
 * F, P of each massive species are stored in Fi, Pi unless NULL */
static void nu_sum(const Bg* b, double a, double* rho4, double* p4, double* mat4, double* Fi,
                   double* Pi) {
  double r = b->nless, p = b->nless / 3.0, m = 0.0;
  for (int i = 0; i < b->nm; ++i) {
    double F, P;
    nu_FP(b->y0[i] * a, &F, &P);
    r += F;
    p += P;
    m += F - 3.0 * P;
    if (Fi) {
      Fi[i] = F;
      Pi[i] = P;
    }
  }
  *rho4 = b->onu1 * r;
  *p4 = b->onu1 * p;
  *mat4 = b->onu1 * m;
}

/* ------------------------------------------------------------------ */
/* Background                                                         */
/* ------------------------------------------------------------------ */
static int make_bg(const cosm_param* cp, Bg* b) {
  if (!cp) return -1;
  const double h = cp->hubble, T = cp->T_cmb, neff = cp->N_eff;
  b->ocb = cp->omega_cdm + cp->omega_b;
  b->ol = cp->omega_de;
  if (!(b->ocb > 0.0) || !isfinite(b->ocb) || !isfinite(b->ol)) return -1;
  if (!(h > 0.0) || !isfinite(h) || !(neff >= 0.0) || !isfinite(neff)) return -1;
  if (!isfinite(cp->omega_r)) return -1;
  const int need_T = (cp->omega_r < 0.0) || (neff > 0.0);
  if (need_T && (!(T > 0.0) || !isfinite(T))) return -1;

  const double og_T = need_T ? OMEGA_GAMMA_H2_REF * pow(T / T_CMB_REF, 4.0) / (h * h) : 0.0;
  b->og = (cp->omega_r >= 0.0) ? cp->omega_r : og_T;
  b->onu1 = og_T * 0.875 * pow(4.0 / 11.0, 4.0 / 3.0) * neff / 3.0;
  const double kT = K_B_EV * cbrt(4.0 / 11.0) * pow(neff / 3.0, 0.25) * T; /* k_B T_nu0 [eV] */
  b->nm = 0;
  b->nless = 0.0;
  for (int i = 0; i < 3; ++i) {
    const double m = cp->m_nu[i];
    if (!(m >= 0.0) || !isfinite(m)) return -1;
    if (neff == 0.0) continue; /* no neutrinos at all */
    if (m == 0.0) {
      b->nless += 1.0;
    } else {
      b->y0[b->nm] = m / kT;
      b->cs0[b->nm] = NU_CS_COEF * kT / m;
      ++b->nm;
    }
  }
  double rho4, p4, mat4;
  nu_sum(b, 1.0, &rho4, &p4, &mat4, NULL, NULL);
  b->onu = rho4;
  b->om = b->ocb + mat4;
  b->ok = 1.0 - b->ocb - b->og - b->onu - b->ol;
  b->aeq = (b->og + b->onu1 * (b->nless + b->nm)) / b->ocb;
  return 0;
}

/* R(a) = a^4 E(a)^2 */
static double R4_a(const Bg* b, double a) {
  double rho4, p4, mat4;
  nu_sum(b, a, &rho4, &p4, &mat4, NULL, NULL);
  return b->og + rho4 + a * (b->ocb + a * (b->ok + b->ol * a * a));
}

static inline double E_a(const Bg* b, double a) { return sqrt(R4_a(b, a)) / (a * a); }

/* integrands in u = sqrt(a) (never evaluated at u = 0) */
static double f_time(double u, const void* p) {
  const double a = u * u;
  return 2.0 * a / sqrt(R4_a((const Bg*)p, a) / a);
}
static double f_conf(double u, const void* p) {
  const double a = u * u;
  return 2.0 / sqrt(R4_a((const Bg*)p, a) / a);
}

static double time_a(const Bg* b, double a) { return integrate(f_time, b, 0.0, sqrt(a)); }
static double conf_a(const Bg* b, double a) { return integrate(f_conf, b, 0.0, sqrt(a)); }
static double DC_a(const Bg* b, double a) { return integrate(f_conf, b, sqrt(a), 1.0); }

static double DM_from_DC(const Bg* b, double dc) {
  const double x = b->ok * dc * dc;
  if (fabs(x) < 1e-6) return dc * (1.0 + x / 6.0 * (1.0 + x / 20.0));
  const double sk = sqrt(fabs(b->ok));
  return (b->ok > 0.0) ? sinh(sk * dc) / sk : sin(sk * dc) / sk;
}

static inline double a_of_z(double z) { return (z > -1.0) ? 1.0 / (1.0 + z) : NAN; }

/* quantities needed by the growth and collapse equations at scale factor a */
typedef struct {
  double R;      /* a^4 E^2 */
  double dlnE;   /* d ln E / d ln a */
  double ocb;    /* Ocb(a) = rho_cb / rho_crit(a) */
  double on[3];  /* (rho_i - 3 p_i) / rho_crit(a) of each massive species */
  double onsum;  /* sum of on[] */
  double smooth; /* a^4 sum (rho + 3p) / rho_c0 over the smooth components (g, nu, OL) */
} Env;

static void env_at(const Bg* b, double a, Env* e) {
  double rho4, p4, mat4, F[3], P[3];
  nu_sum(b, a, &rho4, &p4, &mat4, F, P);
  const double a2 = a * a;
  e->R = b->og + rho4 + a * b->ocb + a2 * (b->ok + b->ol * a2);
  e->dlnE = (-4.0 * b->og - 3.0 * (rho4 + p4) - 3.0 * b->ocb * a - 2.0 * b->ok * a2) / (2.0 * e->R);
  e->ocb = b->ocb * a / e->R;
  e->onsum = 0.0;
  for (int i = 0; i < b->nm; ++i) {
    e->on[i] = b->onu1 * (F[i] - 3.0 * P[i]) / e->R;
    e->onsum += e->on[i];
  }
  e->smooth = 2.0 * b->og + rho4 + 3.0 * p4 - 2.0 * b->ol * a2 * a2;
}

/* ------------------------------------------------------------------ */
/* Inverse of t(a), tau(a) (safeguarded Newton in x = ln a)          */
/* ------------------------------------------------------------------ */
typedef struct {
  const Bg* b;
  int conformal;
  double target;
} InvPar;

static double inv_f(double x, const void* p) {
  const InvPar* ip = p;
  const double a = exp(x);
  return (ip->conformal ? conf_a(ip->b, a) : time_a(ip->b, a)) - ip->target;
}

/* d t/d ln a = 1/E,  d tau/d ln a = 1/(a E) */
static double inv_df(double x, const InvPar* ip) {
  const double a = exp(x);
  const double E = E_a(ip->b, a);
  return ip->conformal ? 1.0 / (a * E) : 1.0 / E;
}

static double invert(const Bg* b, int conformal, double target) {
  if (!(target > 0.0) || !isfinite(target)) return NAN;
  const InvPar ip = {b, conformal, target};
  /* initial guess from the matter-dominated solution
   * t = (2/3) a^{3/2}/sqrt(Ocb), tau = 2 sqrt(a)/sqrt(Ocb) */
  const double so = sqrt(b->ocb);
  double x = conformal ? 2.0 * log(0.5 * so * target) : (2.0 / 3.0) * log(1.5 * so * target);
  x = fmin(x, log(A_MAX));
  double fx = inv_f(x, &ip);
  if (!isfinite(fx)) return NAN;

  /* bracket the root */
  double lo, hi, flo, fhi;
  if (fx < 0.0) {
    lo = x, flo = fx;
    hi = x, fhi = fx;
    while (fhi < 0.0) {
      lo = hi, flo = fhi;
      hi += 1.0;
      if (hi > log(A_MAX)) return NAN; /* e.g. beyond the maximum conformal time (OL > 0) */
      fhi = inv_f(hi, &ip);
      if (!isfinite(fhi)) return NAN;
    }
  } else {
    hi = x, fhi = fx;
    lo = x, flo = fx;
    while (flo > 0.0) {
      hi = lo, fhi = flo;
      lo -= 1.0;
      flo = inv_f(lo, &ip);
      if (!isfinite(flo)) return NAN;
    }
  }
  if (flo == 0.0) return exp(lo);
  if (fhi == 0.0) return exp(hi);

  /* Newton (bisection when the step leaves the bracket) */
  x = 0.5 * (lo + hi);
  fx = inv_f(x, &ip);
  for (int it = 0; it < ROOT_MAXIT; ++it) {
    if (fx < 0.0) lo = x;
    else hi = x;
    double xn = x - fx / inv_df(x, &ip);
    if (!(xn > lo && xn < hi)) xn = 0.5 * (lo + hi);
    const double dx = xn - x;
    x = xn;
    if (fabs(dx) < 4.0 * DBL_EPSILON * fmax(1.0, fabs(x))) break;
    fx = inv_f(x, &ip);
    if (!isfinite(fx)) return NAN;
    if (fx == 0.0) break;
  }
  return exp(x);
}

/* ------------------------------------------------------------------ */
/* Adaptive Dormand-Prince 5(4) ODE integrator                         */
/* ------------------------------------------------------------------ */
typedef void (*OdeFn)(double x, const double* y, double* dy, void* par);

typedef struct {
  OdeFn f;
  void* par;
  int n;
  double x, h, atol;
  double y[ODE_NMAX], k1[ODE_NMAX]; /* state and its derivative (first-same-as-last) */
  long nstep;
} Ode;

static void ode_reset(Ode* o) { o->f(o->x, o->y, o->k1, o->par); }

static void ode_init(Ode* o, OdeFn f, void* par, int n, double x, const double* y, double h,
                     double atol) {
  o->f = f;
  o->par = par;
  o->n = n;
  o->x = x;
  o->h = h;
  o->atol = atol;
  o->nstep = 0;
  memcpy(o->y, y, (size_t)n * sizeof(double));
  ode_reset(o);
}

/* one trial step of size h from (o->x, o->y): fifth-order result in y5, f(x+h, y5) in k7;
 * returns the scaled error norm (NaN if the right-hand side failed) */
static double ode_try(const Ode* o, double h, double* y5, double* k7) {
  static const double a21 = 1.0 / 5.0;
  static const double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
  static const double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
  static const double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0,
                      a54 = -212.0 / 729.0;
  static const double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0,
                      a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;
  static const double b1 = 35.0 / 384.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0,
                      b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0;
  static const double e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0,
                      e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;
  const int n = o->n;
  const double x = o->x, *y = o->y, *k1 = o->k1;
  double k2[ODE_NMAX], k3[ODE_NMAX], k4[ODE_NMAX], k5[ODE_NMAX], k6[ODE_NMAX], t[ODE_NMAX] = {0.0};
  for (int i = 0; i < n; ++i) t[i] = y[i] + h * a21 * k1[i];
  o->f(x + h / 5.0, t, k2, o->par);
  for (int i = 0; i < n; ++i) t[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i]);
  o->f(x + 0.3 * h, t, k3, o->par);
  for (int i = 0; i < n; ++i) t[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
  o->f(x + 0.8 * h, t, k4, o->par);
  for (int i = 0; i < n; ++i)
    t[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
  o->f(x + h * 8.0 / 9.0, t, k5, o->par);
  for (int i = 0; i < n; ++i)
    t[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
  o->f(x + h, t, k6, o->par);
  for (int i = 0; i < n; ++i)
    y5[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);
  o->f(x + h, y5, k7, o->par);
  double err = 0.0;
  for (int i = 0; i < n; ++i) {
    const double ei =
        h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * k7[i]);
    const double sc = o->atol + ODE_RTOL * fmax(fabs(y[i]), fabs(y5[i]));
    err += (ei / sc) * (ei / sc);
  }
  return sqrt(err / n);
}

/* advance by one accepted step without passing x_end (> o->x); returns 0, or -1 on failure */
static int ode_step(Ode* o, double x_end) {
  double y5[ODE_NMAX], k7[ODE_NMAX];
  if (x_end - o->x <= 1e-13 * fmax(1.0, fabs(o->x))) { /* negligible remainder: Euler */
    for (int i = 0; i < o->n; ++i) o->y[i] += (x_end - o->x) * o->k1[i];
    o->x = x_end;
    ode_reset(o);
    return 0;
  }
  for (;;) {
    if (++o->nstep > ODE_MAXSTEP) return -1;
    double h = o->h;
    int last = 0;
    if (o->x + h >= x_end) {
      h = x_end - o->x;
      last = 1;
    }
    if (!(h > 1e-15 * fmax(1.0, fabs(o->x)))) return -1;
    const double err = ode_try(o, h, y5, k7);
    if (err <= 1.0) {
      o->x = last ? x_end : o->x + h;
      memcpy(o->y, y5, (size_t)o->n * sizeof(double));
      memcpy(o->k1, k7, (size_t)o->n * sizeof(double));
      const double fac = (err > 0.0) ? fmin(5.0, fmax(0.2, 0.9 * pow(err, -0.2))) : 5.0;
      if (!last || fac < 1.0) o->h = h * fac;
      return 0;
    }
    o->h = h * (isfinite(err) ? fmax(0.1, 0.9 * pow(err, -0.2)) : 0.1);
  }
}

/* ------------------------------------------------------------------ */
/* Linear growth                                                      */
/* ------------------------------------------------------------------ */
enum { GROW_CB = 0, GROW_M = 1 };

typedef struct {
  const Bg* b;
  int mode;
} GrowPar;

static void grow_rhs(double x, const double* y, double* dy, void* p) {
  const GrowPar* gp = p;
  Env e;
  env_at(gp->b, exp(x), &e);
  const double src = (gp->mode == GROW_M) ? e.ocb + e.onsum : e.ocb;
  dy[0] = y[1];
  dy[1] = -(2.0 + e.dlnE) * y[1] + 1.5 * src * y[0];
}

static double a_init(double a) { return A_INIT * fmin(a, 1.0); }

/* Meszaros growing mode (d, d') at a_i, normalized to d(a_i) = 1 */
static void grow_init(const Bg* b, double ai, double* y) {
  y[0] = 1.0;
  y[1] = ai / (ai + 2.0 / 3.0 * b->aeq);
}

/* linear contrast with d(a_i) = 1 at a1 <= a2 */
static int growth_lin(const Bg* b, int mode, double ai, double a1, double a2, double* d1,
                      double* d2) {
  GrowPar gp = {b, mode};
  double y0[2];
  grow_init(b, ai, y0);
  Ode o;
  ode_init(&o, grow_rhs, &gp, 2, log(ai), y0, 0.01, 1e-14);
  const double x1 = log(a1), x2 = log(a2);
  while (o.x < x1)
    if (ode_step(&o, x1)) return -1;
  *d1 = o.y[0];
  while (o.x < x2)
    if (ode_step(&o, x2)) return -1;
  *d2 = o.y[0];
  return 0;
}

/* D(a) / D(1) */
static double growth_norm(const Bg* b, int mode, double a) {
  const double lo = fmin(a, 1.0), hi = fmax(a, 1.0);
  double dlo, dhi;
  if (growth_lin(b, mode, a_init(lo), lo, hi, &dlo, &dhi)) return NAN;
  return (a <= 1.0) ? dlo / dhi : dhi / dlo;
}

static double growth_cb_a(const Bg* b, double a) { return growth_norm(b, GROW_CB, a); }
static double growth_m_a(const Bg* b, double a) { return growth_norm(b, GROW_M, a); }

/* ---- scale-dependent growth: cb + neutrino fluids ---- */
typedef struct {
  const Bg* b;
  double k2;
  int qs[3]; /* species treated quasi-statically */
} KPar;

/* J_i and the quasi-static contrast of each species for state y */
static void k_aux(const KPar* kp, double a, const Env* e, const double* y, double* J,
                  double* dqs) {
  const Bg* b = kp->b;
  const double aE2 = e->R / (a * a);
  double s = e->ocb * y[0];
  for (int i = 0; i < b->nm; ++i)
    if (!kp->qs[i]) s += e->on[i] * y[2 + 2 * i];
  for (int i = 0; i < b->nm; ++i) {
    const double c2 = (b->cs0[i] / a) * (b->cs0[i] / a);
    J[i] = c2 / (1.0 + 3.0 * c2) * kp->k2 / aE2;
    dqs[i] = 1.5 * s / (J[i] - 1.5 * e->on[i]);
  }
}

static void k_rhs(double x, const double* y, double* dy, void* p) {
  const KPar* kp = p;
  const double a = exp(x);
  Env e;
  env_at(kp->b, a, &e);
  double J[3], dqs[3];
  k_aux(kp, a, &e, y, J, dqs);
  double s = e.ocb * y[0];
  for (int i = 0; i < kp->b->nm; ++i) s += e.on[i] * (kp->qs[i] ? dqs[i] : y[2 + 2 * i]);
  s *= 1.5;
  const double fr = 2.0 + e.dlnE;
  dy[0] = y[1];
  dy[1] = -fr * y[1] + s;
  for (int i = 0; i < kp->b->nm; ++i) {
    if (kp->qs[i]) {
      dy[2 + 2 * i] = dy[3 + 2 * i] = 0.0;
    } else {
      dy[2 + 2 * i] = y[3 + 2 * i];
      dy[3 + 2 * i] = -fr * y[3 + 2 * i] + s - J[i] * y[2 + 2 * i];
    }
  }
}

/* cdm + baryon and total matter contrasts at a for wavenumber k (d(a_i) = 1) */
static int growth_k(const Bg* b, double k, double ai, double a, double* dcb, double* dm) {
  KPar kp = {b, k * k, {0, 0, 0}};
  const int n = 2 + 2 * b->nm;
  double y[ODE_NMAX] = {0.0}, J[3], dqs[3], dqs_prev[3];
  grow_init(b, ai, y);
  Env e;
  env_at(b, ai, &e);
  k_aux(&kp, ai, &e, y, J, dqs);
  for (int i = 0; i < b->nm; ++i) {
    kp.qs[i] = !(J[i] < J_ON);
    y[2 + 2 * i] = y[0]; /* adiabatic: same number-density contrast as cdm + baryons */
    y[3 + 2 * i] = y[1];
  }
  k_aux(&kp, ai, &e, y, J, dqs_prev);
  double x_prev = log(ai);
  Ode o;
  ode_init(&o, k_rhs, &kp, n, log(ai), y, 0.01, 1e-14);
  const double xe = log(a);
  while (o.x < xe) {
    if (ode_step(&o, xe)) return -1;
    const double an = exp(o.x);
    env_at(b, an, &e);
    k_aux(&kp, an, &e, o.y, J, dqs);
    int changed = 0;
    for (int i = 0; i < b->nm; ++i) {
      if (kp.qs[i] && J[i] < J_ON) {
        o.y[2 + 2 * i] = dqs[i];
        o.y[3 + 2 * i] = (dqs[i] - dqs_prev[i]) / (o.x - x_prev);
        kp.qs[i] = 0;
        changed = 1;
      } else if (!kp.qs[i] && J[i] > J_OFF) {
        kp.qs[i] = 1;
        changed = 1;
      }
    }
    if (changed) {
      ode_reset(&o);
      k_aux(&kp, an, &e, o.y, J, dqs);
    }
    memcpy(dqs_prev, dqs, sizeof dqs);
    x_prev = o.x;
  }
  env_at(b, a, &e);
  k_aux(&kp, a, &e, o.y, J, dqs);
  double num = e.ocb * o.y[0], den = e.ocb;
  for (int i = 0; i < b->nm; ++i) {
    num += e.on[i] * (kp.qs[i] ? dqs[i] : o.y[2 + 2 * i]);
    den += e.on[i];
  }
  *dcb = o.y[0];
  *dm = num / den;
  return 0;
}

/* D_cb(k, a) or D_m(k, a), normalized by the large-scale limit at a = 1 */
static double growth_k_norm(const Bg* b, double k, double a, int total) {
  if (!(k >= 0.0) || !isfinite(k)) return NAN;
  const double ai = a_init(fmin(a, 1.0));
  double dcb, dm, d1, dummy;
  if (growth_k(b, k, ai, a, &dcb, &dm)) return NAN;
  if (growth_lin(b, GROW_M, ai, 1.0, 1.0, &d1, &dummy)) return NAN;
  return (total ? dm : dcb) / d1;
}

static double growth_cb_k_a(const Bg* b, double k, double a) { return growth_k_norm(b, k, a, 0); }
static double growth_m_k_a(const Bg* b, double k, double a) { return growth_k_norm(b, k, a, 1); }

/* ------------------------------------------------------------------ */
/* Spherical collapse                                                 */
/* ------------------------------------------------------------------ */
/* state: w = y - 1 and w' (w keeps full precision while the contrast is tiny) */
static void sc_rhs(double x, const double* w, double* dw, void* p) {
  Env e;
  env_at((const Bg*)p, exp(x), &e);
  const double y = 1.0 + w[0];
  dw[0] = w[1];
  /* y^-2 - y = -(3w + 3w^2 + w^3) / y^2 */
  dw[1] = (y > 0.0) ? -(2.0 + e.dlnE) * w[1] +
                          0.5 * e.ocb * w[0] * (3.0 + w[0] * (3.0 + w[0])) / (y * y)
                    : NAN;
}

typedef struct {
  double x_c;        /* ln a at collapse (x_stop + 1 if not collapsed by x_stop) */
  double a_ta, R_ta; /* turnaround scale factor and radius (Lagrangian radius = 1) */
} Collapse;

/* (a y)' / a = y + y' = 1 + w + w' after a trial step of size h from the saved state */
static double ta_f(double h, const void* p) {
  const Ode* o = p;
  double w5[ODE_NMAX], k7[ODE_NMAX];
  if (h <= 0.0) return 1.0 + o->y[0] + o->y[1];
  ode_try(o, h, w5, k7);
  return 1.0 + w5[0] + w5[1];
}

/* evolve a shell with initial contrast eps (linear growing mode) from a_i */
static int sc_run(const Bg* b, double ai, double eps, double x_stop, Collapse* c) {
  double g[2];
  grow_init(b, ai, g);
  const double dp = eps * g[1];
  double w[2];
  w[0] = expm1(-log1p(eps) / 3.0); /* (1 + eps)^{-1/3} - 1 */
  w[1] = -(1.0 + w[0]) * dp / (3.0 * (1.0 + eps));
  Ode o;
  ode_init(&o, sc_rhs, (void*)b, 2, log(ai), w, 0.01, 1e-3 * ODE_RTOL * eps);
  int have_ta = 0;
  c->a_ta = c->R_ta = NAN;
  for (;;) {
    const Ode prev = o;
    if (ode_step(&o, x_stop)) return -1;
    if (!have_ta && 1.0 + o.y[0] + o.y[1] <= 0.0) {
      const double hs = o.x - prev.x;
      const double h = brent(ta_f, &prev, 0.0, hs, 1.0 + prev.y[0] + prev.y[1],
                             1.0 + o.y[0] + o.y[1], 1e-13 * hs);
      double w5[ODE_NMAX], k7[ODE_NMAX];
      ode_try(&prev, h, w5, k7);
      c->a_ta = exp(prev.x + h);
      c->R_ta = c->a_ta * (1.0 + w5[0]);
      have_ta = 1;
    }
    if (1.0 + o.y[0] < Y_COLLAPSE) {
      /* near collapse y ~ (x_c - x)^{2/3}, so x_c = x - (2/3) y / y' */
      c->x_c = o.x - 2.0 / 3.0 * (1.0 + o.y[0]) / o.y[1];
      return have_ta ? 0 : -1;
    }
    if (o.x >= x_stop) {
      c->x_c = x_stop + 1.0;
      return 0;
    }
  }
}

typedef struct {
  const Bg* b;
  double ai, x_target, x_stop;
} CollPar;

static double coll_f(double le, const void* p) {
  const CollPar* cp = p;
  Collapse c;
  if (sc_run(cp->b, cp->ai, exp(le), cp->x_stop, &c)) return NAN;
  return c.x_c - cp->x_target;
}

/* initial contrast ln(eps) of the shell collapsing at a_c, its linear growth d_lin(a_c)
 * (d(a_i) = 1), and the collapse history */
static int collapse_solve(const Bg* b, double a_c, double* dc, Collapse* col) {
  if (!(a_c > 0.0) || !isfinite(a_c)) return -1;
  const double ai = a_init(a_c);
  const CollPar cp = {b, ai, log(a_c), log(a_c) + 2.0};
  double dlin, dummy;
  if (growth_lin(b, GROW_CB, ai, a_c, a_c, &dlin, &dummy)) return -1;
  /* bracket around the Einstein-de Sitter guess delta_c = 1.686 */
  double lo, hi, flo, fhi;
  const double x0 = log(1.686 / dlin), f0 = coll_f(x0, &cp);
  if (!isfinite(f0)) return -1;
  if (f0 > 0.0) { /* collapses too late: increase eps */
    lo = x0, flo = f0;
    for (int i = 0;; ++i) {
      hi = lo + log(2.0);
      fhi = coll_f(hi, &cp);
      if (!isfinite(fhi) || i > 60) return -1;
      if (fhi <= 0.0) break;
      lo = hi, flo = fhi;
    }
  } else {
    hi = x0, fhi = f0;
    for (int i = 0;; ++i) {
      lo = hi - log(2.0);
      flo = coll_f(lo, &cp);
      if (!isfinite(flo) || i > 60) return -1;
      if (flo > 0.0) break;
      hi = lo, fhi = flo;
    }
  }
  const double le = brent(coll_f, &cp, lo, hi, flo, fhi, 1e-12);
  if (!isfinite(le)) return -1;
  if (sc_run(b, ai, exp(le), cp.x_stop, col)) return -1;
  *dc = exp(le) * dlin;
  return 0;
}

static double delta_crit_a(const Bg* b, double a_c) {
  double dc;
  Collapse c;
  return collapse_solve(b, a_c, &dc, &c) ? NAN : dc;
}

static double virial_Delta_a(const Bg* b, double a_c) {
  double dc;
  Collapse c;
  if (collapse_solve(b, a_c, &dc, &c)) return NAN;
  Env ec, et;
  env_at(b, a_c, &ec);
  env_at(b, c.a_ta, &et);
  const double R3 = c.R_ta * c.R_ta * c.R_ta;
  const double ac4 = a_c * a_c * a_c * a_c, at4 = c.a_ta * c.a_ta * c.a_ta * c.a_ta;
  const double eta_c = -ec.smooth / ac4 * R3 / (2.0 * b->ocb);
  const double eta_t = -et.smooth / at4 * R3 / (2.0 * b->ocb);
  /* root of 4 eta_c y^3 - 2 (1 + eta_ta) y + 1 = 0 in (0, 1) */
  if (!(4.0 * eta_c - 2.0 * (1.0 + eta_t) + 1.0 < 0.0)) return NAN;
  double lo = 0.0, hi = 1.0;
  for (int i = 0; i < 200 && hi - lo > 2.0 * DBL_EPSILON; ++i) {
    const double y = 0.5 * (lo + hi);
    if (4.0 * eta_c * y * y * y - 2.0 * (1.0 + eta_t) * y + 1.0 > 0.0) lo = y;
    else hi = y;
  }
  const double Rv = 0.5 * (lo + hi) * c.R_ta;
  double rho4, p4, mat4;
  nu_sum(b, a_c, &rho4, &p4, &mat4, NULL, NULL);
  const double rnu = mat4 / ac4; /* smooth matter-like neutrinos / rho_c0 */
  return (b->ocb / (Rv * Rv * Rv) + rnu) / (b->ocb / (ac4 / a_c) + rnu);
}

/* ------------------------------------------------------------------ */
/* Public functions                                                   */
/* ------------------------------------------------------------------ */
#define BG_OR_NAN(cp, b) \
  Bg b;                  \
  if (make_bg(cp, &b)) return NAN

void FLRW_param_init(cosm_param* cp) {
  if (!cp) return;
  memset(cp, 0, sizeof *cp);
  cp->T_cmb = FLRW_T_CMB_DEFAULT;
  cp->N_eff = FLRW_N_EFF_DEFAULT;
  cp->omega_r = FLRW_AUTO;
}

/* masses (m1, m2, m3) for the lightest mass ml */
static void nu_masses(int ordering, double ml, double* m) {
  if (ordering == FLRW_NU_NORMAL) {
    m[0] = ml;
    m[1] = sqrt(ml * ml + FLRW_DM21_SQ);
    m[2] = sqrt(ml * ml + FLRW_DM31_SQ_NORMAL);
  } else {
    m[2] = ml;
    m[1] = sqrt(ml * ml + FLRW_DM32_SQ_INVERTED);
    m[0] = sqrt(ml * ml + FLRW_DM32_SQ_INVERTED - FLRW_DM21_SQ);
  }
}

int FLRW_set_mnu_from_sum(cosm_param* cp, double sum_eV, int ordering) {
  if (!cp || !(sum_eV >= 0.0) || !isfinite(sum_eV)) return -1;
  if (ordering == FLRW_NU_DEGENERATE) {
    cp->m_nu[0] = cp->m_nu[1] = cp->m_nu[2] = sum_eV / 3.0;
    return 0;
  }
  if (ordering != FLRW_NU_NORMAL && ordering != FLRW_NU_INVERTED) return -1;
  double m[3];
  nu_masses(ordering, 0.0, m);
  if (sum_eV < m[0] + m[1] + m[2]) return -1;
  /* bisection on the lightest mass in [0, sum/3] */
  double lo = 0.0, hi = sum_eV / 3.0;
  for (int i = 0; i < 200 && hi - lo > 2.0 * DBL_EPSILON * hi; ++i) {
    const double ml = 0.5 * (lo + hi);
    nu_masses(ordering, ml, m);
    if (m[0] + m[1] + m[2] > sum_eV) hi = ml;
    else lo = ml;
  }
  nu_masses(ordering, 0.5 * (lo + hi), cp->m_nu);
  return 0;
}

double FLRW_omega_m(const cosm_param* cp) {
  BG_OR_NAN(cp, b);
  return b.om;
}

double FLRW_omega_cb(const cosm_param* cp) {
  BG_OR_NAN(cp, b);
  return b.ocb;
}

double FLRW_omega_nu(const cosm_param* cp) {
  BG_OR_NAN(cp, b);
  return b.onu;
}

double FLRW_omega_r(const cosm_param* cp) {
  BG_OR_NAN(cp, b);
  return b.og;
}

double FLRW_omega_k(const cosm_param* cp) {
  BG_OR_NAN(cp, b);
  return b.ok;
}

double FLRW_hubble_z(const cosm_param* cp, double z) {
  BG_OR_NAN(cp, b);
  return E_a(&b, a_of_z(z));
}

double FLRW_omega_m_z(const cosm_param* cp, double z) {
  BG_OR_NAN(cp, b);
  const double a = a_of_z(z);
  if (!isfinite(a)) return NAN;
  double rho4, p4, mat4;
  nu_sum(&b, a, &rho4, &p4, &mat4, NULL, NULL);
  return (b.ocb * a + mat4) / R4_a(&b, a);
}

double FLRW_hubble_time_Gyr(const cosm_param* cp) {
  /* 1/(100 km/s/Mpc) = 977.7922216807891 / 100 Gyr */
  return (cp && cp->hubble > 0.0) ? 9.777922216807891 / cp->hubble : NAN;
}

double FLRW_hubble_distance_Mpc(const cosm_param* cp) {
  return (cp && cp->hubble > 0.0) ? FLRW_HUBBLE_DISTANCE_MPC_H / cp->hubble : NAN;
}

/* ---- time ---- */
double FLRW_cosmic_time_z(const cosm_param* cp, double z) {
  BG_OR_NAN(cp, b);
  return time_a(&b, a_of_z(z));
}

double FLRW_conformal_time_z(const cosm_param* cp, double z) {
  BG_OR_NAN(cp, b);
  return conf_a(&b, a_of_z(z));
}

double FLRW_scale_factor_t(const cosm_param* cp, double t) {
  BG_OR_NAN(cp, b);
  return invert(&b, 0, t);
}

double FLRW_redshift_t(const cosm_param* cp, double t) {
  return 1.0 / FLRW_scale_factor_t(cp, t) - 1.0;
}

double FLRW_scale_factor_tau(const cosm_param* cp, double tau) {
  BG_OR_NAN(cp, b);
  return invert(&b, 1, tau);
}

double FLRW_redshift_tau(const cosm_param* cp, double tau) {
  return 1.0 / FLRW_scale_factor_tau(cp, tau) - 1.0;
}

/* ---- distances: common part for a given a ---- */
enum { DIST_C, DIST_M, DIST_A, DIST_L };

static double dist_from(const Bg* b, double a, double dc, int kind) {
  if (!isfinite(a) || !isfinite(dc)) return NAN;
  if (kind == DIST_C) return dc;
  const double dm = DM_from_DC(b, dc);
  if (kind == DIST_M) return dm;
  return (kind == DIST_A) ? dm * a : dm / a;
}

static double dist_z(const cosm_param* cp, double z, int kind) {
  BG_OR_NAN(cp, b);
  const double a = a_of_z(z);
  return isfinite(a) ? dist_from(&b, a, DC_a(&b, a), kind) : NAN;
}

static double dist_t(const cosm_param* cp, double t, int kind) {
  BG_OR_NAN(cp, b);
  const double a = invert(&b, 0, t);
  return isfinite(a) ? dist_from(&b, a, DC_a(&b, a), kind) : NAN;
}

static double dist_tau(const cosm_param* cp, double tau, int kind) {
  BG_OR_NAN(cp, b);
  const double a = invert(&b, 1, tau);
  return dist_from(&b, a, conf_a(&b, 1.0) - tau, kind); /* DC = tau0 - tau */
}

double FLRW_DC_z(const cosm_param* cp, double z) { return dist_z(cp, z, DIST_C); }
double FLRW_DM_z(const cosm_param* cp, double z) { return dist_z(cp, z, DIST_M); }
double FLRW_DA_z(const cosm_param* cp, double z) { return dist_z(cp, z, DIST_A); }
double FLRW_DL_z(const cosm_param* cp, double z) { return dist_z(cp, z, DIST_L); }

double FLRW_DC_t(const cosm_param* cp, double t) { return dist_t(cp, t, DIST_C); }
double FLRW_DM_t(const cosm_param* cp, double t) { return dist_t(cp, t, DIST_M); }
double FLRW_DA_t(const cosm_param* cp, double t) { return dist_t(cp, t, DIST_A); }
double FLRW_DL_t(const cosm_param* cp, double t) { return dist_t(cp, t, DIST_L); }

double FLRW_DC_tau(const cosm_param* cp, double tau) { return dist_tau(cp, tau, DIST_C); }
double FLRW_DM_tau(const cosm_param* cp, double tau) { return dist_tau(cp, tau, DIST_M); }
double FLRW_DA_tau(const cosm_param* cp, double tau) { return dist_tau(cp, tau, DIST_A); }
double FLRW_DL_tau(const cosm_param* cp, double tau) { return dist_tau(cp, tau, DIST_L); }

/* ---- functions of a called as functions of z, t, tau ---- */
typedef double (*AFn)(const Bg* b, double a);
typedef double (*AkFn)(const Bg* b, double k, double a);

static double call_z(const cosm_param* cp, double z, AFn f) {
  BG_OR_NAN(cp, b);
  const double a = a_of_z(z);
  return isfinite(a) ? f(&b, a) : NAN;
}

static double call_t(const cosm_param* cp, double t, AFn f) {
  BG_OR_NAN(cp, b);
  const double a = invert(&b, 0, t);
  return isfinite(a) ? f(&b, a) : NAN;
}

static double call_tau(const cosm_param* cp, double tau, AFn f) {
  BG_OR_NAN(cp, b);
  const double a = invert(&b, 1, tau);
  return isfinite(a) ? f(&b, a) : NAN;
}

static double call_k_z(const cosm_param* cp, double k, double z, AkFn f) {
  BG_OR_NAN(cp, b);
  const double a = a_of_z(z);
  return isfinite(a) ? f(&b, k, a) : NAN;
}

static double call_k_t(const cosm_param* cp, double k, double t, AkFn f) {
  BG_OR_NAN(cp, b);
  const double a = invert(&b, 0, t);
  return isfinite(a) ? f(&b, k, a) : NAN;
}

static double call_k_tau(const cosm_param* cp, double k, double tau, AkFn f) {
  BG_OR_NAN(cp, b);
  const double a = invert(&b, 1, tau);
  return isfinite(a) ? f(&b, k, a) : NAN;
}

double FLRW_growth_factor_z(const cosm_param* cp, double z) { return call_z(cp, z, growth_cb_a); }
double FLRW_growth_factor_t(const cosm_param* cp, double t) { return call_t(cp, t, growth_cb_a); }
double FLRW_growth_factor_tau(const cosm_param* cp, double tau) {
  return call_tau(cp, tau, growth_cb_a);
}

double FLRW_growth_factor_cb_z(const cosm_param* cp, double z) {
  return call_z(cp, z, growth_cb_a);
}
double FLRW_growth_factor_cb_t(const cosm_param* cp, double t) {
  return call_t(cp, t, growth_cb_a);
}
double FLRW_growth_factor_cb_tau(const cosm_param* cp, double tau) {
  return call_tau(cp, tau, growth_cb_a);
}

double FLRW_growth_factor_m_z(const cosm_param* cp, double z) { return call_z(cp, z, growth_m_a); }
double FLRW_growth_factor_m_t(const cosm_param* cp, double t) { return call_t(cp, t, growth_m_a); }
double FLRW_growth_factor_m_tau(const cosm_param* cp, double tau) {
  return call_tau(cp, tau, growth_m_a);
}

double FLRW_growth_factor_cb_k_z(const cosm_param* cp, double k, double z) {
  return call_k_z(cp, k, z, growth_cb_k_a);
}
double FLRW_growth_factor_cb_k_t(const cosm_param* cp, double k, double t) {
  return call_k_t(cp, k, t, growth_cb_k_a);
}
double FLRW_growth_factor_cb_k_tau(const cosm_param* cp, double k, double tau) {
  return call_k_tau(cp, k, tau, growth_cb_k_a);
}

double FLRW_growth_factor_m_k_z(const cosm_param* cp, double k, double z) {
  return call_k_z(cp, k, z, growth_m_k_a);
}
double FLRW_growth_factor_m_k_t(const cosm_param* cp, double k, double t) {
  return call_k_t(cp, k, t, growth_m_k_a);
}
double FLRW_growth_factor_m_k_tau(const cosm_param* cp, double k, double tau) {
  return call_k_tau(cp, k, tau, growth_m_k_a);
}

double FLRW_delta_crit_z(const cosm_param* cp, double z) { return call_z(cp, z, delta_crit_a); }
double FLRW_delta_crit_t(const cosm_param* cp, double t) { return call_t(cp, t, delta_crit_a); }
double FLRW_delta_crit_tau(const cosm_param* cp, double tau) {
  return call_tau(cp, tau, delta_crit_a);
}

double FLRW_virial_Delta_z(const cosm_param* cp, double z) {
  return call_z(cp, z, virial_Delta_a);
}
double FLRW_virial_Delta_t(const cosm_param* cp, double t) {
  return call_t(cp, t, virial_Delta_a);
}
double FLRW_virial_Delta_tau(const cosm_param* cp, double tau) {
  return call_tau(cp, tau, virial_Delta_a);
}
