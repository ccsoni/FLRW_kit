/* Tests for FLRW_kit
 *   - classic model (no radiation, no neutrinos): analytic solutions for EdS, open Lambda = 0 and
 *     flat LCDM; growth and delta_c vs direct integration; Delta_vir vs Bryan & Norman (1998)
 *   - parameters: defaults, photon density, neutrino masses from the sum and the ordering,
 *     Onu h^2 vs sum(m)
 *   - radiation + matter: analytic t(a), tau(a) and the Meszaros growth solution
 *   - massive neutrinos: E(z) and the growth limits vs independent Fermi-Dirac integrals,
 *     k -> 0 and k -> infinity limits of the scale-dependent growth, small-scale suppression
 *   - delta_c with radiation vs direct integration of the background and shell equations
 *   - round trips through the inverse functions, z / t / tau consistency, NaN for invalid input
 */
#include <math.h>
#include <stdio.h>

#include "flrw_kit.h"

#define PI 3.14159265358979323846
#define OG_H2 2.472975333e-5 /* Og h^2 at T_cmb = 2.7255 K */
#define KB_EV 8.617333262e-5

static int nfail = 0, ntest = 0;

static void check(const char* name, double got, double ref, double rtol) {
  const double rel = fabs(got - ref) / fmax(fabs(ref), 1e-300);
  const int ok = (rel <= rtol);
  ++ntest;
  if (!ok) ++nfail;
  printf("%-4s %-46s got=% .12e ref=% .12e rel=%.1e\n", ok ? "ok" : "FAIL", name, got, ref, rel);
}

static void check_true(const char* name, int ok) {
  ++ntest;
  if (!ok) ++nfail;
  printf("%-4s %s\n", ok ? "ok" : "FAIL", name);
}

/* ---- independent model for the checks ---- */
typedef struct {
  double ocb, ok, ol, og; /* Ocb, Ok, OL, Og */
  double onu1;            /* one massless species */
  int nm;
  double nless, y0[3];    /* massive species: m / (k_B T_nu0) */
} Model;

/* F(y) or P(y) (P: pressure) by composite Simpson on [0, 60] */
static double fd_int(double y, int pressure, int n) {
  const double h = 60.0 / n;
  double s = 0.0;
  for (int i = 0; i <= n; ++i) {
    const double q = i * h, E = sqrt(q * q + y * y);
    const double f = (E > 0.0) ? q * q * (pressure ? q * q / (3.0 * E) : E) / (exp(q) + 1.0) : 0.0;
    s += f * ((i == 0 || i == n) ? 1.0 : ((i & 1) ? 4.0 : 2.0));
  }
  return s * h / 3.0 * 120.0 / (7.0 * PI * PI * PI * PI);
}

static void model_from(const cosm_param* cp, Model* m, int n) {
  const double h = cp->hubble;
  m->ocb = cp->omega_cdm + cp->omega_b;
  m->ol = cp->omega_de;
  const double ogT = OG_H2 * pow(cp->T_cmb / 2.7255, 4.0) / (h * h);
  m->og = (cp->omega_r >= 0.0) ? cp->omega_r : ogT;
  m->onu1 = ogT * 0.875 * pow(4.0 / 11.0, 4.0 / 3.0) * cp->N_eff / 3.0;
  const double kT = KB_EV * cbrt(4.0 / 11.0) * pow(cp->N_eff / 3.0, 0.25) * cp->T_cmb;
  m->nm = 0;
  m->nless = 0.0;
  for (int i = 0; i < 3; ++i) {
    if (cp->m_nu[i] == 0.0) m->nless += 1.0;
    else m->y0[m->nm++] = cp->m_nu[i] / kT;
  }
  double onu = m->onu1 * m->nless;
  for (int i = 0; i < m->nm; ++i) onu += m->onu1 * fd_int(m->y0[i], 0, n);
  m->ok = 1.0 - m->ocb - m->og - onu - m->ol;
}

/* a^4 rho_nu and a^4 p_nu (units of rho_c0) */
static void model_nu(const Model* m, double a, double* r4, double* p4, int n) {
  *r4 = m->onu1 * m->nless;
  *p4 = m->onu1 * m->nless / 3.0;
  for (int i = 0; i < m->nm; ++i) {
    *r4 += m->onu1 * fd_int(m->y0[i] * a, 0, n);
    *p4 += m->onu1 * fd_int(m->y0[i] * a, 1, n);
  }
}

/* ---- growth equation, RK4 in x = ln a, with the rhs evaluated from tabulated nu densities ---
 * mode 0: small-scale limit (source Ocb), 1: large-scale limit (source Ocb + rho - 3p of nu) */
typedef struct {
  double R, dlnE, src;
} Pt;

static void grow_point(const Model* m, double a, int mode, int n, Pt* pt) {
  double r4, p4;
  model_nu(m, a, &r4, &p4, n);
  const double a2 = a * a;
  pt->R = m->og + r4 + m->ocb * a + a2 * (m->ok + m->ol * a2);
  pt->dlnE = (-4.0 * m->og - 3.0 * (r4 + p4) - 3.0 * m->ocb * a - 2.0 * m->ok * a2) / (2.0 * pt->R);
  pt->src = (m->ocb * a + (mode ? r4 - 3.0 * p4 : 0.0)) / pt->R; /* massless: r - 3p = 0 */
}

static void grow_seg(const Model* m, int mode, int n, double x0, double x1, int ns, double* y) {
  const double h = (x1 - x0) / ns;
  Pt p0, pm, p1;
  grow_point(m, exp(x0), mode, n, &p0);
  for (int s = 0; s < ns; ++s) {
    const double x = x0 + s * h;
    grow_point(m, exp(x + 0.5 * h), mode, n, &pm);
    grow_point(m, exp(x + h), mode, n, &p1);
    double k[4][2], t[2];
    const Pt* P[4] = {&p0, &pm, &pm, &p1};
    const double c[4] = {0.0, 0.5, 0.5, 1.0};
    for (int j = 0; j < 4; ++j) {
      for (int i = 0; i < 2; ++i) t[i] = y[i] + (j ? c[j] * h * k[j - 1][i] : 0.0);
      k[j][0] = t[1];
      k[j][1] = -(2.0 + P[j]->dlnE) * t[1] + 1.5 * P[j]->src * t[0];
    }
    for (int i = 0; i < 2; ++i)
      y[i] += h / 6.0 * (k[0][i] + 2.0 * k[1][i] + 2.0 * k[2][i] + k[3][i]);
    p0 = p1;
  }
}

/* D(a_end)/D(1) for a_end <= 1, Meszaros initial conditions at a_i = 1e-5 */
static double growth_ode(const Model* m, int mode, double a_end, int nsteps, int n) {
  const double ai = 1e-5, aeq = (m->og + m->onu1 * (m->nless + m->nm)) / m->ocb;
  const double xi = log(ai), xe = log(a_end);
  double y[2] = {ai + 2.0 / 3.0 * aeq, ai};
  const int n1 = (int)ceil(nsteps * (xe - xi) / (-xi));
  grow_seg(m, mode, n, xi, xe, n1, y);
  const double d_end = y[0];
  if (xe < 0.0) grow_seg(m, mode, n, xe, 0.0, nsteps > n1 ? nsteps - n1 : 1, y);
  return d_end / y[0];
}

/* ---- spherical collapse by direct integration in time (smooth radiation and OL) ----
 *   a'' = -Om/(2a^2) - Og a^-3 + OL a,   R'' = -Om/(2R^2) - Og a^-4 R + OL R   (r_L = 1),
 *   written for w = R/a - 1 (y = 1 + w):  y'' = -Om/(2a^3) (y^-2 - y) - 2 (a'/a) y' */
static void sc_rhs(const Model* m, const double y[4], double dy[4]) {
  const double a = y[0], w = y[2], yy = 1.0 + w;
  dy[0] = y[1];
  dy[1] = -0.5 * m->ocb / (a * a) - m->og / (a * a * a) + m->ol * a;
  dy[2] = y[3];
  dy[3] = 0.5 * m->ocb / (a * a * a) * w * (3.0 + w * (3.0 + w)) / (yy * yy) -
          2.0 * y[1] / a * y[3];
}

static double collapse_a(const Model* m, double ai, double eps) {
  const double aeq = m->og / m->ocb;
  const double adot = sqrt(m->ocb / ai + m->og / (ai * ai) + m->ok + m->ol * ai * ai);
  const double ddot = (adot / ai) * eps * ai / (ai + 2.0 / 3.0 * aeq); /* dd/dt = H dd/dlna */
  const double w = expm1(-log1p(eps) / 3.0);
  double y[4] = {ai, adot, w, -(1.0 + w) * ddot / (3.0 * (1.0 + eps))};
  double k1[4], k2[4], k3[4], k4[4], t[4];
  while (1.0 + y[2] > 1e-6) {
    if (y[0] > 1e3 || y[1] <= 0.0) return INFINITY; /* not collapsed while expanding */
    double h = y[0] / y[1];
    if (y[3] != 0.0) h = fmin(h, (1.0 + y[2]) / fabs(y[3]));
    h *= 2e-3;
    sc_rhs(m, y, k1);
    for (int j = 0; j < 4; ++j) t[j] = y[j] + 0.5 * h * k1[j];
    sc_rhs(m, t, k2);
    for (int j = 0; j < 4; ++j) t[j] = y[j] + 0.5 * h * k2[j];
    sc_rhs(m, t, k3);
    for (int j = 0; j < 4; ++j) t[j] = y[j] + h * k3[j];
    sc_rhs(m, t, k4);
    for (int j = 0; j < 4; ++j) y[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
  }
  return y[0];
}

/* the same collapse in x = ln a for w = y - 1 (RK4, step control near collapse):
 *   w'' = -(2 + dlnE/dx) w' + (1/2) Ocb(a) w (3 + 3w + w^2) / (1 + w)^2 */
static void scx_rhs(const Model* m, double x, const double w[2], double dw[2]) {
  const double a = exp(x), a2 = a * a;
  const double R = m->og + m->ocb * a + a2 * (m->ok + m->ol * a2);
  const double dlnE = (-4.0 * m->og - 3.0 * m->ocb * a - 2.0 * m->ok * a2) / (2.0 * R);
  const double y = 1.0 + w[0];
  dw[0] = w[1];
  dw[1] = -(2.0 + dlnE) * w[1] +
          0.5 * m->ocb * a / R * w[0] * (3.0 + w[0] * (3.0 + w[0])) / (y * y);
}

static double collapse_a_x(const Model* m, double ai, double eps) {
  const double aeq = m->og / m->ocb;
  double x = log(ai), w[2], k1[2], k2[2], k3[2], k4[2], t[2];
  w[0] = expm1(-log1p(eps) / 3.0);
  w[1] = -(1.0 + w[0]) * eps * ai / (ai + 2.0 / 3.0 * aeq) / (3.0 * (1.0 + eps));
  while (1.0 + w[0] > 1e-6) {
    if (x > log(1e3)) return INFINITY;
    double h = 1e-3;
    if (w[1] != 0.0) h = fmin(h, 1e-3 * (1.0 + w[0]) / fabs(w[1]));
    scx_rhs(m, x, w, k1);
    for (int j = 0; j < 2; ++j) t[j] = w[j] + 0.5 * h * k1[j];
    scx_rhs(m, x + 0.5 * h, t, k2);
    for (int j = 0; j < 2; ++j) t[j] = w[j] + 0.5 * h * k2[j];
    scx_rhs(m, x + 0.5 * h, t, k3);
    for (int j = 0; j < 2; ++j) t[j] = w[j] + h * k3[j];
    scx_rhs(m, x + h, t, k4);
    for (int j = 0; j < 2; ++j) w[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
    x += h;
  }
  return exp(x);
}

/* initial contrast collapsing at a_c (bisection), extrapolated with the linear growth;
 * in_x selects the integration in ln a (used with radiation, where the integration in time
 * is too sensitive to round-off in the radiation era) */
static double delta_crit_ode(const Model* m, double a_c, int in_x) {
  const double ai = 1e-7, aeq = m->og / m->ocb, d0 = ai + 2.0 / 3.0 * aeq;
  double lo = 1e-2 * d0, hi = 100.0 * d0;
  for (int it = 0; it < 52; ++it) {
    const double mid = sqrt(lo * hi);
    if ((in_x ? collapse_a_x(m, ai, mid) : collapse_a(m, ai, mid)) > a_c) lo = mid;
    else hi = mid;
  }
  /* linear growth from the same initial conditions: D(a_c) / D(a_i) */
  const int ns = 40000;
  const double xi = log(ai), h = (log(a_c) - xi) / ns;
  double y[2] = {1.0, ai / d0};
  for (int s = 0; s < ns; ++s) {
    double k[4][2], tt[2];
    const double c[4] = {0.0, 0.5, 0.5, 1.0};
    for (int j = 0; j < 4; ++j) {
      const double a = exp(xi + (s + c[j]) * h), a2 = a * a;
      const double R = m->og + m->ocb * a + a2 * (m->ok + m->ol * a2);
      const double dlnE = (-4.0 * m->og - 3.0 * m->ocb * a - 2.0 * m->ok * a2) / (2.0 * R);
      for (int i = 0; i < 2; ++i) tt[i] = y[i] + (j ? c[j] * h * k[j - 1][i] : 0.0);
      k[j][0] = tt[1];
      k[j][1] = -(2.0 + dlnE) * tt[1] + 1.5 * m->ocb * a / R * tt[0];
    }
    for (int i = 0; i < 2; ++i)
      y[i] += h / 6.0 * (k[0][i] + 2.0 * k[1][i] + 2.0 * k[2][i] + k[3][i]);
  }
  return sqrt(lo * hi) * y[0];
}

/* classic model: no radiation, no neutrinos */
static void set_classic(cosm_param* cp, double om, double ol) {
  FLRW_param_init(cp);
  cp->omega_cdm = 0.85 * om;
  cp->omega_b = 0.15 * om;
  cp->omega_de = ol;
  cp->hubble = 0.7;
  cp->N_eff = 0.0;
  cp->omega_r = 0.0;
}

/* realistic model with radiation and neutrinos */
static void set_full(cosm_param* cp, double ocb, double ol, double sum_mnu, int ordering) {
  FLRW_param_init(cp);
  cp->omega_cdm = 0.84 * ocb;
  cp->omega_b = 0.16 * ocb;
  cp->omega_de = ol;
  cp->hubble = 0.677;
  FLRW_set_mnu_from_sum(cp, sum_mnu, ordering);
}

int main(void) {
  cosm_param cp;
  char name[128];

  /* ================= classic model ================= */
  set_classic(&cp, 1.0, 0.0);
  puts("== classic: Einstein-de Sitter");
  for (double z = 0.0; z <= 10.0; z += 5.0) {
    const double a = 1.0 / (1.0 + z);
    snprintf(name, sizeof name, "t(z=%g)", z);
    check(name, FLRW_cosmic_time_z(&cp, z), 2.0 / 3.0 * pow(a, 1.5), 1e-12);
    snprintf(name, sizeof name, "tau(z=%g)", z);
    check(name, FLRW_conformal_time_z(&cp, z), 2.0 * sqrt(a), 1e-12);
    snprintf(name, sizeof name, "DL(z=%g)", z);
    check(name, FLRW_DL_z(&cp, z), 2.0 * (1.0 + z) * (1.0 - sqrt(a)), 1e-12);
    snprintf(name, sizeof name, "D(z=%g)", z);
    check(name, FLRW_growth_factor_z(&cp, z), a, 1e-10);
    snprintf(name, sizeof name, "delta_c(z=%g)", z);
    check(name, FLRW_delta_crit_z(&cp, z), 0.6 * pow(1.5 * PI, 2.0 / 3.0), 1e-8);
    snprintf(name, sizeof name, "Delta_vir(z=%g)", z);
    check(name, FLRW_virial_Delta_z(&cp, z), 18.0 * PI * PI, 1e-8);
  }

  set_classic(&cp, 0.3, 0.0);
  puts("== classic: open, Om=0.3, OL=0");
  for (double z = 0.5; z <= 8.0; z *= 4.0) {
    const double om = 0.3;
    const double dl = 2.0 * (om * z + (om - 2.0) * (sqrt(1.0 + om * z) - 1.0)) / (om * om);
    snprintf(name, sizeof name, "DL(z=%g) Mattig", z);
    check(name, FLRW_DL_z(&cp, z), dl, 1e-11);
  }
  {
    const double om = 0.3, w = 1.0 - om;
    check("t0", FLRW_cosmic_time_z(&cp, 0.0),
          1.0 / w - om / (2.0 * pow(w, 1.5)) * acosh(2.0 / om - 1.0), 1e-12);
  }

  set_classic(&cp, 0.3, 0.7);
  puts("== classic: flat LCDM, Om=0.3, OL=0.7");
  for (double z = 0.0; z <= 9.0; z += 3.0) {
    const double a = 1.0 / (1.0 + z);
    snprintf(name, sizeof name, "t(z=%g)", z);
    check(name, FLRW_cosmic_time_z(&cp, z),
          2.0 / (3.0 * sqrt(0.7)) * asinh(sqrt(0.7 / 0.3) * pow(a, 1.5)), 1e-12);
  }
  /* values of the closed-form spherical collapse solution (FLRW_kit 1.0) */
  check("delta_c(z=0) vs closed form", FLRW_delta_crit_z(&cp, 0.0), 1.675526794754, 1e-8);
  check("Delta_vir(z=0) vs closed form", FLRW_virial_Delta_z(&cp, 0.0), 339.785156245596, 1e-8);

  {
    const double models[][2] = {{0.3, 0.7}, {0.3, 0.0}, {0.3, 0.5}, {2.0, 0.0}};
    for (size_t k = 0; k < sizeof models / sizeof models[0]; ++k) {
      set_classic(&cp, models[k][0], models[k][1]);
      Model m;
      model_from(&cp, &m, 6000);
      printf("== classic: ODE comparison, Om=%g, OL=%g\n", models[k][0], models[k][1]);
      for (double z = 0.0; z <= 4.0; z += 2.0) {
        const double a = 1.0 / (1.0 + z);
        snprintf(name, sizeof name, "D(z=%g) vs ODE", z);
        check(name, FLRW_growth_factor_z(&cp, z), growth_ode(&m, 0, a, 4000, 2), 1e-8);
        snprintf(name, sizeof name, "delta_c(z=%g) vs ODE", z);
        check(name, FLRW_delta_crit_z(&cp, z), delta_crit_ode(&m, a, 0), 5e-5);
      }
    }
  }

  puts("== classic: Delta_vir vs Bryan & Norman (~1%)");
  set_classic(&cp, 0.3, 0.7);
  for (double z = 0.0; z <= 3.0; z += 1.5) {
    const double omz = FLRW_omega_m_z(&cp, z), x = omz - 1.0;
    snprintf(name, sizeof name, "flat Delta_vir(z=%g)", z);
    check(name, FLRW_virial_Delta_z(&cp, z), (18.0 * PI * PI + 82.0 * x - 39.0 * x * x) / omz,
          1.5e-2);
  }
  set_classic(&cp, 0.3, 0.0);
  for (double z = 0.0; z <= 3.0; z += 1.5) {
    const double omz = FLRW_omega_m_z(&cp, z), x = omz - 1.0;
    snprintf(name, sizeof name, "open Delta_vir(z=%g)", z);
    check(name, FLRW_virial_Delta_z(&cp, z), (18.0 * PI * PI + 60.0 * x - 32.0 * x * x) / omz,
          1.5e-2);
  }

  /* ================= parameters ================= */
  puts("== parameters");
  {
    cosm_param a, b = {FLRW_PARAM_DEFAULTS, .omega_cdm = 0.26, .hubble = 0.7};
    FLRW_param_init(&a);
    a.omega_cdm = 0.26;
    a.hubble = 0.7;
    check_true("FLRW_param_init == FLRW_PARAM_DEFAULTS",
               a.T_cmb == b.T_cmb && a.N_eff == b.N_eff && a.omega_r == b.omega_r &&
                   a.m_nu[0] == 0.0 && b.m_nu[2] == 0.0 && b.omega_de == 0.0);
    check("default T_cmb", a.T_cmb, 2.7255, 0.0);
    check("default N_eff", a.N_eff, 3.044, 0.0);
    check("Og auto = 2.472975e-5 / h^2", FLRW_omega_r(&a), OG_H2 / 0.49, 1e-15);
    a.T_cmb = 2.0;
    check("Og scales as T^4", FLRW_omega_r(&a), OG_H2 / 0.49 * pow(2.0 / 2.7255, 4.0), 1e-14);
    a.omega_r = 1e-4;
    check("explicit omega_r", FLRW_omega_r(&a), 1e-4, 0.0);
    a.omega_r = FLRW_AUTO;
    a.T_cmb = 2.7255;
    check("massless nu: Onu = 0.6913 N_eff/3 Og",
          FLRW_omega_nu(&a), FLRW_omega_r(&a) * 0.875 * pow(4.0 / 11.0, 4.0 / 3.0) * 3.044, 1e-14);
  }
  {
    double m[3];
    for (int ord = 0; ord < 2; ++ord) {
      for (double sum = (ord ? 0.1 : 0.06); sum < 0.5; sum *= 3.0) {
        FLRW_param_init(&cp);
        const int ret = FLRW_set_mnu_from_sum(&cp, sum, ord);
        for (int i = 0; i < 3; ++i) m[i] = cp.m_nu[i];
        snprintf(name, sizeof name, "%s sum=%g: sum of masses", ord ? "IO" : "NO", sum);
        check(name, m[0] + m[1] + m[2], sum, 1e-14);
        snprintf(name, sizeof name, "%s sum=%g: m2^2 - m1^2", ord ? "IO" : "NO", sum);
        check(name, m[1] * m[1] - m[0] * m[0], 7.42e-5, 1e-9);
        snprintf(name, sizeof name, "%s sum=%g: |m3^2 - m%d^2|", ord ? "IO" : "NO", sum,
                 ord ? 2 : 1);
        check(name, ord ? m[1] * m[1] - m[2] * m[2] : m[2] * m[2] - m[0] * m[0],
              ord ? 2.498e-3 : 2.515e-3, 1e-9);
        check_true("  returns 0", ret == 0);
      }
    }
    FLRW_param_init(&cp);
    check_true("NO sum = 0.05 eV rejected", FLRW_set_mnu_from_sum(&cp, 0.05, FLRW_NU_NORMAL) == -1);
    check_true("IO sum = 0.09 eV rejected",
               FLRW_set_mnu_from_sum(&cp, 0.09, FLRW_NU_INVERTED) == -1);
    check_true("unknown ordering rejected", FLRW_set_mnu_from_sum(&cp, 0.3, 7) == -1);
    FLRW_set_mnu_from_sum(&cp, 0.3, FLRW_NU_DEGENERATE);
    check_true("degenerate: m = sum/3",
               cp.m_nu[0] == 0.3 / 3.0 && cp.m_nu[1] == 0.3 / 3.0 && cp.m_nu[2] == 0.3 / 3.0);
    cp.hubble = 0.7;
    cp.omega_cdm = 0.3;
    const double onu_nr = FLRW_omega_nu(&cp) * 0.49;
    check("Onu h^2 = sum m / 93.04 eV (N_eff = 3.044)", onu_nr, 0.3 / 93.0425, 1e-4);
    cp.N_eff = 3.0;
    check("Onu h^2 = sum m / 94.06 eV (N_eff = 3)", FLRW_omega_nu(&cp) * 0.49, 0.3 / 94.0641, 1e-4);
  }
  set_full(&cp, 0.3, 0.7, 0.06, FLRW_NU_NORMAL);
  check("Ok = 1 - Ocb - Og - Onu - OL", FLRW_omega_k(&cp),
        1.0 - FLRW_omega_cb(&cp) - FLRW_omega_r(&cp) - FLRW_omega_nu(&cp) - 0.7, 1e-12);

  /* ================= radiation + matter ================= */
  puts("== radiation + matter (flat, OL = 0, no neutrinos)");
  {
    FLRW_param_init(&cp);
    const double orad = 1e-4, om = 1.0 - orad, aeq = orad / om;
    cp.omega_cdm = om;
    cp.hubble = 0.7;
    cp.N_eff = 0.0;
    cp.omega_r = orad;
    for (double z = 0.0; z <= 1e5; z = (z == 0.0) ? 10.0 : z * 10.0) {
      const double a = 1.0 / (1.0 + z), s = om * a + orad;
      snprintf(name, sizeof name, "t(z=%g)", z);
      check(name, FLRW_cosmic_time_z(&cp, z),
            2.0 / (3.0 * om * om) * (sqrt(s) * (s - 3.0 * orad) + 2.0 * pow(orad, 1.5)), 1e-11);
      snprintf(name, sizeof name, "tau(z=%g)", z);
      check(name, FLRW_conformal_time_z(&cp, z), 2.0 / om * (sqrt(s) - sqrt(orad)), 1e-11);
      if (z > 1e3) continue;
      snprintf(name, sizeof name, "D(z=%g) Meszaros", z);
      check(name, FLRW_growth_factor_z(&cp, z), (a / aeq + 2.0 / 3.0) / (1.0 / aeq + 2.0 / 3.0),
            1e-9);
    }
    Model m;
    model_from(&cp, &m, 6000);
    for (double z = 0.0; z <= 3.0; z += 3.0) {
      snprintf(name, sizeof name, "delta_c(z=%g) vs ODE", z);
      check(name, FLRW_delta_crit_z(&cp, z), delta_crit_ode(&m, 1.0 / (1.0 + z), 1), 1e-5);
    }
    cp.omega_cdm = 0.3 - orad;
    cp.omega_de = 0.7;
    model_from(&cp, &m, 6000);
    check("LCDM + radiation: delta_c(z=0) vs ODE", FLRW_delta_crit_z(&cp, 0.0),
          delta_crit_ode(&m, 1.0, 1), 1e-5);
  }

  /* ================= massive neutrinos ================= */
  {
    FLRW_param_init(&cp);
    cp.omega_cdm = 0.26;
    cp.omega_b = 0.05;
    cp.omega_de = 0.69;
    cp.hubble = 0.68;
    cp.m_nu[0] = 0.001; /* small-y series */
    cp.m_nu[1] = 0.05;  /* table */
    cp.m_nu[2] = 3.0;   /* large-y expansion today */
    Model m;
    model_from(&cp, &m, 200000);
    puts("== massive neutrinos: E(z) vs direct Fermi-Dirac integrals");
    for (double z = 0.0; z <= 1e6; z = (z == 0.0) ? 1.0 : z * 10.0) {
      const double a = 1.0 / (1.0 + z);
      double r4, p4;
      model_nu(&m, a, &r4, &p4, 200000);
      const double E =
          sqrt((m.og + r4) / pow(a, 4.0) + m.ocb / pow(a, 3.0) + m.ok / (a * a) + m.ol);
      snprintf(name, sizeof name, "E(z=%g)", z);
      check(name, FLRW_hubble_z(&cp, z), E, 1e-11);
    }
  }
  {
    set_full(&cp, 0.31, 0.69, 0.3, FLRW_NU_DEGENERATE);
    Model m;
    model_from(&cp, &m, 6000);
    puts("== massive neutrinos (sum = 0.3 eV): growth limits vs ODE");
    for (double z = 0.0; z <= 4.0; z += 2.0) {
      const double a = 1.0 / (1.0 + z);
      snprintf(name, sizeof name, "D_cb(z=%g) vs ODE", z);
      check(name, FLRW_growth_factor_cb_z(&cp, z), growth_ode(&m, 0, a, 1500, 6000), 1e-6);
      snprintf(name, sizeof name, "D_m(z=%g) vs ODE", z);
      check(name, FLRW_growth_factor_m_z(&cp, z), growth_ode(&m, 1, a, 1500, 6000), 1e-6);
    }
    check("growth_factor == growth_factor_cb", FLRW_growth_factor_z(&cp, 1.0),
          FLRW_growth_factor_cb_z(&cp, 1.0), 0.0);
  }
  {
    puts("== massive neutrinos: limits of the scale-dependent growth");
    const double sums[] = {0.06, 0.3};
    for (int s = 0; s < 2; ++s) {
      set_full(&cp, 0.31, 0.69, sums[s], s ? FLRW_NU_DEGENERATE : FLRW_NU_NORMAL);
      const double fnu = FLRW_omega_nu(&cp) / FLRW_omega_m(&cp);
      for (double z = 0.0; z <= 3.0; z += 3.0) {
        const double Dm = FLRW_growth_factor_m_z(&cp, z), Dcb = FLRW_growth_factor_cb_z(&cp, z);
        snprintf(name, sizeof name, "sum=%g: D_cb_k(k=0, z=%g) = D_m", sums[s], z);
        check(name, FLRW_growth_factor_cb_k_z(&cp, 0.0, z), Dm, 1e-9);
        snprintf(name, sizeof name, "sum=%g: D_m_k(k=1e-3, z=%g) = D_m", sums[s], z);
        check(name, FLRW_growth_factor_m_k_z(&cp, 1e-3, z), Dm, 1e-8);
        snprintf(name, sizeof name, "sum=%g: D_cb_k(1e7,z=%g)/D_cb_k(1e7,0) = D_cb", sums[s], z);
        check(name,
              FLRW_growth_factor_cb_k_z(&cp, 1e7, z) / FLRW_growth_factor_cb_k_z(&cp, 1e7, 0.0),
              Dcb, 1e-7);
      }
      double prev = 2.0;
      int mono = 1;
      for (double k = 1.0; k <= 1e5; k *= 10.0) {
        const double d = FLRW_growth_factor_m_k_z(&cp, k, 0.0);
        mono = mono && (d < prev);
        prev = d;
      }
      snprintf(name, sizeof name, "sum=%g: D_m_k(z=0) decreases with k", sums[s]);
      check_true(name, mono);
      const double sup = 1.0 - prev * prev;
      snprintf(name, sizeof name, "sum=%g: 1 - D_m_k^2 = %.2f f_nu (6..10)", sums[s], sup / fnu);
      check_true(name, sup > 6.0 * fnu && sup < 10.0 * fnu);
    }
    /* massless neutrinos only: no scale dependence */
    set_full(&cp, 0.31, 0.69, 0.0, FLRW_NU_DEGENERATE);
    check("massless nu: D_cb_k(1e4) = D_cb", FLRW_growth_factor_cb_k_z(&cp, 1e4, 1.0),
          FLRW_growth_factor_cb_z(&cp, 1.0), 1e-9);
    check("massless nu: D_m = D_cb", FLRW_growth_factor_m_z(&cp, 1.0),
          FLRW_growth_factor_cb_z(&cp, 1.0), 1e-12);
  }
  {
    puts("== massive neutrinos: small-scale growth exponent in the matter era");
    /* Ocb + Onu = 1, OL = 0: D_cb ~ a^p, p = (sqrt(25 - 24 f) - 1) / 4 */
    FLRW_param_init(&cp);
    cp.hubble = 0.7;
    FLRW_set_mnu_from_sum(&cp, 3.0, FLRW_NU_DEGENERATE);
    cp.omega_cdm = 1.0;
    const double onu = FLRW_omega_nu(&cp);
    cp.omega_cdm = 1.0 - onu - FLRW_omega_r(&cp);
    const double f = onu / (cp.omega_cdm + onu);
    const double p = (sqrt(25.0 - 24.0 * f) - 1.0) / 4.0;
    const double slope =
        log(FLRW_growth_factor_cb_z(&cp, 0.0) / FLRW_growth_factor_cb_z(&cp, 9.0)) /
                         log(10.0);
    snprintf(name, sizeof name, "d lnD_cb / d lna (f_nu = %.3f)", f);
    check(name, slope, p, 2e-3);
  }
  {
    puts("== massive neutrinos: spherical collapse");
    set_full(&cp, 0.31, 0.69, 0.0, FLRW_NU_DEGENERATE);
    const double dc0 = FLRW_delta_crit_z(&cp, 0.0), dv0 = FLRW_virial_Delta_z(&cp, 0.0);
    set_full(&cp, 0.31 - 0.3 / 93.0425 / 0.677 / 0.677, 0.69, 0.3, FLRW_NU_DEGENERATE);
    check("delta_c(0): sum = 0.3 eV vs massless", FLRW_delta_crit_z(&cp, 0.0), dc0, 2e-3);
    check("Delta_vir(0): sum = 0.3 eV vs massless", FLRW_virial_Delta_z(&cp, 0.0), dv0, 3e-2);
  }

  /* ================= inverse functions and t, tau versions ================= */
  {
    cosm_param models[3];
    set_classic(&models[0], 0.3, 0.0);
    set_full(&models[1], 0.31, 0.69, 0.06, FLRW_NU_NORMAL);
    set_full(&models[2], 0.4, 0.5, 0.15, FLRW_NU_INVERTED);
    for (int k = 0; k < 3; ++k) {
      cp = models[k];
      printf("== round trip, model %d\n", k);
      const double zs[] = {-0.5, 0.0, 2.0, 30.0, 1e4, 1e7};
      for (size_t i = 0; i < sizeof zs / sizeof zs[0]; ++i) {
        const double z = zs[i], t = FLRW_cosmic_time_z(&cp, z), tau = FLRW_conformal_time_z(&cp, z);
        snprintf(name, sizeof name, "a(t(z=%g))", z);
        check(name, FLRW_scale_factor_t(&cp, t), 1.0 / (1.0 + z), 1e-12);
        snprintf(name, sizeof name, "a(tau(z=%g))", z);
        check(name, FLRW_scale_factor_tau(&cp, tau), 1.0 / (1.0 + z), 1e-12);
        if (z < 0.0 || z > 100.0) continue;
        snprintf(name, sizeof name, "DL_t vs DL_z (z=%g)", z);
        check(name, FLRW_DL_t(&cp, t), FLRW_DL_z(&cp, z), 1e-10);
        snprintf(name, sizeof name, "DA_tau vs DA_z (z=%g)", z);
        check(name, FLRW_DA_tau(&cp, tau), FLRW_DA_z(&cp, z), 1e-10);
        snprintf(name, sizeof name, "D_m_tau vs D_m_z (z=%g)", z);
        check(name, FLRW_growth_factor_m_tau(&cp, tau), FLRW_growth_factor_m_z(&cp, z), 1e-9);
        snprintf(name, sizeof name, "D_cb_k_t vs D_cb_k_z (z=%g)", z);
        check(name, FLRW_growth_factor_cb_k_t(&cp, 300.0, t),
              FLRW_growth_factor_cb_k_z(&cp, 300.0, z), 1e-8);
        snprintf(name, sizeof name, "delta_c_t vs delta_c_z (z=%g)", z);
        check(name, FLRW_delta_crit_t(&cp, t), FLRW_delta_crit_z(&cp, z), 1e-8);
        snprintf(name, sizeof name, "Delta_vir_tau vs _z (z=%g)", z);
        check(name, FLRW_virial_Delta_tau(&cp, tau), FLRW_virial_Delta_z(&cp, z), 1e-8);
      }
    }
  }

  /* ================= invalid input ================= */
  puts("== invalid input (expect NaN)");
  {
    set_full(&cp, 0.31, 0.69, 0.06, FLRW_NU_NORMAL);
    const double tau_inf = 2.0 * FLRW_conformal_time_z(&cp, -1.0 + 1e-12);
    cosm_param bad_h = cp, bad_m = cp, bad_cb = cp, bad_T = cp;
    bad_h.hubble = 0.0;
    bad_m.m_nu[1] = -0.1;
    bad_cb.omega_cdm = bad_cb.omega_b = 0.0;
    bad_T.T_cmb = 0.0;
    const double v[] = {FLRW_cosmic_time_z(&cp, -1.5),
                        FLRW_scale_factor_t(&cp, -1.0),
                        FLRW_scale_factor_tau(&cp, tau_inf),
                        FLRW_growth_factor_cb_k_z(&cp, -1.0, 0.0),
                        FLRW_hubble_z(&bad_h, 0.0),
                        FLRW_hubble_z(&bad_m, 0.0),
                        FLRW_hubble_z(&bad_cb, 0.0),
                        FLRW_hubble_z(&bad_T, 0.0)};
    for (int i = 0; i < (int)(sizeof v / sizeof v[0]); ++i) {
      snprintf(name, sizeof name, "NaN check %d (%g)", i, v[i]);
      check_true(name, isnan(v[i]));
    }
  }

  printf("\n%d / %d tests passed\n", ntest - nfail, ntest);
  return nfail != 0;
}
