/* Tests for FLRW_kit
 *   - analytic solutions for EdS, open Lambda = 0 and flat LCDM
 *   - linear growth factor vs direct RK4 integration of the growth equation
 *   - delta_c vs spherical collapse from direct integration of the background and shell equations
 *   - Delta_vir vs the Bryan & Norman (1998) fit (~1%)
 *   - round trips through the inverse functions (z -> t -> z, z -> tau -> z)
 */
#include <math.h>
#include <stdio.h>

#include "flrw_kit.h"

#define PI 3.14159265358979323846

static int nfail = 0, ntest = 0;

static void check(const char* name, double got, double ref, double rtol) {
  const double rel = fabs(got - ref) / fmax(fabs(ref), 1e-300);
  const int ok = (rel <= rtol);
  ++ntest;
  if (!ok) ++nfail;
  printf("%-4s %-44s got=% .12e ref=% .12e rel=%.1e\n", ok ? "ok" : "FAIL", name, got, ref, rel);
}

/* ---- RK4 integration in ln a of D'' + (2 + dlnE/dlna) D' = (3/2) Om(a) D ---- */
typedef struct {
  double om, ok, ol;
} Model;

static void grow_rhs(const Model* m, double x, const double y[2], double dy[2]) {
  const double a = exp(x);
  const double E2 = m->om / (a * a * a) + m->ok / (a * a) + m->ol;
  const double dlnE = -(1.5 * m->om / (a * a * a) + m->ok / (a * a)) / E2;
  dy[0] = y[1];
  dy[1] = -(2.0 + dlnE) * y[1] + 1.5 * m->om / (a * a * a) / E2 * y[0];
}

static double growth_ode(const Model* m, double a_end) {
  const double xi = log(1e-6), xe = log(a_end);
  const int n = 40000;
  const double h = (xe - xi) / n;
  double y[2] = {exp(xi), exp(xi)}, k1[2], k2[2], k3[2], k4[2], t[2];
  for (int i = 0; i < n; ++i) {
    const double x = xi + i * h;
    grow_rhs(m, x, y, k1);
    for (int j = 0; j < 2; ++j) t[j] = y[j] + 0.5 * h * k1[j];
    grow_rhs(m, x + 0.5 * h, t, k2);
    for (int j = 0; j < 2; ++j) t[j] = y[j] + 0.5 * h * k2[j];
    grow_rhs(m, x + 0.5 * h, t, k3);
    for (int j = 0; j < 2; ++j) t[j] = y[j] + h * k3[j];
    grow_rhs(m, x + h, t, k4);
    for (int j = 0; j < 2; ++j) y[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
  }
  return y[0]; /* D -> a as a -> 0 */
}

/* ---- integrate the background a(t) and the shell R(t) together;
 *      return a when the shell collapses ----
 *   a'' = -Om/(2a^2) + OL a,   R'' = -Om/(2R^2) + OL R   (GM = Om/2, Lagrangian radius 1)
 *   Initial conditions: growing mode in the matter era,
 *     R = a (1 - delta/3), R' = R H (1 - delta/(3(1+delta))) */
static void sc_rhs(const Model* m, const double y[4], double dy[4]) {
  dy[0] = y[1];
  dy[1] = -0.5 * m->om / (y[0] * y[0]) + m->ol * y[0];
  dy[2] = y[3];
  dy[3] = -0.5 * m->om / (y[2] * y[2]) + m->ol * y[2];
}

static double collapse_a(const Model* m, double a_i, double delta_i) {
  const double H = sqrt(m->om / (a_i * a_i * a_i) + m->ok / (a_i * a_i) + m->ol);
  const double R = a_i * cbrt(1.0 / (1.0 + delta_i));
  double y[4] = {a_i, a_i * H, R, R * H * (1.0 - delta_i / (3.0 * (1.0 + delta_i)))};
  double k1[4], k2[4], k3[4], k4[4], t[4], Rmax = R;
  while (y[2] > 1e-6 * Rmax) {
    if (y[0] > 1e3) return INFINITY;
    const double h = 2e-3 * fmin(y[0] / fabs(y[1]), y[2] / fabs(y[3]));
    sc_rhs(m, y, k1);
    for (int j = 0; j < 4; ++j) t[j] = y[j] + 0.5 * h * k1[j];
    sc_rhs(m, t, k2);
    for (int j = 0; j < 4; ++j) t[j] = y[j] + 0.5 * h * k2[j];
    sc_rhs(m, t, k3);
    for (int j = 0; j < 4; ++j) t[j] = y[j] + h * k3[j];
    sc_rhs(m, t, k4);
    for (int j = 0; j < 4; ++j) y[j] += h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
    if (y[2] > Rmax) Rmax = y[2];
  }
  return y[0];
}

/* find by bisection the initial perturbation that collapses at a_c;
 * return the linearly extrapolated delta_c */
static double delta_crit_ode(const Model* m, double a_c) {
  const double a_i = 1e-7;
  double lo = 1e-9, hi = 1e-5; /* delta_i */
  for (int it = 0; it < 48; ++it) {
    const double mid = sqrt(lo * hi);
    if (collapse_a(m, a_i, mid) > a_c) lo = mid;
    else hi = mid;
  }
  const double di = sqrt(lo * hi);
  /* the second-order error of the initial conditions, ~ delta_i ~ 1e-7, is negligible */
  return di / a_i * growth_ode(m, a_c);
}

static void set(cosm_param* cp, double om, double ol) {
  cp->omega_cdm = 0.8 * om;
  cp->omega_b = 0.15 * om;
  cp->omega_nu = 0.05 * om;
  cp->omega_de = ol;
  cp->hubble = 0.7;
}

int main(void) {
  cosm_param cp;
  char name[128];

  /* ---------------- EdS ---------------- */
  set(&cp, 1.0, 0.0);
  puts("== Einstein-de Sitter");
  for (double z = 0.0; z <= 10.0; z += 2.5) {
    const double a = 1.0 / (1.0 + z);
    snprintf(name, sizeof name, "t(z=%g)", z);
    check(name, FLRW_cosmic_time_z(&cp, z), 2.0 / 3.0 * pow(a, 1.5), 1e-12);
    snprintf(name, sizeof name, "tau(z=%g)", z);
    check(name, FLRW_conformal_time_z(&cp, z), 2.0 * sqrt(a), 1e-12);
    snprintf(name, sizeof name, "DL(z=%g)", z);
    check(name, FLRW_DL_z(&cp, z), 2.0 * (1.0 + z) * (1.0 - sqrt(a)), 1e-12);
    snprintf(name, sizeof name, "D(z=%g)", z);
    check(name, FLRW_growth_factor_z(&cp, z), a, 1e-12);
    snprintf(name, sizeof name, "delta_c(z=%g)", z);
    check(name, FLRW_delta_crit_z(&cp, z), 0.6 * pow(1.5 * PI, 2.0 / 3.0), 1e-10);
    snprintf(name, sizeof name, "Delta_vir(z=%g)", z);
    check(name, FLRW_virial_Delta_z(&cp, z), 18.0 * PI * PI, 1e-10);
  }

  /* ---------------- open universe, Lambda = 0 ---------------- */
  set(&cp, 0.3, 0.0);
  puts("== open, Om=0.3, OL=0");
  for (double z = 0.5; z <= 8.0; z *= 2.0) {
    const double om = 0.3;
    /* Mattig's formula */
    const double dl = 2.0 * (om * z + (om - 2.0) * (sqrt(1.0 + om * z) - 1.0)) / (om * om);
    snprintf(name, sizeof name, "DL(z=%g)", z);
    check(name, FLRW_DL_z(&cp, z), dl, 1e-11);
    snprintf(name, sizeof name, "DA(z=%g)", z);
    check(name, FLRW_DA_z(&cp, z), dl / ((1.0 + z) * (1.0 + z)), 1e-11);
  }
  {
    /* t0 = [1/(1-Om) - Om/(2(1-Om)^{3/2}) acosh(2/Om - 1)] */
    const double om = 0.3, w = 1.0 - om;
    check("t0", FLRW_cosmic_time_z(&cp, 0.0),
          1.0 / w - om / (2.0 * pow(w, 1.5)) * acosh(2.0 / om - 1.0), 1e-12);
  }

  /* ---------------- flat LCDM ---------------- */
  set(&cp, 0.3, 0.7);
  puts("== flat LCDM, Om=0.3, OL=0.7");
  for (double z = 0.0; z <= 9.0; z += 3.0) {
    const double a = 1.0 / (1.0 + z);
    const double t = 2.0 / (3.0 * sqrt(0.7)) * asinh(sqrt(0.7 / 0.3) * pow(a, 1.5));
    snprintf(name, sizeof name, "t(z=%g)", z);
    check(name, FLRW_cosmic_time_z(&cp, z), t, 1e-12);
  }

  /* ---------------- growth factor and delta_c vs ODE ---------------- */
  const double models[][2] = {{0.3, 0.7}, {0.3, 0.0}, {0.3, 0.5}, {0.25, 0.9}, {2.0, 0.0}};
  for (size_t k = 0; k < sizeof models / sizeof models[0]; ++k) {
    const Model m = {models[k][0], 1.0 - models[k][0] - models[k][1], models[k][1]};
    set(&cp, m.om, m.ol);
    printf("== ODE comparison, Om=%g, OL=%g\n", m.om, m.ol);
    const double D0 = growth_ode(&m, 1.0);
    for (double z = 0.0; z <= 4.0; z += 2.0) {
      const double a = 1.0 / (1.0 + z);
      snprintf(name, sizeof name, "D(z=%g) vs ODE", z);
      check(name, FLRW_growth_factor_z(&cp, z), growth_ode(&m, a) / D0, 1e-9);
      snprintf(name, sizeof name, "delta_c(z=%g) vs ODE", z);
      check(name, FLRW_delta_crit_z(&cp, z), delta_crit_ode(&m, a), 5e-5);
    }
  }

  /* ---------------- Delta_vir: Bryan & Norman (1998) ---------------- */
  set(&cp, 0.3, 0.7);
  puts("== Delta_vir vs Bryan & Norman (flat, ~1%)");
  for (double z = 0.0; z <= 3.0; z += 1.0) {
    const double omz = FLRW_omega_m_z(&cp, z), x = omz - 1.0;
    snprintf(name, sizeof name, "Delta_vir(z=%g)", z);
    check(name, FLRW_virial_Delta_z(&cp, z), (18.0 * PI * PI + 82.0 * x - 39.0 * x * x) / omz,
          1.5e-2);
  }
  set(&cp, 0.3, 0.0);
  puts("== Delta_vir vs Bryan & Norman (open, ~1%)");
  for (double z = 0.0; z <= 3.0; z += 1.0) {
    const double omz = FLRW_omega_m_z(&cp, z), x = omz - 1.0;
    snprintf(name, sizeof name, "Delta_vir(z=%g)", z);
    check(name, FLRW_virial_Delta_z(&cp, z), (18.0 * PI * PI + 60.0 * x - 32.0 * x * x) / omz,
          1.5e-2);
  }

  /* ---------------- inverse functions and consistency of the t, tau versions ---------------- */
  const double models2[][2] = {{0.3, 0.7}, {0.3, 0.0}, {0.3, 0.5}, {1.5, 0.2}};
  for (size_t k = 0; k < sizeof models2 / sizeof models2[0]; ++k) {
    set(&cp, models2[k][0], models2[k][1]);
    printf("== round trip, Om=%g, OL=%g\n", models2[k][0], models2[k][1]);
    const double zs[] = {-0.5, 0.0, 0.3, 2.0, 30.0, 1e4};
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
      snprintf(name, sizeof name, "D_tau vs D_z (z=%g)", z);
      check(name, FLRW_growth_factor_tau(&cp, tau), FLRW_growth_factor_z(&cp, z), 1e-10);
      snprintf(name, sizeof name, "delta_c_t vs delta_c_z (z=%g)", z);
      check(name, FLRW_delta_crit_t(&cp, t), FLRW_delta_crit_z(&cp, z), 1e-9);
      snprintf(name, sizeof name, "Delta_vir_tau vs _z (z=%g)", z);
      check(name, FLRW_virial_Delta_tau(&cp, tau), FLRW_virial_Delta_z(&cp, z), 1e-9);
    }
  }

  /* ---------------- outside the domain ---------------- */
  set(&cp, 0.3, 0.7);
  puts("== domain errors (expect NaN)");
  {
    const double tau_inf = 2.0 * FLRW_conformal_time_z(&cp, -1.0 + 1e-12);
    const double v[] = {FLRW_cosmic_time_z(&cp, -1.5), FLRW_scale_factor_t(&cp, -1.0),
                        FLRW_scale_factor_tau(&cp, tau_inf)};
    for (int i = 0; i < 3; ++i) {
      ++ntest;
      if (!isnan(v[i])) ++nfail;
      printf("%-4s NaN check %d: %g\n", isnan(v[i]) ? "ok" : "FAIL", i, v[i]);
    }
  }

  printf("\n%d / %d tests passed\n", ntest - nfail, ntest);
  return nfail != 0;
}
