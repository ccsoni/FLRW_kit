/* Tests for the C++ wrapper flrw_kit.hpp
 *   - every member function returns exactly the value of the corresponding C function,
 *     including the (k, X) overloads of growth_factor_cb_X / growth_factor_m_X
 *   - constructor defaults agree with FLRW_param_init, with_mnu_sum agrees with
 *     FLRW_set_mnu_from_sum, and the cosm_param constructor agrees with the value constructor
 *   - std::invalid_argument for invalid parameters, std::domain_error where C returns NaN
 *   - evaluate() over a vector of arguments and over a vector of wavenumbers
 */
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "flrw_kit.hpp"

namespace {

int nfail = 0, ntest = 0;

void check(const char* name, bool ok) {
  ++ntest;
  if (!ok) ++nfail;
  std::printf("%-4s %s\n", ok ? "ok" : "FAIL", name);
}

template <class E, class F>
void expect_throw(const char* name, F f) {
  bool thrown = false;
  try {
    f();
  } catch (const E& e) {
    thrown = true;
    std::printf("     (%s)\n", e.what());
  } catch (...) {
  }
  check(name, thrown);
}

typedef double (*CFn)(const cosm_param*, double);
typedef double (*CFnK)(const cosm_param*, double, double);

struct Pair {
  const char* name;
  flrw::Cosmology::Method cpp;
  CFn c;
  char arg; /* 'z', 't' or 'u' (tau) */
};

struct PairK {
  const char* name;
  flrw::Cosmology::MethodK cpp;
  CFnK c;
  char arg;
};

typedef flrw::Cosmology C;

const Pair pairs[] = {
    {"hubble_z", &C::hubble_z, FLRW_hubble_z, 'z'},
    {"omega_m_z", &C::omega_m_z, FLRW_omega_m_z, 'z'},
    {"cosmic_time_z", &C::cosmic_time_z, FLRW_cosmic_time_z, 'z'},
    {"conformal_time_z", &C::conformal_time_z, FLRW_conformal_time_z, 'z'},
    {"redshift_t", &C::redshift_t, FLRW_redshift_t, 't'},
    {"scale_factor_t", &C::scale_factor_t, FLRW_scale_factor_t, 't'},
    {"redshift_tau", &C::redshift_tau, FLRW_redshift_tau, 'u'},
    {"scale_factor_tau", &C::scale_factor_tau, FLRW_scale_factor_tau, 'u'},
    {"DC_z", &C::DC_z, FLRW_DC_z, 'z'},
    {"DM_z", &C::DM_z, FLRW_DM_z, 'z'},
    {"DA_z", &C::DA_z, FLRW_DA_z, 'z'},
    {"DL_z", &C::DL_z, FLRW_DL_z, 'z'},
    {"DC_t", &C::DC_t, FLRW_DC_t, 't'},
    {"DM_t", &C::DM_t, FLRW_DM_t, 't'},
    {"DA_t", &C::DA_t, FLRW_DA_t, 't'},
    {"DL_t", &C::DL_t, FLRW_DL_t, 't'},
    {"DC_tau", &C::DC_tau, FLRW_DC_tau, 'u'},
    {"DM_tau", &C::DM_tau, FLRW_DM_tau, 'u'},
    {"DA_tau", &C::DA_tau, FLRW_DA_tau, 'u'},
    {"DL_tau", &C::DL_tau, FLRW_DL_tau, 'u'},
    {"growth_factor_z", &C::growth_factor_z, FLRW_growth_factor_z, 'z'},
    {"growth_factor_t", &C::growth_factor_t, FLRW_growth_factor_t, 't'},
    {"growth_factor_tau", &C::growth_factor_tau, FLRW_growth_factor_tau, 'u'},
    {"growth_factor_cb_z(z)", &C::growth_factor_cb_z, FLRW_growth_factor_cb_z, 'z'},
    {"growth_factor_cb_t(t)", &C::growth_factor_cb_t, FLRW_growth_factor_cb_t, 't'},
    {"growth_factor_cb_tau(tau)", &C::growth_factor_cb_tau, FLRW_growth_factor_cb_tau, 'u'},
    {"growth_factor_m_z(z)", &C::growth_factor_m_z, FLRW_growth_factor_m_z, 'z'},
    {"growth_factor_m_t(t)", &C::growth_factor_m_t, FLRW_growth_factor_m_t, 't'},
    {"growth_factor_m_tau(tau)", &C::growth_factor_m_tau, FLRW_growth_factor_m_tau, 'u'},
    {"delta_crit_z", &C::delta_crit_z, FLRW_delta_crit_z, 'z'},
    {"delta_crit_t", &C::delta_crit_t, FLRW_delta_crit_t, 't'},
    {"delta_crit_tau", &C::delta_crit_tau, FLRW_delta_crit_tau, 'u'},
    {"virial_Delta_z", &C::virial_Delta_z, FLRW_virial_Delta_z, 'z'},
    {"virial_Delta_t", &C::virial_Delta_t, FLRW_virial_Delta_t, 't'},
    {"virial_Delta_tau", &C::virial_Delta_tau, FLRW_virial_Delta_tau, 'u'},
};

const PairK pairs_k[] = {
    {"growth_factor_cb_z(k, z)", &C::growth_factor_cb_z, FLRW_growth_factor_cb_k_z, 'z'},
    {"growth_factor_cb_t(k, t)", &C::growth_factor_cb_t, FLRW_growth_factor_cb_k_t, 't'},
    {"growth_factor_cb_tau(k, tau)", &C::growth_factor_cb_tau, FLRW_growth_factor_cb_k_tau, 'u'},
    {"growth_factor_m_z(k, z)", &C::growth_factor_m_z, FLRW_growth_factor_m_k_z, 'z'},
    {"growth_factor_m_t(k, t)", &C::growth_factor_m_t, FLRW_growth_factor_m_k_t, 't'},
    {"growth_factor_m_tau(k, tau)", &C::growth_factor_m_tau, FLRW_growth_factor_m_k_tau, 'u'},
};

double arg_of(const cosm_param& cp, char kind, double z) {
  if (kind == 't') return FLRW_cosmic_time_z(&cp, z);
  if (kind == 'u') return FLRW_conformal_time_z(&cp, z);
  return z;
}

}  // namespace

int main() {
  /* ---------------- wrapper == C API ---------------- */
  const C models[] = {
      C(0.26, 0.05, 0.69, 0.68),                                       /* massless neutrinos */
      C::with_mnu_sum(0.26, 0.05, 0.69, 0.68, 0.06, FLRW_NU_NORMAL),   /* minimal masses */
      C(0.35, 0.05, 0.5, 0.7, {{0.1, 0.1, 0.1}}, 2.7255, 3.0, 0.0),    /* no photons, curved */
  };
  for (const C& c : models) {
    const cosm_param& cp = c.param();
    std::printf("== wrapper vs C, m_nu = (%g, %g, %g), Ok = %.3g\n", cp.m_nu[0], cp.m_nu[1],
                cp.m_nu[2], c.omega_k());
    for (const Pair& p : pairs) {
      bool same = true;
      for (double z : {0.0, 0.5, 3.0}) {
        const double x = arg_of(cp, p.arg, z);
        same = same && ((c.*p.cpp)(x) == p.c(&cp, x));
      }
      check(p.name, same);
    }
    for (const PairK& p : pairs_k) {
      bool same = true;
      for (double k : {0.0, 300.0, 3e4})
        for (double z : {0.0, 2.0}) {
          const double x = arg_of(cp, p.arg, z);
          same = same && ((c.*p.cpp)(k, x) == p.c(&cp, k, x));
        }
      check(p.name, same);
    }
    check("omega_m, omega_cb, omega_nu, omega_r, omega_k",
          c.omega_m() == FLRW_omega_m(&cp) && c.omega_cb() == FLRW_omega_cb(&cp) &&
              c.omega_nu() == FLRW_omega_nu(&cp) && c.omega_r() == FLRW_omega_r(&cp) &&
              c.omega_k() == FLRW_omega_k(&cp));
    check("hubble_time_Gyr, hubble_distance_Mpc",
          c.hubble_time_Gyr() == FLRW_hubble_time_Gyr(&cp) &&
              c.hubble_distance_Mpc() == FLRW_hubble_distance_Mpc(&cp));
  }

  /* ---------------- constructors ---------------- */
  std::puts("== constructors");
  {
    cosm_param ref;
    FLRW_param_init(&ref);
    const C a(0.26, 0.049, 0.69, 0.677);
    check("defaults: T_cmb, N_eff, omega_r, m_nu as FLRW_param_init",
          a.T_cmb() == ref.T_cmb && a.N_eff() == ref.N_eff && a.param().omega_r == ref.omega_r &&
              a.m_nu()[0] == 0.0 && a.m_nu()[1] == 0.0 && a.m_nu()[2] == 0.0);

    ref.omega_cdm = 0.26;
    ref.omega_b = 0.049;
    ref.omega_de = 0.69;
    ref.hubble = 0.677;
    FLRW_set_mnu_from_sum(&ref, 0.1, FLRW_NU_INVERTED);
    const C b = C::with_mnu_sum(0.26, 0.049, 0.69, 0.677, 0.1, FLRW_NU_INVERTED);
    check("with_mnu_sum == FLRW_set_mnu_from_sum", b.m_nu()[0] == ref.m_nu[0] &&
                                                       b.m_nu()[1] == ref.m_nu[1] &&
                                                       b.m_nu()[2] == ref.m_nu[2]);
    const C d(ref);
    check("cosm_param ctor == with_mnu_sum", d.DL_z(1.0) == b.DL_z(1.0) &&
                                                 d.growth_factor_cb_z(300.0, 0.0) ==
                                                     b.growth_factor_cb_z(300.0, 0.0));
    check("growth_factor_cb_z(z) is the small-scale limit",
          b.growth_factor_cb_z(1.0) == FLRW_growth_factor_cb_z(&ref, 1.0) &&
              b.growth_factor_z(1.0) == b.growth_factor_cb_z(1.0));
    const C eds(1.0, 0.0, 0.0, 0.7, {{0.0, 0.0, 0.0}}, 2.7255, 0.0, 0.0);
    check("EdS delta_c", std::fabs(eds.delta_crit_z(1.0) - 1.6864701998411) < 1e-8);
    check("EdS Delta_vir", std::fabs(eds.virial_Delta_z(1.0) - 177.65287921960845) < 1e-6);
  }

  /* ---------------- evaluate ---------------- */
  std::puts("== evaluate");
  {
    const C c = C::with_mnu_sum(0.26, 0.05, 0.69, 0.68, 0.15, FLRW_NU_DEGENERATE);
    const std::vector<double> zs = {0.1, 0.5, 1.0, 3.0};
    const std::vector<double> dl = c.evaluate(&C::DL_z, zs);
    bool same = dl.size() == zs.size();
    for (std::size_t i = 0; same && i < zs.size(); ++i) same = (dl[i] == c.DL_z(zs[i]));
    check("evaluate(&Cosmology::DL_z, zs)", same);
    const std::vector<double> ks = {1.0, 100.0, 1e4};
    const std::vector<double> dk = c.evaluate(&C::growth_factor_cb_z, ks, 0.5);
    same = dk.size() == ks.size();
    for (std::size_t i = 0; same && i < ks.size(); ++i)
      same = (dk[i] == FLRW_growth_factor_cb_k_z(&c.param(), ks[i], 0.5));
    check("evaluate(&Cosmology::growth_factor_cb_z, ks, z)", same);
    const std::vector<double> dz = c.evaluate(&C::growth_factor_cb_z, zs);
    check("evaluate(&Cosmology::growth_factor_cb_z, zs) picks the one-argument overload",
          dz.size() == zs.size() && dz[2] == c.growth_factor_cb_z(1.0));
    check("evaluate on empty vector", c.evaluate(&C::DL_z, std::vector<double>()).empty());
  }

  /* ---------------- exceptions ---------------- */
  std::puts("== exceptions");
  expect_throw<std::invalid_argument>("Ocb = 0", [] { C(0.0, 0.0, 0.7, 0.7); });
  expect_throw<std::invalid_argument>("h = 0", [] { C(0.25, 0.05, 0.7, 0.0); });
  expect_throw<std::invalid_argument>("negative mass",
                                      [] { C(0.25, 0.05, 0.7, 0.7, {{-0.1, 0.0, 0.0}}); });
  expect_throw<std::invalid_argument>("T_cmb = 0",
                                      [] { C(0.25, 0.05, 0.7, 0.7, {{0.0, 0.0, 0.0}}, 0.0); });
  expect_throw<std::invalid_argument>("OL = NaN", [] { C(0.25, 0.05, NAN, 0.7); });
  expect_throw<std::invalid_argument>("with_mnu_sum below the inverted minimum", [] {
    C::with_mnu_sum(0.25, 0.05, 0.7, 0.7, 0.06, FLRW_NU_INVERTED);
  });
  {
    const C c(0.25, 0.05, 0.7, 0.7);
    const double tau_inf = 2.0 * c.conformal_time_z(-1.0 + 1e-12);
    expect_throw<std::domain_error>("cosmic_time_z(-1.5)", [&] { c.cosmic_time_z(-1.5); });
    expect_throw<std::domain_error>("scale_factor_t(-1)", [&] { c.scale_factor_t(-1.0); });
    expect_throw<std::domain_error>("scale_factor_tau(beyond tau_max)",
                                    [&] { c.scale_factor_tau(tau_inf); });
    expect_throw<std::domain_error>("growth_factor_cb_z(k < 0, z)",
                                    [&] { c.growth_factor_cb_z(-1.0, 0.0); });
    expect_throw<std::domain_error>("evaluate propagates", [&] {
      c.evaluate(&C::DL_z, std::vector<double>{0.5, -2.0});
    });
    char name[128];
    std::snprintf(name, sizeof name, "no throw inside the domain (DL_z(1) = %.6f)", c.DL_z(1.0));
    check(name, true);
  }

  std::printf("\n%d / %d tests passed\n", ntest - nfail, ntest);
  return nfail != 0;
}
