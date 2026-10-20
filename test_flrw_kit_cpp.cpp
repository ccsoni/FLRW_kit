/* Tests for the C++ wrapper flrw_kit.hpp
 *   - every member function returns exactly the value of the corresponding C function
 *   - both constructors give the same parameters
 *   - std::invalid_argument for invalid parameters, std::domain_error where C returns NaN
 *   - evaluate() over a vector
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

struct Pair {
  const char* name;
  flrw::Cosmology::Method cpp;
  CFn c;
  char arg; /* 'z', 't' or 'u' (tau) */
};

const Pair pairs[] = {
    {"hubble_z", &flrw::Cosmology::hubble_z, FLRW_hubble_z, 'z'},
    {"omega_m_z", &flrw::Cosmology::omega_m_z, FLRW_omega_m_z, 'z'},
    {"cosmic_time_z", &flrw::Cosmology::cosmic_time_z, FLRW_cosmic_time_z, 'z'},
    {"conformal_time_z", &flrw::Cosmology::conformal_time_z, FLRW_conformal_time_z, 'z'},
    {"redshift_t", &flrw::Cosmology::redshift_t, FLRW_redshift_t, 't'},
    {"scale_factor_t", &flrw::Cosmology::scale_factor_t, FLRW_scale_factor_t, 't'},
    {"redshift_tau", &flrw::Cosmology::redshift_tau, FLRW_redshift_tau, 'u'},
    {"scale_factor_tau", &flrw::Cosmology::scale_factor_tau, FLRW_scale_factor_tau, 'u'},
    {"DC_z", &flrw::Cosmology::DC_z, FLRW_DC_z, 'z'},
    {"DM_z", &flrw::Cosmology::DM_z, FLRW_DM_z, 'z'},
    {"DA_z", &flrw::Cosmology::DA_z, FLRW_DA_z, 'z'},
    {"DL_z", &flrw::Cosmology::DL_z, FLRW_DL_z, 'z'},
    {"DC_t", &flrw::Cosmology::DC_t, FLRW_DC_t, 't'},
    {"DM_t", &flrw::Cosmology::DM_t, FLRW_DM_t, 't'},
    {"DA_t", &flrw::Cosmology::DA_t, FLRW_DA_t, 't'},
    {"DL_t", &flrw::Cosmology::DL_t, FLRW_DL_t, 't'},
    {"DC_tau", &flrw::Cosmology::DC_tau, FLRW_DC_tau, 'u'},
    {"DM_tau", &flrw::Cosmology::DM_tau, FLRW_DM_tau, 'u'},
    {"DA_tau", &flrw::Cosmology::DA_tau, FLRW_DA_tau, 'u'},
    {"DL_tau", &flrw::Cosmology::DL_tau, FLRW_DL_tau, 'u'},
    {"growth_factor_z", &flrw::Cosmology::growth_factor_z, FLRW_growth_factor_z, 'z'},
    {"growth_factor_t", &flrw::Cosmology::growth_factor_t, FLRW_growth_factor_t, 't'},
    {"growth_factor_tau", &flrw::Cosmology::growth_factor_tau, FLRW_growth_factor_tau, 'u'},
    {"delta_crit_z", &flrw::Cosmology::delta_crit_z, FLRW_delta_crit_z, 'z'},
    {"delta_crit_t", &flrw::Cosmology::delta_crit_t, FLRW_delta_crit_t, 't'},
    {"delta_crit_tau", &flrw::Cosmology::delta_crit_tau, FLRW_delta_crit_tau, 'u'},
    {"virial_Delta_z", &flrw::Cosmology::virial_Delta_z, FLRW_virial_Delta_z, 'z'},
    {"virial_Delta_t", &flrw::Cosmology::virial_Delta_t, FLRW_virial_Delta_t, 't'},
    {"virial_Delta_tau", &flrw::Cosmology::virial_Delta_tau, FLRW_virial_Delta_tau, 'u'},
};

}  // namespace

int main() {
  char name[128];

  /* ---------------- wrapper == C API ---------------- */
  const double models[][2] = {{0.3, 0.7}, {0.3, 0.0}, {0.3, 0.5}, {1.5, 0.2}};
  for (const auto& m : models) {
    const flrw::Cosmology c(0.8 * m[0], 0.15 * m[0], m[1], 0.7, 0.05 * m[0]);
    const cosm_param& cp = c.param();
    std::printf("== wrapper vs C, Om=%g, OL=%g\n", m[0], m[1]);
    for (const Pair& p : pairs) {
      bool same = true;
      for (double z : {0.0, 0.5, 2.0, 10.0}) {
        double x = z;
        if (p.arg == 't') x = FLRW_cosmic_time_z(&cp, z);
        if (p.arg == 'u') x = FLRW_conformal_time_z(&cp, z);
        same = same && ((c.*p.cpp)(x) == p.c(&cp, x));
      }
      check(p.name, same);
    }
    check("omega_m", c.omega_m() == FLRW_omega_m(&cp));
    check("omega_k", c.omega_k() == FLRW_omega_k(&cp));
    check("hubble_time_Gyr", c.hubble_time_Gyr() == FLRW_hubble_time_Gyr(&cp));
    check("hubble_distance_Mpc", c.hubble_distance_Mpc() == FLRW_hubble_distance_Mpc(&cp));
  }

  /* ---------------- constructors ---------------- */
  std::puts("== constructors");
  {
    cosm_param cp;
    cp.omega_cdm = 0.26;
    cp.omega_b = 0.049;
    cp.omega_de = 0.69;
    cp.hubble = 0.677;
    cp.omega_nu = 0.0014;
    const flrw::Cosmology a(cp), b(0.26, 0.049, 0.69, 0.677, 0.0014);
    check("cosm_param ctor == value ctor", a.DL_z(1.0) == b.DL_z(1.0) &&
                                               a.omega_m() == b.omega_m() &&
                                               a.hubble() == b.hubble());
    const flrw::Cosmology eds(1.0, 0.0, 0.0, 0.7); /* omega_nu defaults to 0 */
    check("EdS delta_c", std::fabs(eds.delta_crit_z(1.0) - 1.6864701998411) < 1e-10);
    check("EdS Delta_vir", std::fabs(eds.virial_Delta_z(1.0) - 177.65287921960845) < 1e-8);
  }

  /* ---------------- evaluate ---------------- */
  std::puts("== evaluate");
  {
    const flrw::Cosmology c(0.25, 0.05, 0.7, 0.7);
    const std::vector<double> zs = {0.1, 0.5, 1.0, 3.0};
    const std::vector<double> dl = c.evaluate(&flrw::Cosmology::DL_z, zs);
    bool same = dl.size() == zs.size();
    for (std::size_t i = 0; same && i < zs.size(); ++i) same = (dl[i] == c.DL_z(zs[i]));
    check("evaluate(&Cosmology::DL_z, zs)", same);
    check("evaluate on empty vector",
          c.evaluate(&flrw::Cosmology::DL_z, std::vector<double>()).empty());
  }

  /* ---------------- exceptions ---------------- */
  std::puts("== exceptions");
  expect_throw<std::invalid_argument>("Om = 0", [] { flrw::Cosmology(0.0, 0.0, 0.7, 0.7); });
  expect_throw<std::invalid_argument>("Om < 0", [] { flrw::Cosmology(-0.3, 0.0, 0.7, 0.7); });
  expect_throw<std::invalid_argument>("h = 0", [] { flrw::Cosmology(0.25, 0.05, 0.7, 0.0); });
  expect_throw<std::invalid_argument>("OL = NaN",
                                      [] { flrw::Cosmology(0.25, 0.05, NAN, 0.7); });
  {
    const flrw::Cosmology c(0.25, 0.05, 0.7, 0.7);
    const double tau_inf = 2.0 * c.conformal_time_z(-1.0 + 1e-12);
    expect_throw<std::domain_error>("cosmic_time_z(-1.5)", [&] { c.cosmic_time_z(-1.5); });
    expect_throw<std::domain_error>("scale_factor_t(-1)", [&] { c.scale_factor_t(-1.0); });
    expect_throw<std::domain_error>("scale_factor_tau(beyond tau_max)",
                                    [&] { c.scale_factor_tau(tau_inf); });
    expect_throw<std::domain_error>("evaluate propagates", [&] {
      c.evaluate(&flrw::Cosmology::DL_z, std::vector<double>{0.5, -2.0});
    });
    std::snprintf(name, sizeof name, "no throw inside the domain (DL_z(1) = %.6f)", c.DL_z(1.0));
    check(name, true);
  }

  std::printf("\n%d / %d tests passed\n", ntest - nfail, ntest);
  return nfail != 0;
}
