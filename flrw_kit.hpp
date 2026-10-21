#ifndef FLRW_KIT_HPP
#define FLRW_KIT_HPP
/* FLRW_kit: header-only C++ wrapper around the C API in flrw_kit.h
 *
 *   flrw::Cosmology c(0.26, 0.049, 0.69, 0.677, {0.0, 0.0086, 0.050});  // cdm, b, de, h, m_nu
 *   auto c2 = flrw::Cosmology::with_mnu_sum(0.26, 0.049, 0.69, 0.677, 0.06, FLRW_NU_NORMAL);
 *   double t  = c.cosmic_time_z(1.0) * c.hubble_time_Gyr();
 *   double d  = c.growth_factor_cb_z(1.0);         // small-scale limit
 *   double dk = c.growth_factor_cb_z(300.0, 1.0);  // k = 300 H0/c (0.1 h/Mpc)
 *   std::vector<double> dl = c.evaluate(&flrw::Cosmology::DL_z, {0.5, 1.0, 2.0});
 *
 * Differences from the C API
 *   - Defaults are constructor arguments: massless neutrinos, T_cmb = 2.7255 K, N_eff = 3.044,
 *     photon density computed from T_cmb and h.
 *   - The parameters are checked in the constructor; std::invalid_argument is thrown if the
 *     C library rejects them.
 *   - Where a C function would return NaN (argument out of domain, root not found, ...),
 *     the member function throws std::domain_error instead. Call the C functions on param()
 *     directly if NaN is preferred.
 *   - growth_factor_cb_X and growth_factor_m_X (X = z, t, tau) are overloaded: with one argument
 *     they return the small-/large-scale limit, with two arguments (k, X) the scale-dependent
 *     growth factor FLRW_growth_factor_cb_k_X / FLRW_growth_factor_m_k_X.
 *
 * Units and conventions are those of flrw_kit.h (time in 1/H0, distance in c/H0, k in H0/c).
 * A Cosmology object is immutable, so it can be shared between threads.
 * Requires C++11. Link with -lflrw_kit -lm.
 */
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "flrw_kit.h"

namespace flrw {

class Cosmology {
 public:
  explicit Cosmology(const cosm_param& cp) : cp_(cp) { validate(); }

  Cosmology(double omega_cdm, double omega_b, double omega_de, double hubble,
            const std::array<double, 3>& m_nu = std::array<double, 3>{{0.0, 0.0, 0.0}},
            double T_cmb = FLRW_T_CMB_DEFAULT, double N_eff = FLRW_N_EFF_DEFAULT,
            double omega_r = FLRW_AUTO) {
    set(omega_cdm, omega_b, omega_de, hubble, T_cmb, N_eff, omega_r);
    for (int i = 0; i < 3; ++i) cp_.m_nu[i] = m_nu[i];
    validate();
  }

  /* neutrino masses from their sum [eV] and the ordering (FLRW_NU_NORMAL, FLRW_NU_INVERTED,
   * FLRW_NU_DEGENERATE) */
  static Cosmology with_mnu_sum(double omega_cdm, double omega_b, double omega_de, double hubble,
                                double sum_eV, int ordering = FLRW_NU_NORMAL,
                                double T_cmb = FLRW_T_CMB_DEFAULT,
                                double N_eff = FLRW_N_EFF_DEFAULT, double omega_r = FLRW_AUTO) {
    Cosmology c;
    c.set(omega_cdm, omega_b, omega_de, hubble, T_cmb, N_eff, omega_r);
    if (FLRW_set_mnu_from_sum(&c.cp_, sum_eV, ordering) != 0)
      throw std::invalid_argument("flrw::Cosmology::with_mnu_sum: sum " + std::to_string(sum_eV) +
                                  " eV is not allowed for this ordering");
    c.validate();
    return c;
  }

  const cosm_param& param() const { return cp_; }

  /* ---- parameters and auxiliary quantities ---- */
  double omega_m() const { return FLRW_omega_m(&cp_); }
  double omega_cb() const { return FLRW_omega_cb(&cp_); }
  double omega_nu() const { return FLRW_omega_nu(&cp_); }
  double omega_r() const { return FLRW_omega_r(&cp_); }
  double omega_k() const { return FLRW_omega_k(&cp_); }
  double omega_de() const { return cp_.omega_de; }
  double hubble() const { return cp_.hubble; }
  std::array<double, 3> m_nu() const {
    return std::array<double, 3>{{cp_.m_nu[0], cp_.m_nu[1], cp_.m_nu[2]}};
  }
  double T_cmb() const { return cp_.T_cmb; }
  double N_eff() const { return cp_.N_eff; }
  double hubble_z(double z) const { return call(FLRW_hubble_z, z, "hubble_z"); }
  double omega_m_z(double z) const { return call(FLRW_omega_m_z, z, "omega_m_z"); }
  double hubble_time_Gyr() const { return FLRW_hubble_time_Gyr(&cp_); }
  double hubble_distance_Mpc() const { return FLRW_hubble_distance_Mpc(&cp_); }

  /* ---- time [1/H0] ---- */
  double cosmic_time_z(double z) const { return call(FLRW_cosmic_time_z, z, "cosmic_time_z"); }
  double conformal_time_z(double z) const {
    return call(FLRW_conformal_time_z, z, "conformal_time_z");
  }
  double redshift_t(double t) const { return call(FLRW_redshift_t, t, "redshift_t"); }
  double scale_factor_t(double t) const { return call(FLRW_scale_factor_t, t, "scale_factor_t"); }
  double redshift_tau(double tau) const { return call(FLRW_redshift_tau, tau, "redshift_tau"); }
  double scale_factor_tau(double tau) const {
    return call(FLRW_scale_factor_tau, tau, "scale_factor_tau");
  }

  /* ---- distances [c/H0] ---- */
  double DC_z(double z) const { return call(FLRW_DC_z, z, "DC_z"); }
  double DM_z(double z) const { return call(FLRW_DM_z, z, "DM_z"); }
  double DA_z(double z) const { return call(FLRW_DA_z, z, "DA_z"); }
  double DL_z(double z) const { return call(FLRW_DL_z, z, "DL_z"); }

  double DC_t(double t) const { return call(FLRW_DC_t, t, "DC_t"); }
  double DM_t(double t) const { return call(FLRW_DM_t, t, "DM_t"); }
  double DA_t(double t) const { return call(FLRW_DA_t, t, "DA_t"); }
  double DL_t(double t) const { return call(FLRW_DL_t, t, "DL_t"); }

  double DC_tau(double tau) const { return call(FLRW_DC_tau, tau, "DC_tau"); }
  double DM_tau(double tau) const { return call(FLRW_DM_tau, tau, "DM_tau"); }
  double DA_tau(double tau) const { return call(FLRW_DA_tau, tau, "DA_tau"); }
  double DL_tau(double tau) const { return call(FLRW_DL_tau, tau, "DL_tau"); }

  /* ---- linear growth factor: growth_factor = small-scale limit ---- */
  double growth_factor_z(double z) const {
    return call(FLRW_growth_factor_z, z, "growth_factor_z");
  }
  double growth_factor_t(double t) const {
    return call(FLRW_growth_factor_t, t, "growth_factor_t");
  }
  double growth_factor_tau(double tau) const {
    return call(FLRW_growth_factor_tau, tau, "growth_factor_tau");
  }

  /* one argument: small-scale limit; (k, X): scale-dependent growth of the cdm + baryon contrast */
  double growth_factor_cb_z(double z) const {
    return call(FLRW_growth_factor_cb_z, z, "growth_factor_cb_z");
  }
  double growth_factor_cb_t(double t) const {
    return call(FLRW_growth_factor_cb_t, t, "growth_factor_cb_t");
  }
  double growth_factor_cb_tau(double tau) const {
    return call(FLRW_growth_factor_cb_tau, tau, "growth_factor_cb_tau");
  }
  double growth_factor_cb_z(double k, double z) const {
    return call(FLRW_growth_factor_cb_k_z, k, z, "growth_factor_cb_z(k, z)");
  }
  double growth_factor_cb_t(double k, double t) const {
    return call(FLRW_growth_factor_cb_k_t, k, t, "growth_factor_cb_t(k, t)");
  }
  double growth_factor_cb_tau(double k, double tau) const {
    return call(FLRW_growth_factor_cb_k_tau, k, tau, "growth_factor_cb_tau(k, tau)");
  }

  /* one argument: large-scale limit; (k, X): scale-dependent growth of the total matter contrast */
  double growth_factor_m_z(double z) const {
    return call(FLRW_growth_factor_m_z, z, "growth_factor_m_z");
  }
  double growth_factor_m_t(double t) const {
    return call(FLRW_growth_factor_m_t, t, "growth_factor_m_t");
  }
  double growth_factor_m_tau(double tau) const {
    return call(FLRW_growth_factor_m_tau, tau, "growth_factor_m_tau");
  }
  double growth_factor_m_z(double k, double z) const {
    return call(FLRW_growth_factor_m_k_z, k, z, "growth_factor_m_z(k, z)");
  }
  double growth_factor_m_t(double k, double t) const {
    return call(FLRW_growth_factor_m_k_t, k, t, "growth_factor_m_t(k, t)");
  }
  double growth_factor_m_tau(double k, double tau) const {
    return call(FLRW_growth_factor_m_k_tau, k, tau, "growth_factor_m_tau(k, tau)");
  }

  /* ---- spherical collapse: critical overdensity delta_c ---- */
  double delta_crit_z(double z) const { return call(FLRW_delta_crit_z, z, "delta_crit_z"); }
  double delta_crit_t(double t) const { return call(FLRW_delta_crit_t, t, "delta_crit_t"); }
  double delta_crit_tau(double tau) const {
    return call(FLRW_delta_crit_tau, tau, "delta_crit_tau");
  }

  /* ---- spherical collapse: virial overdensity relative to the mean matter density ---- */
  double virial_Delta_z(double z) const { return call(FLRW_virial_Delta_z, z, "virial_Delta_z"); }
  double virial_Delta_t(double t) const { return call(FLRW_virial_Delta_t, t, "virial_Delta_t"); }
  double virial_Delta_tau(double tau) const {
    return call(FLRW_virial_Delta_tau, tau, "virial_Delta_tau");
  }

  /* ---- apply a member function to each element:
   *      evaluate(&Cosmology::DL_z, zs), evaluate(&Cosmology::growth_factor_cb_z, ks, z) ---- */
  typedef double (Cosmology::*Method)(double) const;
  typedef double (Cosmology::*MethodK)(double, double) const;

  std::vector<double> evaluate(Method f, const std::vector<double>& x) const {
    std::vector<double> y;
    y.reserve(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) y.push_back((this->*f)(x[i]));
    return y;
  }

  std::vector<double> evaluate(MethodK f, const std::vector<double>& k, double x) const {
    std::vector<double> y;
    y.reserve(k.size());
    for (std::size_t i = 0; i < k.size(); ++i) y.push_back((this->*f)(k[i], x));
    return y;
  }

 private:
  cosm_param cp_;

  Cosmology() { FLRW_param_init(&cp_); }

  void set(double omega_cdm, double omega_b, double omega_de, double hubble, double T_cmb,
           double N_eff, double omega_r) {
    FLRW_param_init(&cp_);
    cp_.omega_cdm = omega_cdm;
    cp_.omega_b = omega_b;
    cp_.omega_de = omega_de;
    cp_.hubble = hubble;
    cp_.T_cmb = T_cmb;
    cp_.N_eff = N_eff;
    cp_.omega_r = omega_r;
  }

  void validate() const {
    if (std::isnan(FLRW_omega_k(&cp_)))
      throw std::invalid_argument(
          "flrw::Cosmology: invalid parameters (need omega_cdm + omega_b > 0, h > 0, T_cmb > 0, "
          "N_eff >= 0, m_nu >= 0, all finite)");
  }

  double call(double (*f)(const cosm_param*, double), double x, const char* name) const {
    const double y = f(&cp_, x);
    if (std::isnan(y))
      throw std::domain_error(std::string("flrw::Cosmology::") + name +
                              ": no result for argument " + std::to_string(x));
    return y;
  }

  double call(double (*f)(const cosm_param*, double, double), double k, double x,
              const char* name) const {
    const double y = f(&cp_, k, x);
    if (std::isnan(y))
      throw std::domain_error(std::string("flrw::Cosmology::") + name +
                              ": no result for k = " + std::to_string(k) +
                              ", argument " + std::to_string(x));
    return y;
  }
};

}  // namespace flrw

#endif /* FLRW_KIT_HPP */
