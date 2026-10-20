#ifndef FLRW_KIT_HPP
#define FLRW_KIT_HPP
/* FLRW_kit: header-only C++ wrapper around the C API in flrw_kit.h
 *
 *   flrw::Cosmology c(0.26, 0.049, 0.69, 0.677, 0.0014);  // cdm, b, de, h, nu
 *   double t  = c.cosmic_time_z(1.0) * c.hubble_time_Gyr();
 *   std::vector<double> dl = c.evaluate(&flrw::Cosmology::DL_z, {0.5, 1.0, 2.0});
 *
 * Differences from the C API
 *   - The parameters are checked in the constructor (Om > 0, h > 0, all finite);
 *     std::invalid_argument is thrown otherwise.
 *   - Where a C function would return NaN (argument out of domain, root not found, ...),
 *     the member function throws std::domain_error instead. Call the C functions on param()
 *     directly if NaN is preferred.
 *
 * Units and conventions are those of flrw_kit.h (time in 1/H0, distance in c/H0).
 * A Cosmology object is immutable, so it can be shared between threads.
 * Requires C++11. Link with -lflrw_kit -lm.
 */
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
            double omega_nu = 0.0) {
    cp_.omega_cdm = omega_cdm;
    cp_.omega_b = omega_b;
    cp_.omega_de = omega_de;
    cp_.hubble = hubble;
    cp_.omega_nu = omega_nu;
    validate();
  }

  const cosm_param& param() const { return cp_; }

  /* ---- auxiliary quantities ---- */
  double omega_m() const { return FLRW_omega_m(&cp_); }
  double omega_k() const { return FLRW_omega_k(&cp_); }
  double omega_de() const { return cp_.omega_de; }
  double hubble() const { return cp_.hubble; }
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

  /* ---- linear growth factor (D(z=0) = 1) ---- */
  double growth_factor_z(double z) const {
    return call(FLRW_growth_factor_z, z, "growth_factor_z");
  }
  double growth_factor_t(double t) const {
    return call(FLRW_growth_factor_t, t, "growth_factor_t");
  }
  double growth_factor_tau(double tau) const {
    return call(FLRW_growth_factor_tau, tau, "growth_factor_tau");
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

  /* ---- apply a member function to each element, e.g. evaluate(&Cosmology::DL_z, zs) ---- */
  typedef double (Cosmology::*Method)(double) const;

  std::vector<double> evaluate(Method f, const std::vector<double>& x) const {
    std::vector<double> y;
    y.reserve(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) y.push_back((this->*f)(x[i]));
    return y;
  }

 private:
  cosm_param cp_;

  void validate() const {
    const double om = cp_.omega_cdm + cp_.omega_b + cp_.omega_nu;
    if (!(om > 0.0) || !std::isfinite(om))
      throw std::invalid_argument("flrw::Cosmology: omega_m must be positive and finite");
    if (!std::isfinite(cp_.omega_de))
      throw std::invalid_argument("flrw::Cosmology: omega_de must be finite");
    if (!(cp_.hubble > 0.0) || !std::isfinite(cp_.hubble))
      throw std::invalid_argument("flrw::Cosmology: hubble (h) must be positive and finite");
  }

  double call(double (*f)(const cosm_param*, double), double x, const char* name) const {
    const double y = f(&cp_, x);
    if (std::isnan(y))
      throw std::domain_error(std::string("flrw::Cosmology::") + name +
                              ": no result for argument " + std::to_string(x));
    return y;
  }
};

}  // namespace flrw

#endif /* FLRW_KIT_HPP */
