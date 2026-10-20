#ifndef FLRW_KIT_H
#define FLRW_KIT_H
/* FLRW_kit: background, distances, linear growth and spherical collapse in an FLRW universe
 *
 * Cosmological model
 *   E(a)^2 = (H/H0)^2 = Om a^-3 + Ok a^-2 + OL
 *     Om = omega_cdm + omega_b + omega_nu   (neutrinos are treated as non-relativistic matter)
 *     OL = omega_de                         (cosmological constant, w = -1)
 *     Ok = 1 - Om - OL                      (curvature follows from the other parameters)
 *   Radiation is not included.
 *
 * Units
 *   Time t and conformal time tau are in units of 1/H0, distances in units of c/H0
 *   Light travels a comoving distance c dtau in conformal time dtau, so in these units the
 *   comoving distance to a source that emitted at tau is numerically DC = tau0 - tau.
 *   t and tau are measured from the Big Bang (a = 0); today (z = 0) corresponds to t = t0,
 *   tau = tau0.
 *
 * Conventions
 *   - DC is the line-of-sight comoving distance int dz/E(z). The angular diameter and luminosity
 *     distances use the transverse comoving distance DM: DA = DM/(1+z), DL = (1+z) DM.
 *   - Distances as functions of t (tau) are measured today to a source that emitted at t (tau).
 *   - The growth factor is the linear growing mode normalized to D(z=0) = 1.
 *   - The critical overdensity delta_c and the virial overdensity Delta_vir come from exact top-hat
 *     spherical collapse including OL and curvature.
 *     Virialization is taken to occur at collapse (R -> 0). The virial radius follows from energy
 *     conservation and the virial theorem including the Lambda potential. Delta_vir is relative to
 *     the mean matter density of the universe at collapse.
 *   - The universe is assumed to expand monotonically (the contracting phase of a recollapsing
 *     closed universe is not handled).
 *
 * On error (invalid parameters, argument out of domain, root not found, ...) NaN is returned.
 * No internal state is kept, so all functions are thread-safe (callable inside OpenMP regions).
 */

typedef struct cosm_param {
  double omega_cdm; /* density parameter of dark matter */
  double omega_b;   /* density parameter of baryons */
  double omega_de;  /* density parameter of dark energy (cosmological constant) */
  double hubble;    /* h = H0 / (100 km/s/Mpc) */
  double omega_nu;  /* density parameter of neutrinos */
} cosm_param;

#ifdef __cplusplus
extern "C" {
#endif

/* ---- auxiliary quantities ---- */
double FLRW_omega_m(const cosm_param* cp);              /* Om = Ocdm + Ob + Onu */
double FLRW_omega_k(const cosm_param* cp);              /* Ok = 1 - Om - OL */
double FLRW_hubble_z(const cosm_param* cp, double z);   /* E(z) = H(z)/H0 */
double FLRW_omega_m_z(const cosm_param* cp, double z);  /* Om(z) = Om (1+z)^3 / E(z)^2 */
double FLRW_hubble_time_Gyr(const cosm_param* cp);      /* 1/H0 [Gyr] */
double FLRW_hubble_distance_Mpc(const cosm_param* cp);  /* c/H0 [Mpc] */

/* ---- time [1/H0] ---- */
double FLRW_cosmic_time_z(const cosm_param* cp, double z);
double FLRW_conformal_time_z(const cosm_param* cp, double z);

double FLRW_redshift_t(const cosm_param* cp, double t);
double FLRW_scale_factor_t(const cosm_param* cp, double t);
double FLRW_redshift_tau(const cosm_param* cp, double tau);
double FLRW_scale_factor_tau(const cosm_param* cp, double tau);

/* ---- distances [c/H0] ---- */
double FLRW_DC_z(const cosm_param* cp, double z);  /* comoving distance (line of sight) */
double FLRW_DM_z(const cosm_param* cp, double z);  /* transverse comoving distance */
double FLRW_DA_z(const cosm_param* cp, double z);  /* angular diameter distance */
double FLRW_DL_z(const cosm_param* cp, double z);  /* luminosity distance */

double FLRW_DC_t(const cosm_param* cp, double t);
double FLRW_DM_t(const cosm_param* cp, double t);
double FLRW_DA_t(const cosm_param* cp, double t);
double FLRW_DL_t(const cosm_param* cp, double t);

double FLRW_DC_tau(const cosm_param* cp, double tau);
double FLRW_DM_tau(const cosm_param* cp, double tau);
double FLRW_DA_tau(const cosm_param* cp, double tau);
double FLRW_DL_tau(const cosm_param* cp, double tau);

/* ---- linear growth factor (D(z=0) = 1) ---- */
double FLRW_growth_factor_z(const cosm_param* cp, double z);
double FLRW_growth_factor_t(const cosm_param* cp, double t);
double FLRW_growth_factor_tau(const cosm_param* cp, double tau);

/* ---- spherical collapse: critical overdensity delta_c
 *      (linearly extrapolated to virialization) ---- */
double FLRW_delta_crit_z(const cosm_param* cp, double z);
double FLRW_delta_crit_t(const cosm_param* cp, double t);
double FLRW_delta_crit_tau(const cosm_param* cp, double tau);

/* ---- spherical collapse: mean density within the virial radius
 *      / mean density of the universe ---- */
double FLRW_virial_Delta_z(const cosm_param* cp, double z);
double FLRW_virial_Delta_t(const cosm_param* cp, double t);
double FLRW_virial_Delta_tau(const cosm_param* cp, double tau);

#ifdef __cplusplus
}
#endif

#endif /* FLRW_KIT_H */
