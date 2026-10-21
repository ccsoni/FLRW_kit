#ifndef FLRW_KIT_H
#define FLRW_KIT_H
/* FLRW_kit: background, distances, linear growth and spherical collapse in an FLRW universe
 *
 * Cosmological model
 *   E(a)^2 = (H/H0)^2 = Ocb a^-3 + Og a^-4 + rho_nu(a)/rho_c0 + Ok a^-2 + OL
 *     Ocb = omega_cdm + omega_b             (cold dark matter + baryons)
 *     Og  = omega_r                         (photons; FLRW_AUTO computes it from T_cmb and h)
 *     rho_nu(a)                             three neutrino species with masses m_nu[0..2] and a
 *                                           Fermi-Dirac momentum distribution at temperature
 *                                           T_nu = (4/11)^{1/3} (N_eff/3)^{1/4} T_cmb; each one
 *                                           turns from radiation into matter at m_i a ~ 3 k_B T_nu0
 *     OL  = omega_de                        (cosmological constant, w = -1)
 *     Ok  = 1 - Ocb - Og - Onu - OL         (curvature follows from the other parameters)
 *   "Matter" (FLRW_omega_m, Delta_vir) means cdm + baryons + the non-relativistic part
 *   rho_nu - 3 p_nu of the neutrinos, which vanishes for relativistic and massless species.
 *
 * Units
 *   Time t and conformal time tau are in units of 1/H0, distances in units of c/H0, and
 *   wavenumbers k in units of H0/c (k[H0/c] = k[h/Mpc] * FLRW_HUBBLE_DISTANCE_MPC_H).
 *   Light travels a comoving distance c dtau in conformal time dtau, so in these units the
 *   comoving distance to a source that emitted at tau is numerically DC = tau0 - tau.
 *   t and tau are measured from the Big Bang (a = 0); today (z = 0) corresponds to t = t0,
 *   tau = tau0.
 *
 * Conventions
 *   - DC is the line-of-sight comoving distance int dz/E(z). The angular diameter and luminosity
 *     distances use the transverse comoving distance DM: DA = DM/(1+z), DL = (1+z) DM.
 *   - Distances as functions of t (tau) are measured today to a source that emitted at t (tau).
 *   - Linear growth (sub-horizon, Newtonian; cdm and baryons as one pressureless fluid):
 *       growth_factor_cb   small-scale limit: neutrinos do not cluster, D(z=0) = 1
 *       growth_factor_m    large-scale limit: neutrinos cluster with cdm + baryons, D(z=0) = 1
 *       growth_factor      same as growth_factor_cb
 *       growth_factor_cb_k / growth_factor_m_k
 *                          scale-dependent growth of the cdm + baryon / total matter contrast
 *                          from a two-fluid model (neutrinos as a fluid with an effective sound
 *                          speed), normalized so that D_m(z=0) = 1 in the large-scale limit and
 *                          all k share the same initial amplitude; D < 1 at z = 0 for k > k_fs.
 *   - The critical overdensity delta_c and the virial overdensity Delta_vir come from top-hat
 *     spherical collapse of cdm + baryons, with photons, neutrinos and OL as smooth components.
 *     Virialization is taken to occur at collapse (R -> 0). delta_c is the linear contrast
 *     extrapolated with growth_factor_cb. Delta_vir is the mean matter density within the virial
 *     radius relative to the mean matter density of the universe at collapse.
 *   - The universe is assumed to expand monotonically (the contracting phase of a recollapsing
 *     closed universe is not handled).
 *
 * Parameters must be initialized with FLRW_param_init() or FLRW_PARAM_DEFAULTS, then set.
 * On error (invalid parameters, argument out of domain, root not found, ...) NaN is returned.
 * No internal state is kept, so all functions are thread-safe (callable inside OpenMP regions).
 */

typedef struct cosm_param {
  double omega_cdm; /* density parameter of cold dark matter */
  double omega_b;   /* density parameter of baryons */
  double omega_de;  /* density parameter of dark energy (cosmological constant) */
  double hubble;    /* h = H0 / (100 km/s/Mpc) */
  double m_nu[3];   /* neutrino mass eigenvalues [eV] (0 for a massless species) */
  double T_cmb;     /* CMB temperature today [K] */
  double N_eff;     /* effective number of neutrino species */
  double omega_r;   /* density parameter of photons, or FLRW_AUTO to compute it from T_cmb, h */
} cosm_param;

#define FLRW_AUTO (-1.0)
#define FLRW_T_CMB_DEFAULT 2.7255
#define FLRW_N_EFF_DEFAULT 3.044

/* designated initializers for the defaults (C99):
 *   cosm_param cp = {FLRW_PARAM_DEFAULTS, .omega_cdm = 0.26, .omega_b = 0.049, ...};
 * members not listed are zero (massless neutrinos). */
#define FLRW_PARAM_DEFAULTS \
  .T_cmb = FLRW_T_CMB_DEFAULT, .N_eff = FLRW_N_EFF_DEFAULT, .omega_r = FLRW_AUTO

/* c/H0 in h^-1 Mpc: D[h^-1 Mpc] = D[c/H0] * this, k[H0/c] = k[h/Mpc] * this */
#define FLRW_HUBBLE_DISTANCE_MPC_H 2997.92458

/* neutrino mass ordering for FLRW_set_mnu_from_sum() */
enum { FLRW_NU_NORMAL = 0, FLRW_NU_INVERTED = 1, FLRW_NU_DEGENERATE = 2 };

#ifdef __cplusplus
extern "C" {
#endif

/* ---- parameters ---- */
/* all members zero except T_cmb = 2.7255, N_eff = 3.044, omega_r = FLRW_AUTO */
void FLRW_param_init(cosm_param* cp);
/* set m_nu[0..2] = (m1, m2, m3) from the sum of the masses [eV] and the ordering, using the
 * oscillation mass splittings. Returns 0, or -1 if the sum is below the minimum allowed by the
 * ordering (about 0.059 eV normal, 0.099 eV inverted) or the arguments are invalid. */
int FLRW_set_mnu_from_sum(cosm_param* cp, double sum_eV, int ordering);

/* ---- auxiliary quantities (today unless a z argument is given) ---- */
double FLRW_omega_m(const cosm_param* cp);              /* Ocb + matter-like part of Onu */
double FLRW_omega_cb(const cosm_param* cp);             /* Ocdm + Ob */
double FLRW_omega_nu(const cosm_param* cp);             /* all neutrino species */
double FLRW_omega_r(const cosm_param* cp);              /* photons (FLRW_AUTO resolved) */
double FLRW_omega_k(const cosm_param* cp);              /* 1 - Ocb - Og - Onu - OL */
double FLRW_hubble_z(const cosm_param* cp, double z);   /* E(z) = H(z)/H0 */
double FLRW_omega_m_z(const cosm_param* cp, double z);  /* rho_m(z) / rho_crit(z) */
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

/* ---- linear growth factor: small-scale limit (D(z=0) = 1); growth_factor = growth_factor_cb */
double FLRW_growth_factor_z(const cosm_param* cp, double z);
double FLRW_growth_factor_t(const cosm_param* cp, double t);
double FLRW_growth_factor_tau(const cosm_param* cp, double tau);

double FLRW_growth_factor_cb_z(const cosm_param* cp, double z);
double FLRW_growth_factor_cb_t(const cosm_param* cp, double t);
double FLRW_growth_factor_cb_tau(const cosm_param* cp, double tau);

/* ---- linear growth factor: large-scale limit (D(z=0) = 1) ---- */
double FLRW_growth_factor_m_z(const cosm_param* cp, double z);
double FLRW_growth_factor_m_t(const cosm_param* cp, double t);
double FLRW_growth_factor_m_tau(const cosm_param* cp, double tau);

/* ---- scale-dependent linear growth factor, k in units of H0/c
 *      (normalized to the large-scale limit: D_m(k -> 0, z = 0) = 1) ---- */
double FLRW_growth_factor_cb_k_z(const cosm_param* cp, double k, double z);
double FLRW_growth_factor_cb_k_t(const cosm_param* cp, double k, double t);
double FLRW_growth_factor_cb_k_tau(const cosm_param* cp, double k, double tau);

double FLRW_growth_factor_m_k_z(const cosm_param* cp, double k, double z);
double FLRW_growth_factor_m_k_t(const cosm_param* cp, double k, double t);
double FLRW_growth_factor_m_k_tau(const cosm_param* cp, double k, double tau);

/* ---- spherical collapse: critical overdensity delta_c
 *      (linearly extrapolated to virialization with growth_factor_cb) ---- */
double FLRW_delta_crit_z(const cosm_param* cp, double z);
double FLRW_delta_crit_t(const cosm_param* cp, double t);
double FLRW_delta_crit_tau(const cosm_param* cp, double tau);

/* ---- spherical collapse: mean matter density within the virial radius
 *      / mean matter density of the universe ---- */
double FLRW_virial_Delta_z(const cosm_param* cp, double z);
double FLRW_virial_Delta_t(const cosm_param* cp, double t);
double FLRW_virial_Delta_tau(const cosm_param* cp, double tau);

#ifdef __cplusplus
}
#endif

#endif /* FLRW_KIT_H */
