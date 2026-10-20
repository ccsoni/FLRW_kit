# FLRW_kit — a C library for FLRW background, distances, linear growth and spherical collapse

[日本語](README_ja.md)

For a given set of cosmological parameters, FLRW_kit computes the cosmic time, conformal time, cosmological distances, the linear growth factor, the critical overdensity $\delta_c$ and the virial overdensity $\Delta_{\rm vir}$.
It needs only C11 and libm. The library keeps no internal state, so its functions can be called from inside OpenMP parallel regions.

## Build

```sh
make            # builds libflrw_kit.a
make test       # runs the tests
```

In your code, `#include "flrw_kit.h"` and link with `-L<dir> -lflrw_kit -lm`.

## Cosmological model

```c
typedef struct cosm_param {
  double omega_cdm; /* dark matter */
  double omega_b;   /* baryons */
  double omega_de;  /* dark energy (cosmological constant) */
  double hubble;    /* h  (H0 = 100 h km/s/Mpc) */
  double omega_nu;  /* neutrinos */
} cosm_param;
```

$$
E(a)^2 \equiv \left(\frac{H}{H_0}\right)^2 = \Omega_m a^{-3} + \Omega_k a^{-2} + \Omega_\Lambda
$$

- $\Omega_m = \Omega_{\rm cdm} + \Omega_b + \Omega_\nu$ (neutrinos are treated as non-relativistic matter)
- $\Omega_\Lambda$ = `omega_de` ($w=-1$)
- $\Omega_k = 1 - \Omega_m - \Omega_\Lambda$ (the curvature follows from the other parameters)
- Radiation is not included.

## Units and conventions

- Time $t$ and conformal time $\tau$ are in units of $1/H_0$, and distances in units of $c/H_0$. Since $c=1$, $\tau$ and comoving distance share the same unit. $t$ and $\tau$ are measured from the Big Bang ($a=0$).
- $h$ does not appear in the normalized quantities. To convert to physical units, use `FLRW_hubble_time_Gyr()` ($1/H_0$ in Gyr) and `FLRW_hubble_distance_Mpc()` ($c/H_0$ in Mpc).
- The comoving distance is the line-of-sight distance $D_C = \int_0^z dz'/E(z')$. The angular diameter and luminosity distances use the transverse comoving distance $D_M$: $D_A = D_M/(1+z)$ and $D_L = (1+z)D_M$.
  - $D_M = \sinh(\sqrt{\Omega_k}D_C)/\sqrt{\Omega_k}$ for $\Omega_k>0$, and $\sin(\sqrt{-\Omega_k}D_C)/\sqrt{-\Omega_k}$ for $\Omega_k<0$.
- The `_t` and `_tau` distance functions give the distance today to a source that emitted light at time $t$ (conformal time $\tau$). In particular, $D_C = \tau_0 - \tau$.
- On error, functions return NaN. Errors include invalid parameters, $z\le -1$, $t\le 0$, and a conformal time that is never reached.
- The universe is assumed to expand monotonically. The contracting phase of a recollapsing closed universe is not handled. Negative $z$ (the future) is allowed as long as the expansion continues.

## Functions

| Quantity | as a function of redshift $z$ | of time $t$ | of conformal time $\tau$ |
|---|---|---|---|
| cosmic time $t$ | `FLRW_cosmic_time_z` | | |
| conformal time $\tau$ | `FLRW_conformal_time_z` | | |
| redshift | | `FLRW_redshift_t` | `FLRW_redshift_tau` |
| scale factor | | `FLRW_scale_factor_t` | `FLRW_scale_factor_tau` |
| comoving distance $D_C$ | `FLRW_DC_z` | `FLRW_DC_t` | `FLRW_DC_tau` |
| transverse comoving distance $D_M$ | `FLRW_DM_z` | `FLRW_DM_t` | `FLRW_DM_tau` |
| angular diameter distance $D_A$ | `FLRW_DA_z` | `FLRW_DA_t` | `FLRW_DA_tau` |
| luminosity distance $D_L$ | `FLRW_DL_z` | `FLRW_DL_t` | `FLRW_DL_tau` |
| linear growth factor $D$ ($D(z{=}0)=1$) | `FLRW_growth_factor_z` | `FLRW_growth_factor_t` | `FLRW_growth_factor_tau` |
| critical overdensity $\delta_c$ | `FLRW_delta_crit_z` | `FLRW_delta_crit_t` | `FLRW_delta_crit_tau` |
| virial overdensity $\Delta_{\rm vir}$ | `FLRW_virial_Delta_z` | `FLRW_virial_Delta_t` | `FLRW_virial_Delta_tau` |

Every function has the form `double f(const cosm_param* cp, double x)`. Helper functions are also provided: `FLRW_omega_m`, `FLRW_omega_k`, `FLRW_hubble_z` ($E(z)$) and `FLRW_omega_m_z` ($\Omega_m(z)$).

```c
#include <stdio.h>
#include "flrw_kit.h"

int main(void) {
  cosm_param cp = {.omega_cdm = 0.26, .omega_b = 0.049, .omega_de = 0.69,
                   .hubble = 0.677, .omega_nu = 0.0014};
  double z = 1.0;
  printf("t(z)   = %g Gyr\n", FLRW_cosmic_time_z(&cp, z) * FLRW_hubble_time_Gyr(&cp));
  printf("DL(z)  = %g Mpc\n", FLRW_DL_z(&cp, z) * FLRW_hubble_distance_Mpc(&cp));
  printf("D(z)   = %g\n", FLRW_growth_factor_z(&cp, z));
  printf("dc(z)  = %g, Dvir(z) = %g\n", FLRW_delta_crit_z(&cp, z), FLRW_virial_Delta_z(&cp, z));
  return 0;
}
```

## Methods

### Time and distances

The substitution $u=\sqrt a$ removes the singular behavior at $a\to0$. The integrals are then evaluated by globally adaptive Gauss–Kronrod quadrature (7–15 points, relative accuracy $10^{-13}$).

$$
t(a) = \int_0^{\sqrt a}\frac{2u^2\,du}{\sqrt{\Omega_m + \Omega_k u^2 + \Omega_\Lambda u^6}},\qquad
\tau(a) = \int_0^{\sqrt a}\frac{2\,du}{\sqrt{\Omega_m + \Omega_k u^2 + \Omega_\Lambda u^6}}
$$

The inverse maps from $t$ and $\tau$ to $a$ use a bracketed Newton iteration in $\ln a$, with $dt/d\ln a = 1/E$ and $d\tau/d\ln a = 1/(aE)$.

### Linear growth factor

The growth factor uses the integral form of Heath (1977), which is exact for matter + cosmological constant + curvature:

$$
D_1(a) = \frac{5\Omega_m}{2}E(a)\int_0^a\frac{da'}{[a'E(a')]^3}\quad(D_1\to a\ \text{as}\ a\to0),\qquad
D(z) = D_1(a)/D_1(1)
$$

Neutrinos are assumed to cluster together with the rest of the matter. The scale dependence caused by neutrino free streaming is therefore not included.

### Spherical collapse ($\delta_c$, $\Delta_{\rm vir}$)

Instead of fitting formulae, FLRW_kit solves top-hat collapse exactly, including $\Lambda$ and curvature.
In units where the Lagrangian radius is 1 ($GM=\Omega_m/2$), the shell radius $R$ obeys

$$
\ddot R = -\frac{\Omega_m}{2R^2} + \Omega_\Lambda R,\qquad
\dot R^2 = \frac{\Omega_m}{R} + \Omega_\Lambda R^2 + K,\quad K = -\frac{\Omega_m}{R_{\rm ta}} - \Omega_\Lambda R_{\rm ta}^2
$$

For the background ($R=a$), $K=\Omega_k$. The perturbation is therefore described by the "curvature" difference $\delta K = \Omega_k - K$.

1. For a pure growing mode (the shell's Big Bang time equals the background's), the collapse time is twice the turnaround time:
   $$
   t_{\rm ta}(R_{\rm ta}) = \int_0^{\pi/2}\frac{2R_{\rm ta}^{3/2}\sin^2\theta\,d\theta}{\sqrt{\Omega_m - \Omega_\Lambda R_{\rm ta}^3\sin^2\theta(1+\sin^2\theta)}}
   $$
   $R_{\rm ta}$ is found by solving $2t_{\rm ta}(R_{\rm ta}) = t(a_c)$ with Brent's method. For $\Omega_\Lambda>0$, $R_{\rm ta}^3<\Omega_m/2\Omega_\Lambda$.
2. The initial linear density contrast is $\delta = \frac35(\delta K/\Omega_m)\,a$. Linearly extrapolated to the collapse time, it gives
   $$
   \delta_c = \frac35\,\frac{\delta K}{\Omega_m}\,D_1(a_c)
   $$
3. The virialization time is taken to be the collapse time. The virial radius $R_{\rm vir}=yR_{\rm ta}$ follows from energy conservation and the virial theorem, including the $\Lambda$ potential $U_\Lambda\propto -R^2$ (so $2K=-U_G+2U_\Lambda$):
   $$
   4\eta y^3 - 2(1+\eta)y + 1 = 0,\qquad \eta = \frac{\Omega_\Lambda R_{\rm ta}^3}{\Omega_m}
   $$
   The ratio to the mean matter density of the universe is $\Delta_{\rm vir} = (a_c/R_{\rm vir})^3$. For the ratio to the critical density, use $\Delta_{\rm vir}\,\Omega_m(z)$.

For Einstein–de Sitter, these reduce to $\delta_c = \frac35(3\pi/2)^{2/3} = 1.68647$ and $\Delta_{\rm vir} = 18\pi^2 = 177.65$.

## Tests (`test_flrw_kit.c`)

- Analytic solutions for Einstein–de Sitter, open $\Lambda=0$ (Mattig's formula, $t_0$) and flat ΛCDM ($t(a)$): relative error $\lesssim10^{-12}$.
- Growth factor: compared with a direct RK4 integration of the growth equation.
- $\delta_c$: compared with a direct integration of the background and shell equations of motion, with the initial amplitude bisected to match the collapse time. The cases include open, closed and flat models and $\Omega_m>1$.
- $\Delta_{\rm vir}$: compared with the Bryan & Norman (1998) fit, which is accurate to about 1%.
- Round trips through the inverse functions, consistency between the $z$, $t$ and $\tau$ versions, and NaN outside the domain.

## Cost

Each call performs numerical integration, and nothing is cached. Typical times are a few µs for time, distances and the growth factor, a few tens of µs for the inverse maps from $t$ and $\tau$, and a few hundred µs for $\delta_c$ and $\Delta_{\rm vir}$. If you need many evaluations, build a table on your side and interpolate.
