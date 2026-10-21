# FLRW_kit — a C library for FLRW background, distances, linear growth and spherical collapse

[日本語](README_ja.md)

For a given set of cosmological parameters, FLRW_kit computes the cosmic time, conformal time, cosmological distances, the linear growth factor (in the small- and large-scale limits and as a function of wavenumber), the critical overdensity $`\delta_c`$ and the virial overdensity $`\Delta_{\mathrm{vir}}`$.
Radiation and massive neutrinos are included: each neutrino species turns from radiation into matter according to its Fermi–Dirac distribution.
It needs only C11 and libm. The library keeps no internal state, so its functions can be called from inside OpenMP parallel regions.

## Build

```sh
make            # builds libflrw_kit.a
make test       # runs the C and C++ tests (make test_c / make test_cpp for one of them)
```

In your code, `#include "flrw_kit.h"` and link with `-L<dir> -lflrw_kit -lm`.
From C++, either include `flrw_kit.h` directly (it is wrapped in `extern "C"`) or use the C++ wrapper `flrw_kit.hpp` described below.

`flrw_nu_table.h` is a generated table of the neutrino energy density and pressure. It is committed, so numpy is not needed to build. `make table` regenerates it with `tools/gen_nu_table.py` (requires python3 with numpy).

## Cosmological model

```c
typedef struct cosm_param {
  double omega_cdm; /* cold dark matter */
  double omega_b;   /* baryons */
  double omega_de;  /* dark energy (cosmological constant) */
  double hubble;    /* h  (H0 = 100 h km/s/Mpc) */
  double m_nu[3];   /* neutrino mass eigenvalues [eV] (0 for a massless species) */
  double T_cmb;     /* CMB temperature today [K]        default 2.7255 */
  double N_eff;     /* effective number of neutrinos     default 3.044 */
  double omega_r;   /* photon density parameter          default FLRW_AUTO (from T_cmb and h) */
} cosm_param;
```

```math
E(a)^2 \equiv \left(\frac{H}{H_0}\right)^2 = \Omega_{cb}\,a^{-3} + \Omega_\gamma\,a^{-4} + \frac{\rho_\nu(a)}{\rho_{c0}} + \Omega_k\,a^{-2} + \Omega_\Lambda
```

- $`\Omega_{cb} = \Omega_{\mathrm{cdm}} + \Omega_b`$; $`\Omega_\Lambda`$ = `omega_de` ($`w=-1`$).
- Photons: $`\Omega_\gamma h^2 = 2.47298\times10^{-5}\,(T_{\mathrm{cmb}}/2.7255\,\mathrm{K})^4`$ when `omega_r = FLRW_AUTO`; a non-negative `omega_r` is used as given.
- Neutrinos: three species with masses `m_nu[i]` and a Fermi–Dirac distribution at $`T_\nu = (4/11)^{1/3}(N_{\mathrm{eff}}/3)^{1/4}\,T_{\mathrm{cmb}}`$. Each species is radiation-like for $`m_i a \ll 3k_BT_{\nu0}`$ and matter-like afterwards; a massless species stays radiation. $`N_{\mathrm{eff}}=0`$ removes the neutrinos.
- $`\Omega_k = 1 - \Omega_{cb} - \Omega_\gamma - \Omega_\nu - \Omega_\Lambda`$ (the curvature follows from the other parameters). To get a flat model, choose `omega_de` accordingly, e.g. `1 - FLRW_omega_cb - FLRW_omega_r - FLRW_omega_nu`.
- "Matter" ($`\Omega_m`$, `FLRW_omega_m_z`, $`\Delta_{\mathrm{vir}}`$) means $`\rho_m = \rho_{cb} + \sum_i(\rho_{\nu,i} - 3p_{\nu,i})`$: cdm + baryons + the non-relativistic part of the neutrinos. It counts a non-relativistic species fully, and relativistic or massless species not at all.

### Initializing the parameters

A C struct cannot carry default values, so initialize it first with `FLRW_param_init()` or the initializer macro `FLRW_PARAM_DEFAULTS`. Both set `T_cmb = 2.7255`, `N_eff = 3.044`, `omega_r = FLRW_AUTO`, and zero for everything else (massless neutrinos).

```c
cosm_param cp;
FLRW_param_init(&cp);
cp.omega_cdm = 0.26; cp.omega_b = 0.049; cp.omega_de = 0.69; cp.hubble = 0.677;

/* or, with C99 designated initializers */
cosm_param cp2 = {FLRW_PARAM_DEFAULTS,
                  .omega_cdm = 0.26, .omega_b = 0.049, .omega_de = 0.69, .hubble = 0.677};
```

### Neutrino masses from the sum and the ordering

`FLRW_set_mnu_from_sum(&cp, sum_eV, ordering)` sets `m_nu[0..2]` = $`(m_1, m_2, m_3)`$ so that they add up to `sum_eV`, using the oscillation mass splittings $`\Delta m^2_{21} = 7.42\times10^{-5}\,\mathrm{eV}^2`$, $`\Delta m^2_{31} = 2.515\times10^{-3}\,\mathrm{eV}^2`$ (normal) and $`|\Delta m^2_{32}| = 2.498\times10^{-3}\,\mathrm{eV}^2`$ (inverted). The splittings can be changed at compile time with `-DFLRW_DM21_SQ=...`, `-DFLRW_DM31_SQ_NORMAL=...`, `-DFLRW_DM32_SQ_INVERTED=...`.

| `ordering` | masses | minimum sum |
|---|---|---|
| `FLRW_NU_NORMAL` | $`m_1 < m_2 < m_3`$ | 0.0588 eV |
| `FLRW_NU_INVERTED` | $`m_3 < m_1 < m_2`$ | 0.0992 eV |
| `FLRW_NU_DEGENERATE` | $`m_1 = m_2 = m_3 = \Sigma/3`$ | 0 |

It returns 0, or −1 if the sum is below the minimum for the ordering. With $`N_{\mathrm{eff}}=3.044`$ the non-relativistic neutrino density is $`\Omega_\nu h^2 = \Sigma m_\nu / 93.04\,\mathrm{eV}`$.

## Units and conventions

- Time $`t`$ and conformal time $`\tau`$ are in units of $`1/H_0`$, distances in units of $`c/H_0`$, and wavenumbers $`k`$ in units of $`H_0/c`$. Light travels a comoving distance $`c\,d\tau`$ in conformal time $`d\tau`$, so in these units the comoving distance to a source that emitted at $`\tau`$ is numerically $`D_C = \tau_0 - \tau`$. $`t`$ and $`\tau`$ are measured from the Big Bang ($`a=0`$).
- `FLRW_HUBBLE_DISTANCE_MPC_H` = 2997.92458 is $`c/H_0`$ in $`h^{-1}`$ Mpc. Multiply a distance by it to get $`h^{-1}`$ Mpc, and multiply a wavenumber in $`h\,\mathrm{Mpc}^{-1}`$ by it to get $`k`$ in $`H_0/c`$ (e.g. $`k = 0.1\,h\,\mathrm{Mpc}^{-1}`$ is $`k = 299.8`$). This conversion does not depend on $`h`$.
- To convert to physical units, multiply a time by `FLRW_hubble_time_Gyr()` ($`1/H_0 = 9.778\,h^{-1}`$ Gyr) and a distance by `FLRW_hubble_distance_Mpc()` ($`c/H_0 = 2997.9\,h^{-1}`$ Mpc).
- $`h`$ now affects the results, because the photon and neutrino density parameters are fixed by $`T_{\mathrm{cmb}}`$ in physical units ($`\Omega_\gamma \propto h^{-2}`$).
- The comoving distance is the line-of-sight distance $`D_C = \int_0^z dz'/E(z')`$. The angular diameter and luminosity distances use the transverse comoving distance $`D_M`$: $`D_A = D_M/(1+z)`$ and $`D_L = (1+z)D_M`$, with $`D_M = \sinh(\sqrt{\Omega_k}D_C)/\sqrt{\Omega_k}`$ for $`\Omega_k>0`$ and $`\sin(\sqrt{-\Omega_k}D_C)/\sqrt{-\Omega_k}`$ for $`\Omega_k<0`$.
- The `_t` and `_tau` distance functions give the distance today to a source that emitted light at time $`t`$ (conformal time $`\tau`$).
- On error, functions return NaN. Errors include invalid parameters, $`z\le -1`$, $`t\le 0`$, $`k<0`$, and a conformal time that is never reached.
- The universe is assumed to expand monotonically. The contracting phase of a recollapsing closed universe is not handled. Negative $`z`$ (the future) is allowed as long as the expansion continues.

## Functions

| Quantity | as a function of redshift $`z`$ | of time $`t`$ | of conformal time $`\tau`$ |
|---|---|---|---|
| cosmic time $`t`$ | `FLRW_cosmic_time_z` | | |
| conformal time $`\tau`$ | `FLRW_conformal_time_z` | | |
| redshift | | `FLRW_redshift_t` | `FLRW_redshift_tau` |
| scale factor | | `FLRW_scale_factor_t` | `FLRW_scale_factor_tau` |
| comoving distance $`D_C`$ | `FLRW_DC_z` | `FLRW_DC_t` | `FLRW_DC_tau` |
| transverse comoving distance $`D_M`$ | `FLRW_DM_z` | `FLRW_DM_t` | `FLRW_DM_tau` |
| angular diameter distance $`D_A`$ | `FLRW_DA_z` | `FLRW_DA_t` | `FLRW_DA_tau` |
| luminosity distance $`D_L`$ | `FLRW_DL_z` | `FLRW_DL_t` | `FLRW_DL_tau` |
| growth factor, small-scale limit | `FLRW_growth_factor_cb_z` | `FLRW_growth_factor_cb_t` | `FLRW_growth_factor_cb_tau` |
| growth factor, large-scale limit | `FLRW_growth_factor_m_z` | `FLRW_growth_factor_m_t` | `FLRW_growth_factor_m_tau` |
| growth factor (= small-scale limit) | `FLRW_growth_factor_z` | `FLRW_growth_factor_t` | `FLRW_growth_factor_tau` |
| growth factor of $`\delta_{cb}`$ at wavenumber $`k`$ | `FLRW_growth_factor_cb_k_z` | `FLRW_growth_factor_cb_k_t` | `FLRW_growth_factor_cb_k_tau` |
| growth factor of $`\delta_m`$ at wavenumber $`k`$ | `FLRW_growth_factor_m_k_z` | `FLRW_growth_factor_m_k_t` | `FLRW_growth_factor_m_k_tau` |
| critical overdensity $`\delta_c`$ | `FLRW_delta_crit_z` | `FLRW_delta_crit_t` | `FLRW_delta_crit_tau` |
| virial overdensity $`\Delta_{\mathrm{vir}}`$ | `FLRW_virial_Delta_z` | `FLRW_virial_Delta_t` | `FLRW_virial_Delta_tau` |

The functions have the form `double f(const cosm_param* cp, double x)`; the `_k_` functions take `(cp, k, x)` with $`k`$ in units of $`H_0/c`$.
The two limits are normalized to $`D(z{=}0)=1`$; the normalization of the `_k_` functions is described in [Linear growth with massive neutrinos](#linear-growth-with-massive-neutrinos).
Helper functions: `FLRW_omega_m`, `FLRW_omega_cb`, `FLRW_omega_nu`, `FLRW_omega_r`, `FLRW_omega_k` (today), `FLRW_hubble_z` ($`E(z)`$) and `FLRW_omega_m_z` ($`\rho_m/\rho_{\mathrm{crit}}`$ at $`z`$).

```c
#include <stdio.h>
#include "flrw_kit.h"

int main(void) {
  cosm_param cp = {FLRW_PARAM_DEFAULTS,
                   .omega_cdm = 0.26, .omega_b = 0.049, .omega_de = 0.69, .hubble = 0.677};
  FLRW_set_mnu_from_sum(&cp, 0.06, FLRW_NU_NORMAL);
  double z = 1.0, k = 0.1 * FLRW_HUBBLE_DISTANCE_MPC_H; /* 0.1 h/Mpc */
  printf("t(z)     = %g Gyr\n", FLRW_cosmic_time_z(&cp, z) * FLRW_hubble_time_Gyr(&cp));
  printf("DL(z)    = %g Mpc\n", FLRW_DL_z(&cp, z) * FLRW_hubble_distance_Mpc(&cp));
  printf("D_cb(z)  = %g, D_m(z) = %g\n", FLRW_growth_factor_cb_z(&cp, z), FLRW_growth_factor_m_z(&cp, z));
  printf("D(k, z)  = %g\n", FLRW_growth_factor_m_k_z(&cp, k, z));
  printf("dc(z)    = %g, Dvir(z) = %g\n", FLRW_delta_crit_z(&cp, z), FLRW_virial_Delta_z(&cp, z));
  return 0;
}
```

## C++ wrapper (`flrw_kit.hpp`)

`flrw_kit.hpp` is a header-only wrapper (C++11 or later) around the C API. Build the library with the C compiler as usual, include `flrw_kit.hpp`, and link with `-lflrw_kit -lm`.

```cpp
#include <cstdio>
#include "flrw_kit.hpp"

int main() {
  // cdm, b, de, h [, m_nu, T_cmb, N_eff, omega_r]; the optional arguments have the defaults above
  const flrw::Cosmology c(0.26, 0.049, 0.69, 0.677, {{0.0, 0.0086, 0.0502}});
  const auto c2 = flrw::Cosmology::with_mnu_sum(0.26, 0.049, 0.69, 0.677, 0.06, FLRW_NU_NORMAL);
  std::printf("t(z=1)   = %g Gyr\n", c.cosmic_time_z(1.0) * c.hubble_time_Gyr());
  std::printf("D_cb(1)  = %g\n", c2.growth_factor_cb_z(1.0));         // small-scale limit
  std::printf("D_cb(k,1)= %g\n", c2.growth_factor_cb_z(299.8, 1.0));  // k = 0.1 h/Mpc
  std::vector<double> dl = c.evaluate(&flrw::Cosmology::DL_z, {0.5, 1.0, 2.0});  // in c/H0
  return 0;
}
```

```sh
g++ -std=c++11 -I<dir> main.cpp -L<dir> -lflrw_kit -lm
```

- `flrw::Cosmology` takes `(omega_cdm, omega_b, omega_de, hubble, m_nu = {0,0,0}, T_cmb = 2.7255, N_eff = 3.044, omega_r = FLRW_AUTO)` or a `cosm_param`. `Cosmology::with_mnu_sum(omega_cdm, omega_b, omega_de, hubble, sum_eV, ordering, ...)` sets the masses from their sum.
- The member functions have the C names without the `FLRW_` prefix, with the same units and conventions.
- `growth_factor_cb_X` and `growth_factor_m_X` (X = `z`, `t`, `tau`) are overloaded: with one argument they return the small-/large-scale limit, with two arguments `(k, X)` the scale-dependent growth factor (`FLRW_growth_factor_cb_k_X` / `FLRW_growth_factor_m_k_X`).
- The constructor throws `std::invalid_argument` if the C library rejects the parameters ($`\Omega_{cb}\le0`$, $`h\le0`$, $`T_{\mathrm{cmb}}\le0`$, $`N_{\mathrm{eff}}<0`$, a negative mass or a non-finite value); `with_mnu_sum` also throws if the sum is below the minimum.
- Where the C function would return NaN, the member function throws `std::domain_error`. If you prefer NaN, call the C functions with `&c.param()`.
- `evaluate(&flrw::Cosmology::f, xs)` applies a member function to every element of a `std::vector<double>`; `evaluate(&flrw::Cosmology::growth_factor_cb_z, ks, z)` does the same over wavenumbers.
- A `Cosmology` object is immutable, so one object can be shared between threads.

## Linear growth with massive neutrinos

This section explains the three kinds of growth factor, the method used for the scale-dependent one and its normalization.

### Why the growth depends on scale

After becoming non-relativistic, neutrinos still move with a large thermal velocity, $`v_{\mathrm{th}} \simeq 158\,(1+z)\,(1\,\mathrm{eV}/m_\nu)\ \mathrm{km/s}`$. The distance they travel in a Hubble time defines the free-streaming wavenumber

```math
k_{\mathrm{fs}}(a) = \sqrt{\tfrac32}\,\frac{aH(a)}{v_{\mathrm{th}}(a)} \simeq 0.82\,\frac{\sqrt{\Omega_\Lambda + \Omega_m(1+z)^3}}{(1+z)^2}\left(\frac{m_\nu}{1\,\mathrm{eV}}\right) h\,\mathrm{Mpc}^{-1}
```

It plays the role of a Jeans scale: the velocity dispersion acts like a pressure.

- On large scales ($`k\ll k_{\mathrm{fs}}`$) gravity wins and the neutrinos cluster together with cdm and baryons: $`\delta_\nu = \delta_{cb}`$.
- On small scales ($`k\gg k_{\mathrm{fs}}`$) the neutrinos stream out of the perturbations and stay almost uniform: $`\delta_\nu \simeq 0`$.

Since $`k_{\mathrm{fs}}`$ grows with time (as $`a^{1/2}`$ in the matter era), a given scale can move from the second regime to the first.

### Small- and large-scale limits

With $`x=\ln a`$ and $`' = d/dx`$, the density contrast of cdm + baryons obeys

```math
\delta_{cb}'' + \left(2 + \frac{d\ln E}{d\ln a}\right)\delta_{cb}' = \frac32\left[\Omega_{cb}(a)\,\delta_{cb} + \sum_i \tilde\Omega_{\nu,i}(a)\,\delta_{\nu,i}\right],
\qquad \tilde\Omega_{\nu,i}(a) = \frac{\rho_{\nu,i} - 3p_{\nu,i}}{\rho_{\mathrm{crit}}(a)}
```

where $`\tilde\Omega_{\nu,i}`$ is the clustering (non-relativistic) part of each species: relativistic neutrinos do not cluster on sub-horizon scales.

- **Small-scale limit** `FLRW_growth_factor_cb_*` (and `FLRW_growth_factor_*`): $`\delta_\nu = 0`$, so only $`\Omega_{cb}(a)`$ sources the growth while $`H`$ still contains the neutrinos. In the matter era $`\delta_{cb}\propto a^p`$ with $`p = (\sqrt{25-24f_\nu}-1)/4 \simeq 1 - \frac35 f_\nu`$, $`f_\nu = \Omega_\nu/\Omega_m`$: the growth is slower. This limit applies to halo scales and is the one used for $`\delta_c`$.
- **Large-scale limit** `FLRW_growth_factor_m_*`: $`\delta_\nu = \delta_{cb}`$ for the non-relativistic part, so the source is the full $`\Omega_m(a)`$. Without neutrino masses this is the usual growth factor.

Both limits are normalized to $`D(z{=}0)=1`$. Because the small-scale growth is slower, its $`D(z)`$ at high redshift is slightly larger than the large-scale one (for $`\Sigma m_\nu = 0.3`$ eV: $`D_{cb}(1) = 0.6099`$, $`D_m(1) = 0.6057`$).

### Intermediate scales: two-fluid approximation

`FLRW_growth_factor_cb_k_*` and `FLRW_growth_factor_m_k_*` treat each massive species as a fluid with an effective sound speed (cf. Shoji & Komatsu 2010) and integrate it together with cdm + baryons:

```math
\delta_{\nu,i}'' + \left(2 + \frac{d\ln E}{d\ln a}\right)\delta_{\nu,i}' = \frac32\left[\Omega_{cb}\,\delta_{cb} + \sum_j \tilde\Omega_{\nu,j}\,\delta_{\nu,j}\right] - J_i\,\delta_{\nu,i},
\qquad J_i = \frac{c_{s,i}^2 k^2}{(aH)^2}
```

```math
c_{s,i}^2 = \frac{c_{\mathrm{nr}}^2}{1 + 3c_{\mathrm{nr}}^2},\qquad
c_{\mathrm{nr}} = \sqrt{\frac{25\,\zeta(5)}{3\,\zeta(3)}}\,\frac{k_BT_{\nu0}}{m_i\,a} = 2.681\,\frac{k_BT_{\nu0}}{m_i\,a}\quad(\text{in units of } c)
```

- $`c_{\mathrm{nr}}`$ is the non-relativistic sound speed of a Fermi–Dirac gas, $`c_s^2 = \frac59\langle q^4\rangle/\langle q^2\rangle\,(T/m)^2`$; the interpolation tends to $`c/\sqrt3`$ for relativistic neutrinos.
- $`J_i`$ is the Jeans term: $`J_i\gg1`$ means the scale is below the free-streaming scale. When $`J_i > 400`$ the species is treated quasi-statically ($`\delta_{\nu,i}''`$ and $`\delta_{\nu,i}'`$ dropped, so $`\delta_{\nu,i}`$ follows from the right-hand side); it is integrated again when $`J_i < 100`$. This avoids following fast acoustic oscillations, which real neutrinos damp by phase mixing anyway.
- `FLRW_growth_factor_cb_k_*` returns the growth of $`\delta_{cb}`$, and `FLRW_growth_factor_m_k_*` the growth of the total matter contrast $`\delta_m = (\rho_{cb}\delta_{cb} + \sum_i(\rho_{\nu,i}-3p_{\nu,i})\delta_{\nu,i})/\rho_m`$. The total matter power spectrum is $`P_m\propto D_m(k)^2`$; $`P_{cb}\propto D_{cb}(k)^2`$ is the one relevant for halos (see below).
- $`k\to0`$ gives the large-scale limit for both; $`k\to\infty`$ gives the small-scale limit for $`\delta_{cb}`$ (up to the normalization below). Massless neutrinos never cluster, so with massless neutrinos only there is no scale dependence.

The fluid approximation describes the neutrino perturbation itself only approximately. Because the neutrinos are a small fraction $`f_\nu`$ of the matter, the error in $`\delta_{cb}`$ and $`\delta_m`$ is much smaller than the suppression. For precision work, compare with a Boltzmann code (CLASS, CAMB).

### Normalization of the scale-dependent growth factor

All wavenumbers start from the same initial amplitude, and the result is divided by the large-scale limit today:

```math
D_{cb}(k,z) = \frac{\delta_{cb}(k,z)}{\delta_m(k{\to}0,\,z{=}0)},\qquad
D_m(k,z) = \frac{\delta_m(k,z)}{\delta_m(k{\to}0,\,z{=}0)}
```

so $`D_m(k{\to}0, 0) = 1`$ and, on scales where neutrinos free-stream, $`D(k,0) < 1`$. The suppression caused by the neutrino masses is kept in the value itself. Examples ($`\Omega_{cb}=0.309`$, $`\Omega_\Lambda=0.69`$, $`h=0.677`$):

| $`k`$ [$`h\,\mathrm{Mpc}^{-1}`$] | $`D_{cb}(k,0)`$, 0.06 eV (normal) | $`D_m(k,0)`$, 0.06 eV | $`D_{cb}(k,0)`$, 0.3 eV (degenerate) | $`D_m(k,0)`$, 0.3 eV |
|---|---|---|---|---|
| 0.001 | 0.9996 | 0.9994 | 0.9992 | 0.9991 |
| 0.01 | 0.9918 | 0.9898 | 0.9650 | 0.9610 |
| 0.1 | 0.9882 | 0.9840 | 0.9324 | 0.9146 |
| 1 | 0.9881 | 0.9837 | 0.9307 | 0.9100 |

$`1 - D_m(k,0)^2`$ on small scales is about $`7f_\nu`$–$`8f_\nu`$, the familiar suppression of the matter power spectrum by massive neutrinos.

How to use it:

- **Evolving a spectrum in redshift**: $`P(k,z) = P(k,0)\,[D(k,z)/D(k,0)]^2`$. The ratio does not depend on the normalization.
- **Suppression relative to clustering neutrinos**: $`D_m(k,0)^2`$ is the ratio of the matter power spectrum to that of a model in which the neutrinos cluster like cdm with the same $`\Omega_m`$ and the same early-time amplitude.
- Since $`D_{cb}(k\to\infty,z)/D_{cb}(k\to\infty,0) = D_{cb}(z)`$, the small-scale limit normalized at $`z=0`$ is recovered from the `_k_` function.

### Assumptions and limitations

- The equations are Newtonian and valid on sub-horizon scales ($`k\gg aH`$); for larger scales the growth factor is a convention.
- cdm and baryons are one pressureless fluid; the coupling of baryons to photons before recombination is ignored.
- Integration starts in the radiation era ($`a_i = 10^{-9}`$) on the Meszaros growing mode $`\delta\propto a + \frac23a_{\mathrm{eq}}`$ (radiation and relativistic neutrinos smooth). Dynamical neutrinos start with $`\delta_\nu=\delta_{cb}`$ (adiabatic).

## Spherical collapse ($`\delta_c`$, $`\Delta_{\mathrm{vir}}`$)

A top-hat sphere of cdm + baryons collapses in a background with smooth photons, neutrinos and $`\Lambda`$ (neutrinos do not cluster on halo scales). With $`y = R/(a\,r_L)`$ ($`r_L`$: Lagrangian radius, $`1+\delta = y^{-3}`$), the radius obeys

```math
y'' + \left(2 + \frac{d\ln E}{d\ln a}\right)y' = -\frac12\,\Omega_{cb}(a)\left(y^{-2} - y\right)
```

The smooth components enter only through $`H(a)`$. The equation is integrated for $`w = y-1`$, which keeps the precision of the small initial contrast, starting on the linear growing mode in the radiation era.

**1.** The initial amplitude is adjusted (Brent's method) so that $`y\to0`$ at $`a_c`$. Near collapse $`y\propto(x_c - x)^{2/3}`$, so the collapse time is extrapolated from $`y = 10^{-4}`$ as $`x_c = x - \frac23 y/y'`$.

**2.** $`\delta_c`$ is that initial amplitude evolved with the linear equation in the small-scale limit (`FLRW_growth_factor_cb`) up to $`a_c`$.

**3.** The virialization time is the collapse time. Turnaround is where $`(ay)'=0`$. The virial radius $`R_{\mathrm{vir}} = y_vR_{\mathrm{ta}}`$ follows from energy conservation between turnaround and collapse and the virial theorem at collapse. A smooth component with density $`\rho_s`$ and pressure $`p_s`$ adds the potential energy $`\propto(\rho_s+3p_s)R^2`$ (for $`\Lambda`$, $`\rho+3p = -2\rho_\Lambda`$), evaluated at each epoch:

```math
4\eta_c\,y_v^3 - 2(1+\eta_{\mathrm{ta}})\,y_v + 1 = 0,\qquad
\eta(a) = -\frac{R_{\mathrm{ta}}^3}{2\,\Omega_{cb}}\sum_s\frac{\rho_s(a) + 3p_s(a)}{\rho_{c0}}
```

With $`\Lambda`$ alone, $`\eta_c=\eta_{\mathrm{ta}}=\Omega_\Lambda R_{\mathrm{ta}}^3/\Omega_{cb}`$, the classical result. With massive neutrinos and radiation, energy is not exactly conserved because their densities change with time; evaluating their contribution at each epoch is an approximation, and its effect on $`\Delta_{\mathrm{vir}}`$ is small because these components are a small fraction of the density at late times.

### Which mean density $`\Delta_{\mathrm{vir}}`$ refers to

$`\Delta_{\mathrm{vir}}`$ is the mean **matter** density within the virial radius divided by the mean **matter** density of the universe at collapse:

```math
\Delta_{\mathrm{vir}} = \frac{\bar\rho_{cb,\mathrm{halo}} + \bar\rho_{\nu,m}(a_c)}{\bar\rho_{cb}(a_c) + \bar\rho_{\nu,m}(a_c)},\qquad
\bar\rho_{\nu,m} = \sum_i(\rho_{\nu,i} - 3p_{\nu,i})
```

- The halo contains the collapsed cdm + baryons and the smooth neutrinos at their mean density.
- The reference density includes the non-relativistic neutrinos, i.e. $`\rho_m`$ as defined above (`FLRW_omega_m_z`).
- For another reference density, rescale:
  - critical density: $`\Delta_{\mathrm{c}} = \Delta_{\mathrm{vir}}\,\Omega_m(z)`$ with $`\Omega_m(z)`$ = `FLRW_omega_m_z(cp, z)`;
  - cdm + baryon density: $`\Delta_{cb} = \Delta_{\mathrm{vir}}\,\Omega_m(z)/\Omega_{cb}(z)`$ with $`\Omega_{cb}(z) = \Omega_{cb}(1+z)^3/E(z)^2`$.
- At fixed $`\Omega_{cb}`$, adding neutrino mass raises $`\rho_m`$ and therefore lowers $`\Delta_{\mathrm{vir}}`$ (e.g. 333.7 → 325.7 at $`z=0`$ for $`\Sigma m_\nu = 0`$ → 0.3 eV with $`\Omega_{cb}=0.309`$).

For Einstein–de Sitter, $`\delta_c = \frac35(3\pi/2)^{2/3} = 1.68647`$ and $`\Delta_{\mathrm{vir}} = 18\pi^2 = 177.65`$.

## Methods

### Time and distances

With $`R(a) = a^4E(a)^2`$ and $`u=\sqrt a`$, which keeps the integrands smooth at $`a\to0`$ with or without radiation,

```math
t(a) = \int_0^{\sqrt a}\frac{2u^3\,du}{\sqrt{R(u^2)}},\qquad
\tau(a) = \int_0^{\sqrt a}\frac{2u\,du}{\sqrt{R(u^2)}}
```

evaluated by globally adaptive Gauss–Kronrod quadrature (7–15 points, relative accuracy $`10^{-13}`$). The inverse maps from $`t`$ and $`\tau`$ to $`a`$ use a bracketed Newton iteration in $`\ln a`$.

### Neutrino energy density and pressure

For one species, $`\rho_\nu`$ and $`p_\nu`$ in units of the massless energy density are

```math
F(y) = \frac{120}{7\pi^4}\int_0^\infty\frac{q^2\sqrt{q^2+y^2}}{e^q+1}\,dq,\qquad
P(y) = \frac{120}{7\pi^4}\int_0^\infty\frac{q^4}{3\sqrt{q^2+y^2}\,(e^q+1)}\,dq,\qquad
y = \frac{m_\nu a}{k_BT_{\nu0}}
```

$`\ln F`$ and $`\ln P`$ are tabulated with their first and second derivatives in $`\ln y`$ for $`10^{-3}\le y\le 10^4`$ (spacing 0.05) and interpolated by quintic Hermite polynomials (relative error ~$`10^{-12}`$). Outside the table the small-$`y`$ series $`F = 1 + 5y^2/7\pi^2`$, $`P = 1/3 - 5y^2/21\pi^2`$ and the large-$`y`$ expansion are used.

### Growth and collapse equations

The ordinary differential equations are integrated with an adaptive Dormand–Prince 5(4) method (relative tolerance $`10^{-11}`$) in $`x=\ln a`$.

## Tests

`test_flrw_kit.c` checks the C library:

- Classic model (no radiation, no neutrinos): analytic solutions for Einstein–de Sitter, open $`\Lambda=0`$ (Mattig's formula, $`t_0`$) and flat ΛCDM ($`t(a)`$), relative error $`\lesssim10^{-12}`$. $`\delta_c`$ and $`\Delta_{\mathrm{vir}}`$ agree with the closed-form spherical collapse solution to $`10^{-8}`$, the growth factor with a direct integration to $`10^{-8}`$, $`\delta_c`$ with an independent integration of the background and shell equations in time to a few $`10^{-5}`$, and $`\Delta_{\mathrm{vir}}`$ with the Bryan & Norman (1998) fit (~1%).
- Parameters: defaults, $`\Omega_\gamma\propto T^4/h^2`$, masses from the sum and ordering (sum and splittings), rejection below the minimum sum, $`\Omega_\nu h^2 = \Sigma m_\nu/93.04`$ eV.
- Radiation + matter: analytic $`t(a)`$, $`\tau(a)`$ and the Meszaros solution $`D\propto a + \frac23a_{\mathrm{eq}}`$ ($`10^{-9}`$); $`\delta_c`$ against an independent integration ($`10^{-5}`$).
- Massive neutrinos: $`E(z)`$ against direct Fermi–Dirac integrals ($`10^{-11}`$); both growth limits against an independent integration with these integrals ($`10^{-6}`$); $`k\to0`$ and $`k\to\infty`$ limits of the scale-dependent growth; monotonic decrease with $`k`$; small-scale suppression of 6–10 $`f_\nu`$; the matter-era exponent $`p`$; $`\delta_c`$ and $`\Delta_{\mathrm{vir}}`$ close to the massless case.
- Round trips through the inverse functions, consistency between the $`z`$, $`t`$ and $`\tau`$ versions, and NaN for invalid input.

`test_flrw_kit_cpp.cpp` checks the C++ wrapper: every member function, including the `(k, X)` overloads, returns exactly the value of the C function; the constructor defaults and `with_mnu_sum` agree with the C functions; `evaluate()` works over arguments and wavenumbers; the expected exceptions are thrown.

## Cost

Each call performs numerical integration, and nothing is cached. Typical times are 0.05 ms for time and distances, 0.2 ms for the inverse maps from $`t`$ and $`\tau`$, 0.5 ms for a growth-factor limit, 3 ms for a scale-dependent growth factor, and 7 ms for $`\delta_c`$ or $`\Delta_{\mathrm{vir}}`$. If you need many evaluations, build a table on your side and interpolate; the functions can be called in parallel.
