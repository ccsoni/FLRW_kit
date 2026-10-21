# FLRW_kit — FLRW 宇宙の背景量・距離・線形成長・球対称崩壊を計算する C ライブラリ

[English](README.md)

与えられた宇宙論パラメータに対して、宇宙年齢・共形時刻・各種距離・線形成長因子（小スケール極限、大スケール極限、波数の関数）・臨界密度揺らぎ $`\delta_c`$・ビリアル過密度 $`\Delta_{\mathrm{vir}}`$ を計算する。
放射と質量のあるニュートリノを含む。ニュートリノは種ごとにフェルミ・ディラック分布に従って、放射から物質へと移り変わる。
C11 + libm だけで動き、内部状態を持たないので OpenMP の並列領域からそのまま呼べる。

## ビルド

```sh
make            # libflrw_kit.a
make test       # C と C++ のテスト（片方だけなら make test_c / make test_cpp）
```

利用側では `#include "flrw_kit.h"` して `-L<dir> -lflrw_kit -lm` をリンクする。
C++ からは `flrw_kit.h` をそのまま include してもよい（`extern "C"` で囲んである）し、後述の C++ ラッパー `flrw_kit.hpp` を使ってもよい。

`flrw_nu_table.h` はニュートリノのエネルギー密度と圧力の表を生成したファイルで、リポジトリに含めてあるのでビルドに numpy は要らない。`make table` で `tools/gen_nu_table.py` から作り直せる（numpy の入った python3 が必要）。

## 宇宙モデル

```c
typedef struct cosm_param {
  double omega_cdm; /* 冷たいダークマター */
  double omega_b;   /* バリオン */
  double omega_de;  /* ダークエネルギー（宇宙定数） */
  double hubble;    /* h  (H0 = 100 h km/s/Mpc) */
  double m_nu[3];   /* ニュートリノの質量固有値 [eV]（質量のない種は 0） */
  double T_cmb;     /* 現在の CMB 温度 [K]              デフォルト 2.7255 */
  double N_eff;     /* ニュートリノの有効種数           デフォルト 3.044 */
  double omega_r;   /* 光子の密度パラメータ             デフォルト FLRW_AUTO（T_cmb と h から計算） */
} cosm_param;
```

```math
E(a)^2 \equiv \left(\frac{H}{H_0}\right)^2 = \Omega_{cb}\,a^{-3} + \Omega_\gamma\,a^{-4} + \frac{\rho_\nu(a)}{\rho_{c0}} + \Omega_k\,a^{-2} + \Omega_\Lambda
```

- $`\Omega_{cb} = \Omega_{\mathrm{cdm}} + \Omega_b`$。$`\Omega_\Lambda`$ = `omega_de`（$`w=-1`$）。
- 光子: `omega_r = FLRW_AUTO` のとき $`\Omega_\gamma h^2 = 2.47298\times10^{-5}\,(T_{\mathrm{cmb}}/2.7255\,\mathrm{K})^4`$。0 以上の値を与えればその値を使う。
- ニュートリノ: 質量 `m_nu[i]` の 3 種で、温度 $`T_\nu = (4/11)^{1/3}(N_{\mathrm{eff}}/3)^{1/4}\,T_{\mathrm{cmb}}`$ のフェルミ・ディラック分布に従う。$`m_i a \ll 3k_BT_{\nu0}`$ では放射のように、その後は物質のように振る舞う。質量のない種はずっと放射のまま。$`N_{\mathrm{eff}}=0`$ でニュートリノを除ける。
- $`\Omega_k = 1 - \Omega_{cb} - \Omega_\gamma - \Omega_\nu - \Omega_\Lambda`$（曲率は自動的に決まる）。平坦にしたいときは、`omega_de` を `1 - FLRW_omega_cb - FLRW_omega_r - FLRW_omega_nu` のように選ぶ。
- 「物質」（$`\Omega_m`$、`FLRW_omega_m_z`、$`\Delta_{\mathrm{vir}}`$）は $`\rho_m = \rho_{cb} + \sum_i(\rho_{\nu,i} - 3p_{\nu,i})`$、つまり cdm + バリオン + ニュートリノの非相対論的な部分を指す。非相対論的になった種は全量を数え、相対論的な種や質量のない種は数えない。

### パラメータの初期化

C の構造体にはデフォルト値を書けないので、まず `FLRW_param_init()` か初期化子マクロ `FLRW_PARAM_DEFAULTS` で初期化する。どちらも `T_cmb = 2.7255`、`N_eff = 3.044`、`omega_r = FLRW_AUTO` とし、ほかはすべて 0（質量のないニュートリノ）にする。

```c
cosm_param cp;
FLRW_param_init(&cp);
cp.omega_cdm = 0.26; cp.omega_b = 0.049; cp.omega_de = 0.69; cp.hubble = 0.677;

/* C99 の指示付き初期化子を使う場合 */
cosm_param cp2 = {FLRW_PARAM_DEFAULTS,
                  .omega_cdm = 0.26, .omega_b = 0.049, .omega_de = 0.69, .hubble = 0.677};
```

### 質量の和と質量階層から質量固有値を決める

`FLRW_set_mnu_from_sum(&cp, sum_eV, ordering)` は、和が `sum_eV` になるように `m_nu[0..2]` = $`(m_1, m_2, m_3)`$ を決める。振動実験の質量二乗差 $`\Delta m^2_{21} = 7.42\times10^{-5}\,\mathrm{eV}^2`$、$`\Delta m^2_{31} = 2.515\times10^{-3}\,\mathrm{eV}^2`$（順階層）、$`|\Delta m^2_{32}| = 2.498\times10^{-3}\,\mathrm{eV}^2`$（逆階層）を使う。これらの値はコンパイル時に `-DFLRW_DM21_SQ=...`、`-DFLRW_DM31_SQ_NORMAL=...`、`-DFLRW_DM32_SQ_INVERTED=...` で変えられる。

| `ordering` | 質量 | 和の最小値 |
|---|---|---|
| `FLRW_NU_NORMAL` | $`m_1 < m_2 < m_3`$ | 0.0588 eV |
| `FLRW_NU_INVERTED` | $`m_3 < m_1 < m_2`$ | 0.0992 eV |
| `FLRW_NU_DEGENERATE` | $`m_1 = m_2 = m_3 = \Sigma/3`$ | 0 |

戻り値は 0。和がその階層の最小値より小さいときは −1 を返す。$`N_{\mathrm{eff}}=3.044`$ では、非相対論的なニュートリノの密度は $`\Omega_\nu h^2 = \Sigma m_\nu / 93.04\,\mathrm{eV}`$ になる。

## 単位と規約

- 時刻 $`t`$、共形時刻 $`\tau`$ は $`1/H_0`$、距離は $`c/H_0`$、波数 $`k`$ は $`H_0/c`$ を単位とする。光は共形時刻 $`d\tau`$ の間に共動距離 $`c\,d\tau`$ だけ進むので、この単位では共形時刻 $`\tau`$ に光を出した天体までの共動距離が数値的に $`D_C = \tau_0 - \tau`$ になる。$`t,\tau`$ はビッグバン（$`a=0`$）から測る。
- `FLRW_HUBBLE_DISTANCE_MPC_H` = 2997.92458 は $`h^{-1}`$ Mpc 単位の $`c/H_0`$。距離に掛けると $`h^{-1}`$ Mpc に、$`h\,\mathrm{Mpc}^{-1}`$ 単位の波数に掛けると $`H_0/c`$ 単位の $`k`$ になる（例: $`k = 0.1\,h\,\mathrm{Mpc}^{-1}`$ は $`k = 299.8`$）。この換算は $`h`$ によらない。
- 物理単位への換算は、時刻に `FLRW_hubble_time_Gyr()`（$`1/H_0 = 9.778\,h^{-1}`$ Gyr）を、距離に `FLRW_hubble_distance_Mpc()`（$`c/H_0 = 2997.9\,h^{-1}`$ Mpc）を掛ける。
- 光子とニュートリノの密度パラメータは物理単位の $`T_{\mathrm{cmb}}`$ で決まる（$`\Omega_\gamma \propto h^{-2}`$）ので、$`h`$ は結果に影響する。
- 共動距離は視線方向の $`D_C = \int_0^z dz'/E(z')`$。角径距離・光度距離は横方向共動距離 $`D_M`$ を使って $`D_A = D_M/(1+z)`$、$`D_L = (1+z)D_M`$。$`D_M`$ は $`\Omega_k>0`$ で $`\sinh(\sqrt{\Omega_k}D_C)/\sqrt{\Omega_k}`$、$`\Omega_k<0`$ で $`\sin(\sqrt{-\Omega_k}D_C)/\sqrt{-\Omega_k}`$。
- `_t`, `_tau` 版の距離は、時刻 $`t`$（共形時刻 $`\tau`$）に光を出した天体を現在観測したときの距離。
- エラー（不正なパラメータ、$`z\le -1`$、$`t\le 0`$、$`k<0`$、到達できない共形時刻など）は NaN を返す。
- 単調に膨張する宇宙を仮定する（再収縮する閉じた宇宙の収縮期は扱わない）。$`z<0`$（未来）も膨張が続く範囲で計算できる。

## 関数一覧

| 量 | 赤方偏移 $`z`$ の関数 | 時刻 $`t`$ の関数 | 共形時刻 $`\tau`$ の関数 |
|---|---|---|---|
| 宇宙年齢 $`t`$ | `FLRW_cosmic_time_z` | | |
| 共形時刻 $`\tau`$ | `FLRW_conformal_time_z` | | |
| 赤方偏移 | | `FLRW_redshift_t` | `FLRW_redshift_tau` |
| スケール因子 | | `FLRW_scale_factor_t` | `FLRW_scale_factor_tau` |
| 共動距離 $`D_C`$ | `FLRW_DC_z` | `FLRW_DC_t` | `FLRW_DC_tau` |
| 横方向共動距離 $`D_M`$ | `FLRW_DM_z` | `FLRW_DM_t` | `FLRW_DM_tau` |
| 角径距離 $`D_A`$ | `FLRW_DA_z` | `FLRW_DA_t` | `FLRW_DA_tau` |
| 光度距離 $`D_L`$ | `FLRW_DL_z` | `FLRW_DL_t` | `FLRW_DL_tau` |
| 線形成長因子・小スケール極限 | `FLRW_growth_factor_cb_z` | `FLRW_growth_factor_cb_t` | `FLRW_growth_factor_cb_tau` |
| 線形成長因子・大スケール極限 | `FLRW_growth_factor_m_z` | `FLRW_growth_factor_m_t` | `FLRW_growth_factor_m_tau` |
| 線形成長因子（= 小スケール極限） | `FLRW_growth_factor_z` | `FLRW_growth_factor_t` | `FLRW_growth_factor_tau` |
| 波数 $`k`$ での $`\delta_{cb}`$ の成長因子 | `FLRW_growth_factor_cb_k_z` | `FLRW_growth_factor_cb_k_t` | `FLRW_growth_factor_cb_k_tau` |
| 波数 $`k`$ での $`\delta_m`$ の成長因子 | `FLRW_growth_factor_m_k_z` | `FLRW_growth_factor_m_k_t` | `FLRW_growth_factor_m_k_tau` |
| 臨界密度揺らぎ $`\delta_c`$ | `FLRW_delta_crit_z` | `FLRW_delta_crit_t` | `FLRW_delta_crit_tau` |
| ビリアル過密度 $`\Delta_{\mathrm{vir}}`$ | `FLRW_virial_Delta_z` | `FLRW_virial_Delta_t` | `FLRW_virial_Delta_tau` |

関数は `double f(const cosm_param* cp, double x)` の形で、`_k_` の付く関数だけは `(cp, k, x)` を取る（$`k`$ は $`H_0/c`$ 単位）。
2 つの極限は $`D(z{=}0)=1`$ に規格化してある。`_k_` の関数の規格化は「[質量のあるニュートリノを含む線形成長](#質量のあるニュートリノを含む線形成長)」で説明する。
補助関数として `FLRW_omega_m`, `FLRW_omega_cb`, `FLRW_omega_nu`, `FLRW_omega_r`, `FLRW_omega_k`（現在の値）、`FLRW_hubble_z`（$`E(z)`$）、`FLRW_omega_m_z`（赤方偏移 $`z`$ での $`\rho_m/\rho_{\mathrm{crit}}`$）がある。

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

## C++ ラッパー（`flrw_kit.hpp`）

`flrw_kit.hpp` は C の API を包むヘッダだけの C++ ラッパー（C++11 以降）。ライブラリはこれまでどおり C コンパイラでビルドし、`flrw_kit.hpp` を include して `-lflrw_kit -lm` をリンクする。

```cpp
#include <cstdio>
#include "flrw_kit.hpp"

int main() {
  // cdm, b, de, h [, m_nu, T_cmb, N_eff, omega_r]。省略した引数は上のデフォルト値になる
  const flrw::Cosmology c(0.26, 0.049, 0.69, 0.677, {{0.0, 0.0086, 0.0502}});
  const auto c2 = flrw::Cosmology::with_mnu_sum(0.26, 0.049, 0.69, 0.677, 0.06, FLRW_NU_NORMAL);
  std::printf("t(z=1)   = %g Gyr\n", c.cosmic_time_z(1.0) * c.hubble_time_Gyr());
  std::printf("D_cb(1)  = %g\n", c2.growth_factor_cb_z(1.0));         // 小スケール極限
  std::printf("D_cb(k,1)= %g\n", c2.growth_factor_cb_z(299.8, 1.0));  // k = 0.1 h/Mpc
  std::vector<double> dl = c.evaluate(&flrw::Cosmology::DL_z, {0.5, 1.0, 2.0});  // c/H0 単位
  return 0;
}
```

```sh
g++ -std=c++11 -I<dir> main.cpp -L<dir> -lflrw_kit -lm
```

- `flrw::Cosmology` は `(omega_cdm, omega_b, omega_de, hubble, m_nu = {0,0,0}, T_cmb = 2.7255, N_eff = 3.044, omega_r = FLRW_AUTO)` か `cosm_param` から作る。`Cosmology::with_mnu_sum(omega_cdm, omega_b, omega_de, hubble, sum_eV, ordering, ...)` は質量の和から質量固有値を決める。
- メンバ関数名は C の関数名から `FLRW_` を除いたもので、単位と規約も同じ。
- `growth_factor_cb_X` と `growth_factor_m_X`（X = `z`, `t`, `tau`）はオーバーロードしてある。引数が 1 つなら小スケール極限・大スケール極限を、`(k, X)` の 2 つなら波数依存の成長因子（`FLRW_growth_factor_cb_k_X` / `FLRW_growth_factor_m_k_X`）を返す。
- C のライブラリがパラメータを受け付けない場合（$`\Omega_{cb}\le0`$、$`h\le0`$、$`T_{\mathrm{cmb}}\le0`$、$`N_{\mathrm{eff}}<0`$、負の質量、有限でない値）、コンストラクタは `std::invalid_argument` を投げる。`with_mnu_sum` は和が最小値より小さいときも投げる。
- C の関数が NaN を返す場合、メンバ関数は `std::domain_error` を投げる。NaN のほうがよければ `&c.param()` を渡して C の関数を呼ぶ。
- `evaluate(&flrw::Cosmology::f, xs)` は `std::vector<double>` の各要素にメンバ関数を適用する。`evaluate(&flrw::Cosmology::growth_factor_cb_z, ks, z)` は同じことを波数の列について行う。
- `Cosmology` オブジェクトは変更できないので、1 つのオブジェクトを複数のスレッドで共有してよい。

## 質量のあるニュートリノを含む線形成長

ここでは 3 種類の成長因子、波数依存の成長因子の計算方法と規格化を説明する。

### 成長がスケールに依存する理由

ニュートリノは非相対論的になった後も大きな熱速度 $`v_{\mathrm{th}} \simeq 158\,(1+z)\,(1\,\mathrm{eV}/m_\nu)\ \mathrm{km/s}`$ を持つ。1 ハッブル時間に進む距離から、自由流（free-streaming）波数が決まる。

```math
k_{\mathrm{fs}}(a) = \sqrt{\tfrac32}\,\frac{aH(a)}{v_{\mathrm{th}}(a)} \simeq 0.82\,\frac{\sqrt{\Omega_\Lambda + \Omega_m(1+z)^3}}{(1+z)^2}\left(\frac{m_\nu}{1\,\mathrm{eV}}\right) h\,\mathrm{Mpc}^{-1}
```

これは Jeans スケールにあたり、速度分散が圧力の役割をする。

- 大スケール（$`k\ll k_{\mathrm{fs}}`$）では重力が勝ち、ニュートリノも cdm・バリオンと一緒に集まる: $`\delta_\nu = \delta_{cb}`$。
- 小スケール（$`k\gg k_{\mathrm{fs}}`$）ではニュートリノが揺らぎから逃げ出し、ほぼ一様なまま: $`\delta_\nu \simeq 0`$。

$`k_{\mathrm{fs}}`$ は時間とともに大きくなる（物質優勢期には $`a^{1/2}`$ に比例）ので、あるスケールが後者から前者に移ることがある。

### 小スケール極限と大スケール極限

$`x=\ln a`$、$`' = d/dx`$ として、cdm + バリオンの密度揺らぎは次の式に従う。

```math
\delta_{cb}'' + \left(2 + \frac{d\ln E}{d\ln a}\right)\delta_{cb}' = \frac32\left[\Omega_{cb}(a)\,\delta_{cb} + \sum_i \tilde\Omega_{\nu,i}(a)\,\delta_{\nu,i}\right],
\qquad \tilde\Omega_{\nu,i}(a) = \frac{\rho_{\nu,i} - 3p_{\nu,i}}{\rho_{\mathrm{crit}}(a)}
```

$`\tilde\Omega_{\nu,i}`$ は各種の集まりうる（非相対論的な）部分で、相対論的なニュートリノはホライズン内では集まらない。

- **小スケール極限** `FLRW_growth_factor_cb_*`（と `FLRW_growth_factor_*`）: $`\delta_\nu = 0`$。成長の源は $`\Omega_{cb}(a)`$ だけになるが、$`H`$ にはニュートリノが含まれたままなので成長が遅くなる。物質優勢期には $`f_\nu = \Omega_\nu/\Omega_m`$ として $`\delta_{cb}\propto a^p`$、$`p = (\sqrt{25-24f_\nu}-1)/4 \simeq 1 - \frac35 f_\nu`$。ハローのスケールにはこちらが当てはまり、$`\delta_c`$ もこちらを使う。
- **大スケール極限** `FLRW_growth_factor_m_*`: 非相対論的な部分について $`\delta_\nu = \delta_{cb}`$ とし、成長の源は $`\Omega_m(a)`$ 全体。ニュートリノに質量がなければ通常の成長因子と同じ。

どちらの極限も $`D(z{=}0)=1`$ に規格化してある。小スケール極限のほうが成長が遅いので、高赤方偏移での $`D(z)`$ は大スケール極限より少し大きくなる（$`\Sigma m_\nu = 0.3`$ eV では $`D_{cb}(1) = 0.6099`$、$`D_m(1) = 0.6057`$）。

### 中間スケール: 2 流体近似

`FLRW_growth_factor_cb_k_*` と `FLRW_growth_factor_m_k_*` は、質量のある各種を有効音速を持つ流体として扱い（Shoji & Komatsu 2010 と同様の方法）、cdm + バリオンと連立させて積分する。

```math
\delta_{\nu,i}'' + \left(2 + \frac{d\ln E}{d\ln a}\right)\delta_{\nu,i}' = \frac32\left[\Omega_{cb}\,\delta_{cb} + \sum_j \tilde\Omega_{\nu,j}\,\delta_{\nu,j}\right] - J_i\,\delta_{\nu,i},
\qquad J_i = \frac{c_{s,i}^2 k^2}{(aH)^2}
```

```math
c_{s,i}^2 = \frac{c_{\mathrm{nr}}^2}{1 + 3c_{\mathrm{nr}}^2},\qquad
c_{\mathrm{nr}} = \sqrt{\frac{25\,\zeta(5)}{3\,\zeta(3)}}\,\frac{k_BT_{\nu0}}{m_i\,a} = 2.681\,\frac{k_BT_{\nu0}}{m_i\,a}\quad(c\ \text{単位})
```

- $`c_{\mathrm{nr}}`$ はフェルミ・ディラック気体の非相対論的な音速 $`c_s^2 = \frac59\langle q^4\rangle/\langle q^2\rangle\,(T/m)^2`$。上の内挿式は、相対論的なニュートリノで $`c/\sqrt3`$ に近づく。
- $`J_i`$ は Jeans 項で、$`J_i\gg1`$ はそのスケールが自由流スケールより小さいことを意味する。$`J_i > 400`$ になった種は準静的に扱い（$`\delta_{\nu,i}''`$ と $`\delta_{\nu,i}'`$ を落として、右辺から $`\delta_{\nu,i}`$ を決める）、$`J_i < 100`$ になったら再び積分する。速い音響振動を追いかけずに済むようにするためで、実際のニュートリノではこの振動は位相混合で減衰する。
- `FLRW_growth_factor_cb_k_*` は $`\delta_{cb}`$ の成長を、`FLRW_growth_factor_m_k_*` は全物質の揺らぎ $`\delta_m = (\rho_{cb}\delta_{cb} + \sum_i(\rho_{\nu,i}-3p_{\nu,i})\delta_{\nu,i})/\rho_m`$ の成長を返す。全物質のパワースペクトルは $`P_m\propto D_m(k)^2`$。ハローに関係するのは $`P_{cb}\propto D_{cb}(k)^2`$ のほう（後述）。
- $`k\to0`$ ではどちらも大スケール極限に、$`k\to\infty`$ では $`\delta_{cb}`$ が小スケール極限（下の規格化を除いて）になる。質量のないニュートリノは決して集まらないので、質量のないニュートリノしかなければスケール依存性はない。

流体近似はニュートリノの揺らぎ自体を近似的にしか表さない。ただしニュートリノは物質の小さな割合 $`f_\nu`$ なので、$`\delta_{cb}`$ と $`\delta_m`$ への誤差は抑制の大きさよりずっと小さい。精密な計算では Boltzmann コード（CLASS、CAMB）と比べること。

### 波数依存の成長因子の規格化

すべての波数を同じ初期振幅から始め、結果を現在の大スケール極限で割る。

```math
D_{cb}(k,z) = \frac{\delta_{cb}(k,z)}{\delta_m(k{\to}0,\,z{=}0)},\qquad
D_m(k,z) = \frac{\delta_m(k,z)}{\delta_m(k{\to}0,\,z{=}0)}
```

したがって $`D_m(k{\to}0, 0) = 1`$ で、ニュートリノが自由流で広がるスケールでは $`D(k,0) < 1`$ になる。ニュートリノの質量による抑制が、値そのものに残る。例（$`\Omega_{cb}=0.309`$、$`\Omega_\Lambda=0.69`$、$`h=0.677`$）:

| $`k`$ [$`h\,\mathrm{Mpc}^{-1}`$] | $`D_{cb}(k,0)`$、0.06 eV（順階層） | $`D_m(k,0)`$、0.06 eV | $`D_{cb}(k,0)`$、0.3 eV（縮退） | $`D_m(k,0)`$、0.3 eV |
|---|---|---|---|---|
| 0.001 | 0.9996 | 0.9994 | 0.9992 | 0.9991 |
| 0.01 | 0.9918 | 0.9898 | 0.9650 | 0.9610 |
| 0.1 | 0.9882 | 0.9840 | 0.9324 | 0.9146 |
| 1 | 0.9881 | 0.9837 | 0.9307 | 0.9100 |

小スケールでの $`1 - D_m(k,0)^2`$ は $`7f_\nu`$〜$`8f_\nu`$ 程度で、質量のあるニュートリノによる物質パワースペクトルのよく知られた抑制にあたる。

使い方:

- **スペクトルを赤方偏移方向に進める**: $`P(k,z) = P(k,0)\,[D(k,z)/D(k,0)]^2`$。この比は規格化によらない。
- **ニュートリノが集まる場合と比べた抑制**: $`D_m(k,0)^2`$ は、同じ $`\Omega_m`$ と同じ初期振幅で、ニュートリノが cdm と同じように集まるモデルに対する物質パワースペクトルの比。
- $`D_{cb}(k\to\infty,z)/D_{cb}(k\to\infty,0) = D_{cb}(z)`$ なので、`_k_` の関数からも $`z=0`$ で規格化した小スケール極限が得られる。

### 仮定と適用範囲

- 方程式はニュートン近似で、ホライズン内（$`k\gg aH`$）で成り立つ。それより大きなスケールでの成長因子は便宜上の定義になる。
- cdm とバリオンは 1 つの圧力のない流体として扱う。再結合前のバリオンと光子の結合は無視する。
- 積分は放射優勢期（$`a_i = 10^{-9}`$）に、Meszaros の成長モード $`\delta\propto a + \frac23a_{\mathrm{eq}}`$（放射と相対論的ニュートリノは一様）から始める。積分するニュートリノは $`\delta_\nu=\delta_{cb}`$（断熱的）から始める。

## 球対称崩壊（$`\delta_c`$, $`\Delta_{\mathrm{vir}}`$）

cdm + バリオンのトップハット球が、一様な光子・ニュートリノ・$`\Lambda`$ の背景の中で崩壊するとする（ハローのスケールではニュートリノは集まらない）。$`y = R/(a\,r_L)`$（$`r_L`$: Lagrange 半径、$`1+\delta = y^{-3}`$）として、半径は次の式に従う。

```math
y'' + \left(2 + \frac{d\ln E}{d\ln a}\right)y' = -\frac12\,\Omega_{cb}(a)\left(y^{-2} - y\right)
```

一様な成分は $`H(a)`$ を通してだけ入る。小さな初期揺らぎの精度を保つため $`w = y-1`$ について積分し、放射優勢期の線形成長モードから始める。

**1.** $`a_c`$ で $`y\to0`$ となるように初期振幅を調整する（Brent 法）。崩壊の直前では $`y\propto(x_c - x)^{2/3}`$ なので、崩壊時刻は $`y = 10^{-4}`$ の点から $`x_c = x - \frac23 y/y'`$ で外挿する。

**2.** $`\delta_c`$ は、その初期振幅を小スケール極限の線形方程式（`FLRW_growth_factor_cb`）で $`a_c`$ まで進めた値。

**3.** ビリアル化の時刻は崩壊時刻とする。ターンアラウンドは $`(ay)'=0`$ となる点。ビリアル半径 $`R_{\mathrm{vir}} = y_vR_{\mathrm{ta}}`$ は、ターンアラウンドから崩壊までのエネルギー保存と、崩壊時のビリアル定理から決める。密度 $`\rho_s`$、圧力 $`p_s`$ の一様な成分は $`(\rho_s+3p_s)R^2`$ に比例するポテンシャルエネルギーを加え（$`\Lambda`$ では $`\rho+3p = -2\rho_\Lambda`$）、これをそれぞれの時刻で評価する。

```math
4\eta_c\,y_v^3 - 2(1+\eta_{\mathrm{ta}})\,y_v + 1 = 0,\qquad
\eta(a) = -\frac{R_{\mathrm{ta}}^3}{2\,\Omega_{cb}}\sum_s\frac{\rho_s(a) + 3p_s(a)}{\rho_{c0}}
```

$`\Lambda`$ だけなら $`\eta_c=\eta_{\mathrm{ta}}=\Omega_\Lambda R_{\mathrm{ta}}^3/\Omega_{cb}`$ で、古典的な結果に一致する。質量のあるニュートリノや放射があると、それらの密度が時間変化するのでエネルギーは厳密には保存しない。各時刻で評価するのは近似だが、これらの成分は後期には密度のごく一部なので、$`\Delta_{\mathrm{vir}}`$ への影響は小さい。

### $`\Delta_{\mathrm{vir}}`$ の基準となる平均密度

$`\Delta_{\mathrm{vir}}`$ は、ビリアル半径内の平均**物質**密度を、崩壊時の宇宙の平均**物質**密度で割ったもの。

```math
\Delta_{\mathrm{vir}} = \frac{\bar\rho_{cb,\mathrm{halo}} + \bar\rho_{\nu,m}(a_c)}{\bar\rho_{cb}(a_c) + \bar\rho_{\nu,m}(a_c)},\qquad
\bar\rho_{\nu,m} = \sum_i(\rho_{\nu,i} - 3p_{\nu,i})
```

- ハローの中身は、崩壊した cdm + バリオンと、平均密度のままの一様なニュートリノ。
- 基準の密度には非相対論的なニュートリノを含む。つまり上で定義した $`\rho_m`$（`FLRW_omega_m_z`）。
- 別の基準密度にしたいときは換算する。
  - 臨界密度: $`\Delta_{\mathrm{c}} = \Delta_{\mathrm{vir}}\,\Omega_m(z)`$、$`\Omega_m(z)`$ = `FLRW_omega_m_z(cp, z)`。
  - cdm + バリオンの密度: $`\Delta_{cb} = \Delta_{\mathrm{vir}}\,\Omega_m(z)/\Omega_{cb}(z)`$、$`\Omega_{cb}(z) = \Omega_{cb}(1+z)^3/E(z)^2`$。
- $`\Omega_{cb}`$ を固定してニュートリノの質量を足すと $`\rho_m`$ が増えるので、$`\Delta_{\mathrm{vir}}`$ は小さくなる（例: $`\Omega_{cb}=0.309`$、$`z=0`$ で $`\Sigma m_\nu = 0`$ → 0.3 eV のとき 333.7 → 325.7）。

EdS では $`\delta_c = \frac35(3\pi/2)^{2/3} = 1.68647`$、$`\Delta_{\mathrm{vir}} = 18\pi^2 = 177.65`$ になる。

## 計算方法

### 時刻・距離

$`R(a) = a^4E(a)^2`$ とし、放射の有無によらず $`a\to0`$ で被積分関数がなめらかになるように $`u=\sqrt a`$ と変数変換する。

```math
t(a) = \int_0^{\sqrt a}\frac{2u^3\,du}{\sqrt{R(u^2)}},\qquad
\tau(a) = \int_0^{\sqrt a}\frac{2u\,du}{\sqrt{R(u^2)}}
```

これを大域適応型 Gauss–Kronrod（7–15 点、相対精度 $`10^{-13}`$）で積分する。$`t`$、$`\tau`$ からの逆変換は $`\ln a`$ についての括弧付き Newton 法。

### ニュートリノのエネルギー密度と圧力

1 種の $`\rho_\nu`$ と $`p_\nu`$ を、質量のない場合のエネルギー密度を単位として表すと次のとおり。

```math
F(y) = \frac{120}{7\pi^4}\int_0^\infty\frac{q^2\sqrt{q^2+y^2}}{e^q+1}\,dq,\qquad
P(y) = \frac{120}{7\pi^4}\int_0^\infty\frac{q^4}{3\sqrt{q^2+y^2}\,(e^q+1)}\,dq,\qquad
y = \frac{m_\nu a}{k_BT_{\nu0}}
```

$`\ln F`$ と $`\ln P`$ を、$`\ln y`$ についての 1 階・2 階微分とともに $`10^{-3}\le y\le 10^4`$ で表にし（間隔 0.05）、5 次 Hermite 多項式で補間する（相対誤差 ~$`10^{-12}`$）。表の外では、小さな $`y`$ の級数 $`F = 1 + 5y^2/7\pi^2`$、$`P = 1/3 - 5y^2/21\pi^2`$ と、大きな $`y`$ での展開を使う。

### 成長と崩壊の方程式

常微分方程式は、$`x=\ln a`$ について適応型 Dormand–Prince 5(4) 法（相対許容誤差 $`10^{-11}`$）で積分する。

## テスト

`test_flrw_kit.c` は C のライブラリを確かめる。

- 古典的なモデル（放射・ニュートリノなし）: EdS・開いた $`\Lambda=0`$（Mattig の式、$`t_0`$）・平坦 ΛCDM（$`t(a)`$）の解析解と相対誤差 $`\lesssim10^{-12}`$ で一致。$`\delta_c`$ と $`\Delta_{\mathrm{vir}}`$ は球対称崩壊の厳密解と $`10^{-8}`$、成長因子は直接積分と $`10^{-8}`$、$`\delta_c`$ は背景とシェルの運動方程式を独立に時間積分した結果と数 $`10^{-5}`$、$`\Delta_{\mathrm{vir}}`$ は Bryan & Norman (1998) のフィットと ~1% で一致。
- パラメータ: デフォルト値、$`\Omega_\gamma\propto T^4/h^2`$、和と階層から決めた質量（和と質量二乗差）、最小値より小さい和を受け付けないこと、$`\Omega_\nu h^2 = \Sigma m_\nu/93.04`$ eV。
- 放射 + 物質: $`t(a)`$、$`\tau(a)`$ の解析解と Meszaros の解 $`D\propto a + \frac23a_{\mathrm{eq}}`$（$`10^{-9}`$）。$`\delta_c`$ を独立な積分と比較（$`10^{-5}`$）。
- 質量のあるニュートリノ: $`E(z)`$ をフェルミ・ディラック積分の直接計算と比較（$`10^{-11}`$）。2 つの成長の極限を、この積分を使った独立な積分と比較（$`10^{-6}`$）。波数依存の成長の $`k\to0`$ と $`k\to\infty`$ の極限、$`k`$ について単調に減ること、小スケールでの抑制が $`6`$〜$`10\,f_\nu`$、物質優勢期の指数 $`p`$、$`\delta_c`$ と $`\Delta_{\mathrm{vir}}`$ が質量のない場合に近いこと。
- 逆関数の往復、$`z`$ 版と $`t,\tau`$ 版の整合性、不正な入力で NaN。

`test_flrw_kit_cpp.cpp` は C++ ラッパーを確かめる。`(k, X)` のオーバーロードを含むすべてのメンバ関数が C の関数とまったく同じ値を返すこと、コンストラクタのデフォルト値と `with_mnu_sum` が C の関数と一致すること、`evaluate()` が引数の列と波数の列の両方で動くこと、想定どおりに例外が投げられることを調べる。

## 計算コスト

1 回の呼び出しで数値積分を行う（キャッシュはしない）。目安は時刻・距離が 0.05 ms、$`t,\tau`$ からの逆変換が 0.2 ms、成長因子の極限が 0.5 ms、波数依存の成長因子が 3 ms、$`\delta_c`$・$`\Delta_{\mathrm{vir}}`$ が 7 ms。大量に呼ぶ場合は利用側で表を作って補間するとよい（関数は並列に呼べる）。
