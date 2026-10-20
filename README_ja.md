# FLRW_kit — FLRW 宇宙の背景量・距離・線形成長・球対称崩壊を計算する C ライブラリ

[English](README.md)

与えられた宇宙論パラメータに対して、宇宙年齢・共形時刻・各種距離・線形成長因子・臨界密度揺らぎ $`\delta_c`$・ビリアル過密度 $`\Delta_{\mathrm{vir}}`$ を計算する。
C11 + libm だけで動き、内部状態を持たないので OpenMP の並列領域からそのまま呼べる。

## ビルド

```sh
cd FLRW_kit
make            # libflrw_kit.a
make test       # テスト
```

利用側では `#include "flrw_kit.h"` して `-L<dir> -lflrw_kit -lm` をリンクする。

## 宇宙モデル

```c
typedef struct cosm_param {
  double omega_cdm; /* ダークマター */
  double omega_b;   /* バリオン */
  double omega_de;  /* ダークエネルギー (宇宙定数) */
  double hubble;    /* h  (H0 = 100 h km/s/Mpc) */
  double omega_nu;  /* ニュートリノ */
} cosm_param;
```

```math
E(a)^2 \equiv \left(\frac{H}{H_0}\right)^2 = \Omega_m a^{-3} + \Omega_k a^{-2} + \Omega_\Lambda
```

- $`\Omega_m = \Omega_{\mathrm{cdm}} + \Omega_b + \Omega_\nu`$（ニュートリノは非相対論的物質として扱う）
- $`\Omega_\Lambda`$ = `omega_de`（$`w=-1`$）
- $`\Omega_k = 1 - \Omega_m - \Omega_\Lambda`$（曲率は自動的に決まる）
- 放射は含まない。

## 単位と規約

- 時刻 $`t`$、共形時刻 $`\tau`$ は $`1/H_0`$、距離は $`c/H_0`$ で規格化（$`c=1`$ なので $`\tau`$ と共動距離は同じ単位）。$`t,\tau`$ はビッグバン（$`a=0`$）から測る。
- $`h`$ は規格化された量には現れない。物理単位への換算に `FLRW_hubble_time_Gyr()`（$`1/H_0`$ [Gyr]）と `FLRW_hubble_distance_Mpc()`（$`c/H_0`$ [Mpc]）を使う。
- 共動距離 $`D_C = \int_0^z dz'/E(z')`$（視線方向）。角径距離・光度距離は横方向共動距離 $`D_M`$ を使い $`D_A = D_M/(1+z)`$、$`D_L = (1+z)D_M`$。
  - $`D_M = \sinh(\sqrt{\Omega_k}D_C)/\sqrt{\Omega_k}`$（$`\Omega_k>0`$）、$`\sin(\sqrt{-\Omega_k}D_C)/\sqrt{-\Omega_k}`$（$`\Omega_k<0`$）
- `_t`, `_tau` 版の距離は「時刻 $`t`$（共形時刻 $`\tau`$）に光を出した天体を現在観測したときの距離」。$`D_C = \tau_0 - \tau`$。
- エラー（不正なパラメータ、$`z\le -1`$、$`t\le 0`$、到達できない共形時刻など）は NaN を返す。
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
| 線形成長因子 $`D`$（$`D(z{=}0)=1`$） | `FLRW_growth_factor_z` | `FLRW_growth_factor_t` | `FLRW_growth_factor_tau` |
| 臨界密度揺らぎ $`\delta_c`$ | `FLRW_delta_crit_z` | `FLRW_delta_crit_t` | `FLRW_delta_crit_tau` |
| ビリアル過密度 $`\Delta_{\mathrm{vir}}`$ | `FLRW_virial_Delta_z` | `FLRW_virial_Delta_t` | `FLRW_virial_Delta_tau` |

すべて `double f(const cosm_param* cp, double x)` の形。補助関数として `FLRW_omega_m`, `FLRW_omega_k`, `FLRW_hubble_z`（$`E(z)`$）, `FLRW_omega_m_z`（$`\Omega_m(z)`$）もある。

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

## 計算方法

### 時刻・距離

$`u=\sqrt a`$ と変数変換して $`a\to0`$ の特異性を除き、大域適応型 Gauss–Kronrod（7–15 点、相対精度 $`10^{-13}`$）で積分する。

```math
t(a) = \int_0^{\sqrt a}\frac{2u^2\,du}{\sqrt{\Omega_m + \Omega_k u^2 + \Omega_\Lambda u^6}},\qquad
\tau(a) = \int_0^{\sqrt a}\frac{2\,du}{\sqrt{\Omega_m + \Omega_k u^2 + \Omega_\Lambda u^6}}
```

$`t`$、$`\tau`$ からの逆変換は $`\ln a`$ についての括弧付き Newton 法（$`dt/d\ln a = 1/E`$、$`d\tau/d\ln a = 1/(aE)`$）。

### 線形成長因子

物質＋宇宙定数＋曲率で厳密な Heath (1977) の積分表示を使う。

```math
D_1(a) = \frac{5\Omega_m}{2}E(a)\int_0^a\frac{da'}{[a'E(a')]^3}\quad(a\to0\ \text{で}\ D_1\to a),\qquad
D(z) = D_1(a)/D_1(1)
```

ニュートリノも物質として一緒に成長するとみなしているので、自由流減衰によるスケール依存性は含まない。

### 球対称崩壊（$`\delta_c`$, $`\Delta_{\mathrm{vir}}`$）

フィッティング公式ではなく、$`\Lambda`$ と曲率を含むトップハット崩壊を厳密に解く。
Lagrange 半径を 1 とする単位（$`GM=\Omega_m/2`$）でシェル半径 $`R`$ は

```math
\ddot R = -\frac{\Omega_m}{2R^2} + \Omega_\Lambda R,\qquad
\dot R^2 = \frac{\Omega_m}{R} + \Omega_\Lambda R^2 + K,\quad K = -\frac{\Omega_m}{R_{\mathrm{ta}}} - \Omega_\Lambda R_{\mathrm{ta}}^2
```

背景（$`R=a`$）では $`K=\Omega_k`$ なので、揺らぎは「曲率」の差 $`\delta K = \Omega_k - K`$ で表される。

**1.** 成長モードのみ（ビッグバン時刻が背景と同じ）なら、崩壊時刻はターンアラウンド時刻の 2 倍。

```math
t_{\mathrm{ta}}(R_{\mathrm{ta}}) = \int_0^{\pi/2}\frac{2R_{\mathrm{ta}}^{3/2}\sin^2\theta\,d\theta}{\sqrt{\Omega_m - \Omega_\Lambda R_{\mathrm{ta}}^3\sin^2\theta(1+\sin^2\theta)}}
```

$`2t_{\mathrm{ta}}(R_{\mathrm{ta}}) = t(a_c)`$ を Brent 法で解いて $`R_{\mathrm{ta}}`$ を求める（$`\Omega_\Lambda>0`$ では $`R_{\mathrm{ta}}^3<\Omega_m/2\Omega_\Lambda`$）。

**2.** 初期の線形揺らぎは $`\delta = \frac35(\delta K/\Omega_m)\,a`$ なので、崩壊時刻に線形外挿した値は

```math
\delta_c = \frac35\,\frac{\delta K}{\Omega_m}\,D_1(a_c)
```

**3.** ビリアル化の時刻は崩壊時刻とする。ビリアル半径 $`R_{\mathrm{vir}}=yR_{\mathrm{ta}}`$ は、エネルギー保存とビリアル定理（$`U_\Lambda\propto -R^2`$ も含めて $`2K=-U_G+2U_\Lambda`$）から

```math
4\eta y^3 - 2(1+\eta)y + 1 = 0,\qquad \eta = \frac{\Omega_\Lambda R_{\mathrm{ta}}^3}{\Omega_m}
```

宇宙の平均物質密度に対する比は $`\Delta_{\mathrm{vir}} = (a_c/R_{\mathrm{vir}})^3`$。臨界密度に対する比が欲しいときは $`\Delta_{\mathrm{vir}}\,\Omega_m(z)`$。

EdS では $`\delta_c = \frac35(3\pi/2)^{2/3} = 1.68647`$、$`\Delta_{\mathrm{vir}} = 18\pi^2 = 177.65`$ に一致する。

## テスト（`test_flrw_kit.c`）

- EdS・開いた $`\Lambda=0`$（Mattig の式、$`t_0`$）・平坦 ΛCDM（$`t(a)`$）の解析解: 相対誤差 $`\lesssim10^{-12}`$
- 成長因子: 成長方程式を RK4 で直接積分した結果と比較
- $`\delta_c`$: 背景とシェルの運動方程式を直接積分し、二分法で崩壊時刻を合わせた結果と比較（曲率あり・なし、$`\Omega_m>1`$ も含む）
- $`\Delta_{\mathrm{vir}}`$: Bryan & Norman (1998) のフィット（精度 ~1%）と比較
- 逆関数の往復、$`z`$ 版と $`t,\tau`$ 版の整合性、定義域外で NaN

## 計算コスト

1 回の呼び出しで数値積分を行う（キャッシュはしない）。目安は時刻・距離・成長因子が数 µs、$`t,\tau`$ からの逆変換が数十 µs、$`\delta_c`$・$`\Delta_{\mathrm{vir}}`$ が数百 µs。大量に呼ぶ場合は利用側で表を作って補間するとよい。
