#ifndef FLRW_KIT_H
#define FLRW_KIT_H
/* FLRW_kit: 一様等方宇宙 (FLRW) の背景量・距離・線形成長・球対称崩壊
 *
 * 宇宙モデル
 *   E(a)^2 = (H/H0)^2 = Ωm a^-3 + Ωk a^-2 + ΩΛ
 *     Ωm = omega_cdm + omega_b + omega_nu   (ニュートリノは非相対論的物質として扱う)
 *     ΩΛ = omega_de                         (宇宙定数, w = -1)
 *     Ωk = 1 - Ωm - ΩΛ                      (曲率は自動的に決まる)
 *   放射は含まない.
 *
 * 単位
 *   時刻 t, 共形時刻 tau は 1/H0, 距離は c/H0 で規格化する (c = 1 なので tau と共動距離は同じ単位).
 *   t, tau はビッグバン (a = 0) から測る. 現在 (z = 0) は t = t0, tau = tau0.
 *
 * 規約
 *   - 共動距離 DC は視線方向の共動距離 ∫ dz/E(z). 角径距離・光度距離は横方向共動距離
 *     DM (曲率で DC から変換) を使い, DA = DM/(1+z), DL = (1+z) DM.
 *   - t, tau を引数とする距離は「時刻 t (共形時刻 tau) に光を出した天体を現在観測したとき」の距離.
 *   - 成長因子は D(z=0) = 1 に規格化した線形成長モード.
 *   - 臨界密度揺らぎ δc, ビリアル過密度 Δvir はトップハット球対称崩壊 (ΩΛ と曲率を含む厳密解) から求める.
 *     ビリアル化の時刻 = 崩壊 (半径 → 0) の時刻とし, ビリアル半径はエネルギー保存とビリアル定理
 *     (宇宙定数のポテンシャルを含む) から決める. Δvir は崩壊時の宇宙平均 (物質) 密度に対する比.
 *   - 宇宙は単調に膨張しているものとする (再収縮する閉じた宇宙の収縮期は扱わない).
 *
 * エラー時 (不正なパラメータ, 定義域外の引数, 根が見つからない等) は NaN を返す.
 * 内部状態を持たないのでスレッドセーフ (OpenMP の並列領域から呼んでよい).
 */

typedef struct cosm_param {
  double omega_cdm; /* ダークマターの密度パラメータ */
  double omega_b;   /* バリオンの密度パラメータ */
  double omega_de;  /* ダークエネルギー (宇宙定数) の密度パラメータ */
  double hubble;    /* h = H0 / (100 km/s/Mpc) */
  double omega_nu;  /* ニュートリノの密度パラメータ */
} cosm_param;

#ifdef __cplusplus
extern "C" {
#endif

/* ---- 補助的な量 ---- */
double FLRW_omega_m(const cosm_param* cp);              /* Ωm = Ωcdm + Ωb + Ων */
double FLRW_omega_k(const cosm_param* cp);              /* Ωk = 1 - Ωm - ΩΛ */
double FLRW_hubble_z(const cosm_param* cp, double z);   /* E(z) = H(z)/H0 */
double FLRW_omega_m_z(const cosm_param* cp, double z);  /* Ωm(z) = Ωm (1+z)^3 / E(z)^2 */
double FLRW_hubble_time_Gyr(const cosm_param* cp);      /* 1/H0 [Gyr] */
double FLRW_hubble_distance_Mpc(const cosm_param* cp);  /* c/H0 [Mpc] */

/* ---- 時刻 [1/H0] ---- */
double FLRW_cosmic_time_z(const cosm_param* cp, double z);
double FLRW_conformal_time_z(const cosm_param* cp, double z);

double FLRW_redshift_t(const cosm_param* cp, double t);
double FLRW_scale_factor_t(const cosm_param* cp, double t);
double FLRW_redshift_tau(const cosm_param* cp, double tau);
double FLRW_scale_factor_tau(const cosm_param* cp, double tau);

/* ---- 距離 [c/H0] ---- */
double FLRW_DC_z(const cosm_param* cp, double z);  /* 共動距離 (視線方向) */
double FLRW_DM_z(const cosm_param* cp, double z);  /* 横方向共動距離 */
double FLRW_DA_z(const cosm_param* cp, double z);  /* 角径距離 */
double FLRW_DL_z(const cosm_param* cp, double z);  /* 光度距離 */

double FLRW_DC_t(const cosm_param* cp, double t);
double FLRW_DM_t(const cosm_param* cp, double t);
double FLRW_DA_t(const cosm_param* cp, double t);
double FLRW_DL_t(const cosm_param* cp, double t);

double FLRW_DC_tau(const cosm_param* cp, double tau);
double FLRW_DM_tau(const cosm_param* cp, double tau);
double FLRW_DA_tau(const cosm_param* cp, double tau);
double FLRW_DL_tau(const cosm_param* cp, double tau);

/* ---- 線形成長因子 (D(z=0) = 1) ---- */
double FLRW_growth_factor_z(const cosm_param* cp, double z);
double FLRW_growth_factor_t(const cosm_param* cp, double t);
double FLRW_growth_factor_tau(const cosm_param* cp, double tau);

/* ---- 球対称崩壊: 臨界密度揺らぎ δc (ビリアル化の時点に線形外挿した値) ---- */
double FLRW_delta_crit_z(const cosm_param* cp, double z);
double FLRW_delta_crit_t(const cosm_param* cp, double t);
double FLRW_delta_crit_tau(const cosm_param* cp, double tau);

/* ---- 球対称崩壊: ビリアル半径内の平均密度 / 宇宙平均密度 ---- */
double FLRW_virial_Delta_z(const cosm_param* cp, double z);
double FLRW_virial_Delta_t(const cosm_param* cp, double t);
double FLRW_virial_Delta_tau(const cosm_param* cp, double tau);

#ifdef __cplusplus
}
#endif

#endif /* FLRW_KIT_H */
