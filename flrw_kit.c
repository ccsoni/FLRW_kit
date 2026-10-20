/* FLRW_kit: 実装
 *
 * 積分はすべて u = sqrt(a) で変数変換して a -> 0 の特異性を除いてから
 * 適応型 Gauss-Kronrod (7-15 点) で行う. E(a)^2 a^3 = Ωm + Ωk a + ΩΛ a^3 として
 *
 *   t(a)   = ∫_0^a da/(a E)   = ∫_0^√a 2u^2 du / sqrt(Ωm + Ωk u^2 + ΩΛ u^6)
 *   tau(a) = ∫_0^a da/(a^2 E) = ∫_0^√a 2    du / sqrt(Ωm + Ωk u^2 + ΩΛ u^6)
 *   DC(a)  = ∫_a^1 da/(a^2 E) = ∫_√a^1  (tau と同じ被積分関数)
 *   D1(a)  = (5Ωm/2) E(a) ∫_0^a da/(a E)^3
 *          = (5Ωm/2) E(a) ∫_0^√a 2u^4 du / (Ωm + Ωk u^2 + ΩΛ u^6)^{3/2}
 *            (Heath 1977. 物質 + 宇宙定数 + 曲率で厳密. a -> 0 で D1 -> a)
 *
 * 球対称崩壊 (質量 M のシェル, 半径 R. Lagrange 半径 r_L = 1 になるよう質量を規格化し,
 * GM = Ωm/2 [H0 = 1, a(t0) = 1]):
 *   d^2R/dt^2 = -Ωm/(2R^2) + ΩΛ R
 *   エネルギー: (dR/dt)^2 = Ωm/R + ΩΛ R^2 + K,  K = -Ωm/R_ta - ΩΛ R_ta^2
 *   背景 (R = a) では K = Ωk なので, 摂動は「曲率」の差 δK = Ωk - K で表される.
 *   成長モードのみ (ビッグバン時刻が背景と同じ) なら崩壊時刻 = 2 × ターンアラウンド時刻:
 *     t_ta(R_ta) = ∫_0^R_ta dR / sqrt(Ωm (1/R - 1/R_ta) + ΩΛ (R^2 - R_ta^2))
 *                = ∫_0^{π/2} 2 R_ta^{3/2} s^2 dθ / sqrt(Ωm - ΩΛ R_ta^3 s^2 (1+s^2)),  s = sinθ
 *   初期の線形密度揺らぎは δ = (3/5)(δK/Ωm) a なので, 崩壊時に線形外挿した値は
 *     δc = (3/5) (δK/Ωm) D1(a_c)
 *   ビリアル半径 R_vir = y R_ta はエネルギー保存とビリアル定理
 *   (U_G ∝ -1/R, U_Λ ∝ -R^2 なので 2K = -U_G + 2U_Λ) から
 *     4η y^3 - 2(1+η) y + 1 = 0,   η = ΩΛ R_ta^3 / Ωm
 *   ビリアル半径内の平均密度 / 宇宙の平均物質密度 = (a_c / R_vir)^3.
 */
#include "flrw_kit.h"

#include <float.h>
#include <math.h>

#define INT_EPSREL 1e-13
#define INT_MAXSEG 400
#define ROOT_XTOL 1e-14
#define ROOT_MAXIT 200
#define A_MAX 1e12 /* 逆関数探索の a の上限 */
#define PI 3.14159265358979323846

typedef struct {
  double om, ok, ol; /* Ωm, Ωk, ΩΛ */
} Bg;

typedef double (*IntFn)(double x, const void* par);

/* ------------------------------------------------------------------ */
/* 適応型 Gauss-Kronrod 7-15 点積分                                     */
/* ------------------------------------------------------------------ */
static const double xgk[8] = {
    0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
    0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
    0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
    0.207784955007898467600689403773245, 0.000000000000000000000000000000000};
static const double wgk[8] = {
    0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
    0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
    0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
    0.204432940075298892414161999234649, 0.209482141084727828012999174891714};
static const double wg[4] = {
    0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
    0.381830050505118944950369775488975, 0.417959183673469387755102040816327};

static double gk15(IntFn f, const void* par, double a, double b, double* err) {
  const double c = 0.5 * (a + b), h = 0.5 * (b - a);
  const double fc = f(c, par);
  double rk = fc * wgk[7], rg = fc * wg[3];
  for (int j = 0; j < 7; ++j) {
    const double dx = h * xgk[j];
    const double fs = f(c - dx, par) + f(c + dx, par);
    rk += wgk[j] * fs;
    if (j & 1) rg += wg[j / 2] * fs;
  }
  *err = fabs((rk - rg) * h);
  return rk * h;
}

/* 大域適応型: 誤差最大の区間を二分することを, 誤差の総和が許容値を下回るか
 * 丸め誤差の限界に達するまで繰り返す (QUADPACK の QAG と同じ方針) */
static double integrate(IntFn f, const void* par, double a, double b) {
  if (a == b) return 0.0;
  double lo[INT_MAXSEG], hi[INT_MAXSEG], val[INT_MAXSEG], err[INT_MAXSEG];
  int n = 1;
  lo[0] = a;
  hi[0] = b;
  val[0] = gk15(f, par, a, b, &err[0]);
  double sum = val[0], esum = err[0];
  while (isfinite(sum) && n < INT_MAXSEG) {
    if (esum <= fmax(INT_EPSREL * fabs(sum), 1e3 * DBL_MIN)) break;
    int k = 0;
    for (int i = 1; i < n; ++i)
      if (err[i] > err[k]) k = i;
    if (err[k] <= 50.0 * DBL_EPSILON * fabs(val[k])) break; /* 丸め誤差の限界 */
    const double m = 0.5 * (lo[k] + hi[k]);
    if (!(m > fmin(lo[k], hi[k]) && m < fmax(lo[k], hi[k]))) break;
    double el, er;
    const double vl = gk15(f, par, lo[k], m, &el);
    const double vr = gk15(f, par, m, hi[k], &er);
    sum += vl + vr - val[k];
    esum += el + er - err[k];
    lo[n] = m, hi[n] = hi[k], val[n] = vr, err[n] = er;
    hi[k] = m, val[k] = vl, err[k] = el;
    ++n;
  }
  /* 和の丸め誤差の蓄積を避けるため最後に足し直す */
  sum = 0.0;
  for (int i = 0; i < n; ++i) sum += val[i];
  return sum;
}

/* ------------------------------------------------------------------ */
/* Brent 法 (f(lo), f(hi) は異符号であること)                            */
/* ------------------------------------------------------------------ */
typedef double (*RootFn)(double x, const void* par);

static double brent(RootFn f, const void* par, double a, double b, double fa, double fb) {
  if (!(fa * fb <= 0.0)) return NAN;
  double c = b, fc = fb, d = 0.0, e = 0.0;
  for (int it = 0; it < ROOT_MAXIT; ++it) {
    if ((fb > 0.0) == (fc > 0.0)) {
      c = a;
      fc = fa;
      e = d = b - a;
    }
    if (fabs(fc) < fabs(fb)) {
      a = b, b = c, c = a;
      fa = fb, fb = fc, fc = fa;
    }
    const double tol = 2.0 * DBL_EPSILON * fabs(b) + 0.5 * ROOT_XTOL * fabs(c);
    const double xm = 0.5 * (c - b);
    if (fabs(xm) <= tol || fb == 0.0) return b;
    if (fabs(e) >= tol && fabs(fa) > fabs(fb)) {
      double p, q, r;
      const double s = fb / fa;
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0) q = -q;
      p = fabs(p);
      if (2.0 * p < fmin(3.0 * xm * q - fabs(tol * q), fabs(e * q))) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    b += (fabs(d) > tol) ? d : copysign(tol, xm);
    fb = f(b, par);
    if (!isfinite(fb)) return NAN;
  }
  return b;
}

/* ------------------------------------------------------------------ */
/* 背景                                                                */
/* ------------------------------------------------------------------ */
static int make_bg(const cosm_param* cp, Bg* b) {
  if (!cp) return -1;
  b->om = cp->omega_cdm + cp->omega_b + cp->omega_nu;
  b->ol = cp->omega_de;
  b->ok = 1.0 - b->om - b->ol;
  if (!(b->om > 0.0) || !isfinite(b->om) || !isfinite(b->ol)) return -1;
  return 0;
}

/* a^3 E(a)^2 = Ωm + Ωk a + ΩΛ a^3 を u = sqrt(a) で書いたもの */
static inline double a3E2_u(const Bg* b, double u) {
  const double u2 = u * u;
  return b->om + u2 * (b->ok + b->ol * u2 * u2);
}

static inline double E_a(const Bg* b, double a) {
  return sqrt(b->om / (a * a * a) + b->ok / (a * a) + b->ol);
}

static double f_time(double u, const void* p) {
  return 2.0 * u * u / sqrt(a3E2_u((const Bg*)p, u));
}
static double f_conf(double u, const void* p) {
  return 2.0 / sqrt(a3E2_u((const Bg*)p, u));
}
static double f_grow(double u, const void* p) {
  const double q = a3E2_u((const Bg*)p, u);
  const double u2 = u * u;
  return 2.0 * u2 * u2 / (q * sqrt(q));
}

static double time_a(const Bg* b, double a) { return integrate(f_time, b, 0.0, sqrt(a)); }
static double conf_a(const Bg* b, double a) { return integrate(f_conf, b, 0.0, sqrt(a)); }
static double DC_a(const Bg* b, double a) { return integrate(f_conf, b, sqrt(a), 1.0); }

/* a -> 0 で D1 -> a となる規格化の成長因子 */
static double D1_a(const Bg* b, double a) {
  return 2.5 * b->om * E_a(b, a) * integrate(f_grow, b, 0.0, sqrt(a));
}

static double DM_from_DC(const Bg* b, double dc) {
  const double x = b->ok * dc * dc;
  if (fabs(x) < 1e-6) return dc * (1.0 + x / 6.0 * (1.0 + x / 20.0));
  const double sk = sqrt(fabs(b->ok));
  return (b->ok > 0.0) ? sinh(sk * dc) / sk : sin(sk * dc) / sk;
}

static inline double a_of_z(double z) { return (z > -1.0) ? 1.0 / (1.0 + z) : NAN; }

/* ------------------------------------------------------------------ */
/* t(a), tau(a) の逆関数 (x = ln a で safeguarded Newton)               */
/* ------------------------------------------------------------------ */
typedef struct {
  const Bg* b;
  int conformal;
  double target;
} InvPar;

static double inv_f(double x, const void* p) {
  const InvPar* ip = p;
  const double a = exp(x);
  return (ip->conformal ? conf_a(ip->b, a) : time_a(ip->b, a)) - ip->target;
}

/* d t/d ln a = 1/E,  d tau/d ln a = 1/(a E) */
static double inv_df(double x, const InvPar* ip) {
  const double a = exp(x);
  const double E = E_a(ip->b, a);
  return ip->conformal ? 1.0 / (a * E) : 1.0 / E;
}

static double invert(const Bg* b, int conformal, double target) {
  if (!(target > 0.0) || !isfinite(target)) return NAN;
  const InvPar ip = {b, conformal, target};
  /* 物質優勢の解 t = (2/3) a^{3/2}/√Ωm, tau = 2 √a/√Ωm を初期値にする */
  const double so = sqrt(b->om);
  double x = conformal ? 2.0 * log(0.5 * so * target) : (2.0 / 3.0) * log(1.5 * so * target);
  x = fmin(x, log(A_MAX));
  double fx = inv_f(x, &ip);
  if (!isfinite(fx)) return NAN;

  /* 括弧付け */
  double lo, hi, flo, fhi;
  if (fx < 0.0) {
    lo = x, flo = fx;
    hi = x, fhi = fx;
    while (fhi < 0.0) {
      lo = hi, flo = fhi;
      hi += 1.0;
      if (hi > log(A_MAX)) return NAN; /* 共形時刻の上限 (ΩΛ > 0) 等 */
      fhi = inv_f(hi, &ip);
      if (!isfinite(fhi)) return NAN;
    }
  } else {
    hi = x, fhi = fx;
    lo = x, flo = fx;
    while (flo > 0.0) {
      hi = lo, fhi = flo;
      lo -= 1.0;
      flo = inv_f(lo, &ip);
      if (!isfinite(flo)) return NAN;
    }
  }
  if (flo == 0.0) return exp(lo);
  if (fhi == 0.0) return exp(hi);

  /* Newton (括弧の外に出る・収束が遅いときは二分法) */
  x = 0.5 * (lo + hi);
  fx = inv_f(x, &ip);
  for (int it = 0; it < ROOT_MAXIT; ++it) {
    if (fx < 0.0) lo = x;
    else hi = x;
    double xn = x - fx / inv_df(x, &ip);
    if (!(xn > lo && xn < hi)) xn = 0.5 * (lo + hi);
    const double dx = xn - x;
    x = xn;
    if (fabs(dx) < 4.0 * DBL_EPSILON * fmax(1.0, fabs(x))) break;
    fx = inv_f(x, &ip);
    if (!isfinite(fx)) return NAN;
    if (fx == 0.0) break;
  }
  return exp(x);
}

/* ------------------------------------------------------------------ */
/* 球対称崩壊                                                          */
/* ------------------------------------------------------------------ */
typedef struct {
  const Bg* b;
  double R3, R32; /* R_ta^3, R_ta^{3/2} */
} TaPar;

static double f_ta(double th, const void* p) {
  const TaPar* tp = p;
  const double s2 = sin(th) * sin(th);
  return 2.0 * tp->R32 * s2 / sqrt(tp->b->om - tp->b->ol * tp->R3 * s2 * (1.0 + s2));
}

static double t_turnaround(const Bg* b, double R) {
  const TaPar tp = {b, R * R * R, R * sqrt(R)};
  return integrate(f_ta, &tp, 0.0, 0.5 * PI);
}

typedef struct {
  const Bg* b;
  double t_c;
} CollPar;

static double coll_f(double R, const void* p) {
  const CollPar* cp = p;
  return 2.0 * t_turnaround(cp->b, R) - cp->t_c;
}

/* a_c で崩壊するシェルのターンアラウンド半径 (Lagrange 半径 = 1 の単位) */
static double R_turnaround(const Bg* b, double a_c) {
  if (!(a_c > 0.0)) return NAN;
  const CollPar par = {b, time_a(b, a_c)};
  if (!isfinite(par.t_c)) return NAN;
  double hi, fhi;
  if (b->ol > 0.0) {
    /* R_ta^3 < Ωm/(2ΩΛ) でないとシェルは引き返さない */
    hi = cbrt(b->om / (2.0 * b->ol)) * (1.0 - 1e-12);
    fhi = coll_f(hi, &par);
    if (!(fhi >= 0.0)) return NAN;
  } else {
    hi = 2.0 * cbrt(2.0 / (3.0 * PI) * 2.0 / (3.0 * PI)) * a_c; /* EdS 解の 2 倍 */
    fhi = coll_f(hi, &par);
    for (int i = 0; fhi < 0.0; ++i) {
      if (i > 200 || !isfinite(fhi)) return NAN;
      hi *= 2.0;
      fhi = coll_f(hi, &par);
    }
  }
  return brent(coll_f, &par, 0.0, hi, -par.t_c, fhi);
}

static double delta_crit_a(const Bg* b, double a_c) {
  const double R = R_turnaround(b, a_c);
  if (!isfinite(R)) return NAN;
  const double dK = b->ok + b->om / R + b->ol * R * R;
  return 0.6 * dK / b->om * D1_a(b, a_c);
}

static double virial_Delta_a(const Bg* b, double a_c) {
  const double R = R_turnaround(b, a_c);
  if (!isfinite(R)) return NAN;
  /* 4η y^3 - 2(1+η) y + 1 = 0 の (0,1) の根 (η < 1/2 なら唯一) */
  const double eta = b->ol * R * R * R / b->om;
  double lo = 0.0, hi = 1.0;
  for (int i = 0; i < 200 && hi - lo > 2.0 * DBL_EPSILON; ++i) {
    const double y = 0.5 * (lo + hi);
    if (4.0 * eta * y * y * y - 2.0 * (1.0 + eta) * y + 1.0 > 0.0) lo = y;
    else hi = y;
  }
  const double y = 0.5 * (lo + hi);
  const double r = a_c / (y * R);
  return r * r * r;
}

/* ------------------------------------------------------------------ */
/* 公開関数                                                            */
/* ------------------------------------------------------------------ */
#define BG_OR_NAN(cp, b) \
  Bg b;                  \
  if (make_bg(cp, &b)) return NAN

double FLRW_omega_m(const cosm_param* cp) {
  BG_OR_NAN(cp, b);
  return b.om;
}

double FLRW_omega_k(const cosm_param* cp) {
  BG_OR_NAN(cp, b);
  return b.ok;
}

double FLRW_hubble_z(const cosm_param* cp, double z) {
  BG_OR_NAN(cp, b);
  return E_a(&b, a_of_z(z));
}

double FLRW_omega_m_z(const cosm_param* cp, double z) {
  BG_OR_NAN(cp, b);
  const double a = a_of_z(z);
  return b.om / a3E2_u(&b, sqrt(a));
}

double FLRW_hubble_time_Gyr(const cosm_param* cp) {
  /* 1/(100 km/s/Mpc) = 977.7922216807891 / 100 Gyr */
  return (cp && cp->hubble > 0.0) ? 9.777922216807891 / cp->hubble : NAN;
}

double FLRW_hubble_distance_Mpc(const cosm_param* cp) {
  return (cp && cp->hubble > 0.0) ? 2997.92458 / cp->hubble : NAN;
}

/* ---- 時刻 ---- */
double FLRW_cosmic_time_z(const cosm_param* cp, double z) {
  BG_OR_NAN(cp, b);
  return time_a(&b, a_of_z(z));
}

double FLRW_conformal_time_z(const cosm_param* cp, double z) {
  BG_OR_NAN(cp, b);
  return conf_a(&b, a_of_z(z));
}

double FLRW_scale_factor_t(const cosm_param* cp, double t) {
  BG_OR_NAN(cp, b);
  return invert(&b, 0, t);
}

double FLRW_redshift_t(const cosm_param* cp, double t) {
  return 1.0 / FLRW_scale_factor_t(cp, t) - 1.0;
}

double FLRW_scale_factor_tau(const cosm_param* cp, double tau) {
  BG_OR_NAN(cp, b);
  return invert(&b, 1, tau);
}

double FLRW_redshift_tau(const cosm_param* cp, double tau) {
  return 1.0 / FLRW_scale_factor_tau(cp, tau) - 1.0;
}

/* ---- 距離: a を与えたときの共通処理 ---- */
enum { DIST_C, DIST_M, DIST_A, DIST_L };

static double dist_from(const Bg* b, double a, double dc, int kind) {
  if (!isfinite(a) || !isfinite(dc)) return NAN;
  if (kind == DIST_C) return dc;
  const double dm = DM_from_DC(b, dc);
  if (kind == DIST_M) return dm;
  return (kind == DIST_A) ? dm * a : dm / a;
}

static double dist_z(const cosm_param* cp, double z, int kind) {
  BG_OR_NAN(cp, b);
  const double a = a_of_z(z);
  return isfinite(a) ? dist_from(&b, a, DC_a(&b, a), kind) : NAN;
}

static double dist_t(const cosm_param* cp, double t, int kind) {
  BG_OR_NAN(cp, b);
  const double a = invert(&b, 0, t);
  return isfinite(a) ? dist_from(&b, a, DC_a(&b, a), kind) : NAN;
}

static double dist_tau(const cosm_param* cp, double tau, int kind) {
  BG_OR_NAN(cp, b);
  const double a = invert(&b, 1, tau);
  return dist_from(&b, a, conf_a(&b, 1.0) - tau, kind); /* DC = tau0 - tau */
}

double FLRW_DC_z(const cosm_param* cp, double z) { return dist_z(cp, z, DIST_C); }
double FLRW_DM_z(const cosm_param* cp, double z) { return dist_z(cp, z, DIST_M); }
double FLRW_DA_z(const cosm_param* cp, double z) { return dist_z(cp, z, DIST_A); }
double FLRW_DL_z(const cosm_param* cp, double z) { return dist_z(cp, z, DIST_L); }

double FLRW_DC_t(const cosm_param* cp, double t) { return dist_t(cp, t, DIST_C); }
double FLRW_DM_t(const cosm_param* cp, double t) { return dist_t(cp, t, DIST_M); }
double FLRW_DA_t(const cosm_param* cp, double t) { return dist_t(cp, t, DIST_A); }
double FLRW_DL_t(const cosm_param* cp, double t) { return dist_t(cp, t, DIST_L); }

double FLRW_DC_tau(const cosm_param* cp, double tau) { return dist_tau(cp, tau, DIST_C); }
double FLRW_DM_tau(const cosm_param* cp, double tau) { return dist_tau(cp, tau, DIST_M); }
double FLRW_DA_tau(const cosm_param* cp, double tau) { return dist_tau(cp, tau, DIST_A); }
double FLRW_DL_tau(const cosm_param* cp, double tau) { return dist_tau(cp, tau, DIST_L); }

/* ---- a の関数を z, t, tau の関数として呼ぶ ---- */
typedef double (*AFn)(const Bg* b, double a);

static double growth_a(const Bg* b, double a) { return D1_a(b, a) / D1_a(b, 1.0); }

static double call_z(const cosm_param* cp, double z, AFn f) {
  BG_OR_NAN(cp, b);
  const double a = a_of_z(z);
  return isfinite(a) ? f(&b, a) : NAN;
}

static double call_t(const cosm_param* cp, double t, AFn f) {
  BG_OR_NAN(cp, b);
  const double a = invert(&b, 0, t);
  return isfinite(a) ? f(&b, a) : NAN;
}

static double call_tau(const cosm_param* cp, double tau, AFn f) {
  BG_OR_NAN(cp, b);
  const double a = invert(&b, 1, tau);
  return isfinite(a) ? f(&b, a) : NAN;
}

double FLRW_growth_factor_z(const cosm_param* cp, double z) { return call_z(cp, z, growth_a); }
double FLRW_growth_factor_t(const cosm_param* cp, double t) { return call_t(cp, t, growth_a); }
double FLRW_growth_factor_tau(const cosm_param* cp, double tau) {
  return call_tau(cp, tau, growth_a);
}

double FLRW_delta_crit_z(const cosm_param* cp, double z) { return call_z(cp, z, delta_crit_a); }
double FLRW_delta_crit_t(const cosm_param* cp, double t) { return call_t(cp, t, delta_crit_a); }
double FLRW_delta_crit_tau(const cosm_param* cp, double tau) {
  return call_tau(cp, tau, delta_crit_a);
}

double FLRW_virial_Delta_z(const cosm_param* cp, double z) {
  return call_z(cp, z, virial_Delta_a);
}
double FLRW_virial_Delta_t(const cosm_param* cp, double t) {
  return call_t(cp, t, virial_Delta_a);
}
double FLRW_virial_Delta_tau(const cosm_param* cp, double tau) {
  return call_tau(cp, tau, virial_Delta_a);
}
