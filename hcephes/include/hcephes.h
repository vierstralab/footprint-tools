/*******************************************************************************
 * Project: hcephes
 * Purpose: Netlib Cephes updated
 * Author: Danilo Horta, danilo.horta@gmail.com
 * Language: C
 *******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2018 Danilo Horta
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/
#ifndef HCEPHES_H
#define HCEPHES_H

#define HCEPHES_VERSION "0.4.1"
#define HCEPHES_VERSION_MAJOR 0
#define HCEPHES_VERSION_MINOR 4
#define HCEPHES_VERSION_PATCH 1

#ifdef __cplusplus
extern "C" {
#endif

#ifdef _WIN32
#ifdef HCEPHES_EXPORTS
#define HCEPHES_API __declspec(dllexport)
#else
#define HCEPHES_API __declspec(dllimport)
#endif
#else
#define HCEPHES_API
#endif

#include <float.h>
#include <math.h>

#ifndef NAN
#if _MSC_VER == 1500
#define NAN sqrt(-1)
#else
#include <bits/nan.h>
#endif
#endif

/* Constant definitions for math error conditions
 */
#define HCEPHES_DOMAIN 1    /* argument domain error */
#define HCEPHES_SING 2      /* argument singularity */
#define HCEPHES_OVERFLOW 3  /* overflow range error */
#define HCEPHES_UNDERFLOW 4 /* underflow range error */
#define HCEPHES_TLOSS 5     /* total loss of precision */
#define HCEPHES_PLOSS 6     /* partial loss of precision */
#define HCEPHES_EDOM 33
#define HCEPHES_ERANGE 34

#define HCEPHES_LOGE2 6.93147180559945309417E-1   /* log(2) */
#define HCEPHES_LOGSQ2 3.46573590279972654709E-1  /* log(2)/2 */
#define HCEPHES_MACHEP 1.11022302462515654042E-16 /* 2**-53 */
#define HCEPHES_MAXLOG 7.09782712893383996732E2   /* log(DBL_MAX) */
#define HCEPHES_MAXNUM HUGE_VAL
#define HCEPHES_MINLOG -7.451332191019412076235E2 /* log(2**-1075) */
#define HCEPHES_NEGZERO -0.0
#define HCEPHES_PI 3.14159265358979323846
#define HCEPHES_PIO2 (3.14159265358979323846 / 2)
#define HCEPHES_PIO4 (3.14159265358979323846 / 4)
#define HCEPHES_SQ2OPI 7.9788456080286535587989E-1  /* sqrt( 2/pi ) */
#define HCEPHES_SQRTH 7.07106781186547524401E-1     /* sqrt(2)/2 */
#define HCEPHES_THPIO4 2.35619449019234492885       /* 3*pi/4 */
#define HCEPHES_TWOOPI 6.36619772367581343075535E-1 /* 2/pi */

typedef struct {
    double n; /* numerator */
    double d; /* denominator */
} fract;

typedef struct {
    long double r;
    long double i;
} cmplxl;

typedef struct {
    float r;
    float i;
} cmplxf;

typedef struct {
    double r;
    double i;
} cmplx;

HCEPHES_API double hcephes_bdtr(int k, int n, double p);
HCEPHES_API double hcephes_bdtrc(int k, int n, double p);
HCEPHES_API double hcephes_bdtri(int k, int n, double y);
HCEPHES_API double hcephes_beta(double a, double b);
HCEPHES_API double hcephes_btdtr(double a, double b, double x);
HCEPHES_API double hcephes_cabs(cmplx *z);
HCEPHES_API double hcephes_cbrt(double x);
HCEPHES_API double hcephes_chbevl(double x, double array[], int n);
HCEPHES_API double hcephes_chdtr(double df, double x);
HCEPHES_API double hcephes_chdtrc(double df, double x);
HCEPHES_API double hcephes_chdtri(double df, double y);
HCEPHES_API double hcephes_cosm1(double x);
HCEPHES_API double hcephes_dawsn(double xx);
HCEPHES_API double hcephes_ei(double x);
HCEPHES_API double hcephes_ellie(double phi, double m);
HCEPHES_API double hcephes_ellik(double phi, double m);
HCEPHES_API double hcephes_ellpe(double x);
HCEPHES_API double hcephes_ellpk(double x);
HCEPHES_API double hcephes_erf(double a);
HCEPHES_API double hcephes_erfc(double a);
HCEPHES_API double hcephes_erfce(double x);
HCEPHES_API double hcephes_euclid(double *num, double *den);
HCEPHES_API double hcephes_expm1(double x);
HCEPHES_API double hcephes_expn(int n, double x);
HCEPHES_API double hcephes_expx2(double x, int sign);
HCEPHES_API double hcephes_fac(int i);
HCEPHES_API double hcephes_fdtr(int ia, int ib, double x);
HCEPHES_API double hcephes_fdtrc(int ia, int ib, double x);
HCEPHES_API double hcephes_fdtri(int ia, int ib, double y);
HCEPHES_API double hcephes_gamma(double x);
HCEPHES_API double hcephes_gdtr(double a, double b, double x);
HCEPHES_API double hcephes_gdtrc(double a, double b, double x);
HCEPHES_API double hcephes_hyp2f0(double a, double b, double x, int type, double *err);
HCEPHES_API double hcephes_hyp2f1(double a, double b, double c, double x);
HCEPHES_API double hcephes_hyperg(double a, double b, double x);
HCEPHES_API double hcephes_hypot(double x, double y);
HCEPHES_API double hcephes_i0(double x);
HCEPHES_API double hcephes_i0e(double x);
HCEPHES_API double hcephes_i1(double x);
HCEPHES_API double hcephes_i1e(double x);
HCEPHES_API double hcephes_igam(double a, double x);
HCEPHES_API double hcephes_igamc(double a, double x);
HCEPHES_API double hcephes_igami(double a, double y0);
HCEPHES_API double hcephes_incbet(double aa, double bb, double xx);
HCEPHES_API double hcephes_incbi(double aa, double bb, double yy0);
HCEPHES_API double hcephes_iv(double v, double x);
HCEPHES_API double hcephes_j0(double x);
HCEPHES_API double hcephes_j1(double x);
HCEPHES_API double hcephes_jn(int n, double x);
HCEPHES_API double hcephes_jv(double n, double x);
HCEPHES_API double hcephes_k0(double x);
HCEPHES_API double hcephes_k0e(double x);
HCEPHES_API double hcephes_k1(double x);
HCEPHES_API double hcephes_k1e(double x);
HCEPHES_API double hcephes_kn(int nn, double x);
HCEPHES_API double hcephes_kolmogi(double p);
HCEPHES_API double hcephes_kolmogorov(double y);
HCEPHES_API double hcephes_lbeta(double a, double b);
HCEPHES_API double hcephes_lgam(double x);
HCEPHES_API double hcephes_lgam_sgn(double x, int *sign);
HCEPHES_API double hcephes_log1p(double x);
HCEPHES_API double hcephes_nbdtr(int k, int n, double p);
HCEPHES_API double hcephes_nbdtrc(int k, int n, double p);
HCEPHES_API double hcephes_nbdtri(int k, int n, double p);
HCEPHES_API double hcephes_ndtr(double a);
HCEPHES_API double hcephes_ndtri(double y0);
HCEPHES_API double hcephes_onef2(double a, double b, double c, double x, double *err);
HCEPHES_API double hcephes_p1evl(double x, double coef[], int N);
HCEPHES_API double hcephes_pdtr(int k, double m);
HCEPHES_API double hcephes_pdtrc(int k, double m);
HCEPHES_API double hcephes_pdtri(int k, double y);
HCEPHES_API double hcephes_planckc(double w, double T);
HCEPHES_API double hcephes_planckd(double w, double T);
HCEPHES_API double hcephes_plancki(double w, double T);
HCEPHES_API double hcephes_planckw(double T);
HCEPHES_API double hcephes_poleva(double a[], int na, double x);
HCEPHES_API double hcephes_polevl(double x, double coef[], int N);
HCEPHES_API double hcephes_polylog(int n, double x);
HCEPHES_API double hcephes_powi(double x, int nn);
HCEPHES_API double hcephes_psi(double x);
HCEPHES_API double hcephes_simpsn(double f[], double delta);
HCEPHES_API double hcephes_smirnov(int n, double e);
HCEPHES_API double hcephes_smirnovi(int n, double p);
HCEPHES_API double hcephes_spence(double x);
HCEPHES_API double hcephes_stdtr(int k, double t);
HCEPHES_API double hcephes_stdtri(int k, double p);
HCEPHES_API double hcephes_struve(double v, double x);
HCEPHES_API double hcephes_threef0(double a, double b, double c, double x, double *err);
HCEPHES_API double hcephes_y0(double x);
HCEPHES_API double hcephes_y1(double x);
HCEPHES_API double hcephes_yn(int n, double x);
HCEPHES_API double hcephes_yv(double v, double x);
HCEPHES_API double hcephes_zeta(double x, double q);
HCEPHES_API double hcephes_zetac(double x);
HCEPHES_API double rcephes_gamma(double x);
HCEPHES_API int hcephes_airy(double x, double *ai, double *aip, double *bi, double *bip);
HCEPHES_API int hcephes_drand(double *a);
HCEPHES_API int hcephes_ellpj(double u, double m, double *sn, double *cn, double *dn, double *ph);
HCEPHES_API int hcephes_fresnl(double xxa, double *ssa, double *cca);
HCEPHES_API int hcephes_mtherr(char *name, int code);
HCEPHES_API int hcephes_poldiv(double *, int, double *, int, double *);
HCEPHES_API int hcephes_polrt(double xcof[], double cof[], int m, cmplx root[]);
HCEPHES_API int hcephes_shichi(double x, double *si, double *ci);
HCEPHES_API int hcephes_sici(double x, double *si, double *ci);
HCEPHES_API void hcephes_cadd(cmplx *a, cmplx *b, cmplx *c);
HCEPHES_API void hcephes_cdiv(cmplx *, cmplx *, cmplx *);
HCEPHES_API void hcephes_cmov(void *a, void *b);
HCEPHES_API void hcephes_cmul(cmplx *a, cmplx *b, cmplx *c);
HCEPHES_API void hcephes_cneg(cmplx *a);
HCEPHES_API void hcephes_csqrt(cmplx *z, cmplx *w);
HCEPHES_API void hcephes_csub(cmplx *a, cmplx *b, cmplx *c);
HCEPHES_API void hcephes_poladd(double a[], int na, double b[], int nb, double c[]);
HCEPHES_API void hcephes_polatn(double num[], double den[], double ans[], int nn);
HCEPHES_API void hcephes_polclr(double *, int);
HCEPHES_API void hcephes_polcos(double x[], double y[], int nn);
HCEPHES_API void hcephes_polini(int maxdeg);
HCEPHES_API void hcephes_polmov(double *, int, double *);
HCEPHES_API void hcephes_polmul(double *, int, double *, int, double *);
HCEPHES_API void hcephes_polprt(double a[], int na, int d);
HCEPHES_API void hcephes_polsbt(double a[], int na, double b[], int nb, double c[]);
HCEPHES_API void hcephes_polsin(double x[], double y[], int nn);
HCEPHES_API void hcephes_polsqt(double pol[], double ans[], int nn);
HCEPHES_API void hcephes_polsub(double a[], int na, double b[], int nb, double c[]);
HCEPHES_API void hcephes_radd(fract *f1, fract *f2, fract *f3);
HCEPHES_API void hcephes_rdiv(fract *ff1, fract *ff2, fract *ff3);
HCEPHES_API void hcephes_revers(double y[], double x[], int n);
HCEPHES_API void hcephes_rmul(fract *ff1, fract *ff2, fract *ff3);
HCEPHES_API void hcephes_rsub(fract *f1, fract *f2, fract *f3);

#ifdef __cplusplus
}
#endif

#endif /* end of include guard: HCEPHES_HCEPHES_H */
