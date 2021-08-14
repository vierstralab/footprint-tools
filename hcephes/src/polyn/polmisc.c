#include "hcephes.h"

#include <stdio.h>
#include <stdlib.h>

/* Highest degree of polynomial to be handled
   by the polyn.c subroutine package.  */
#define N 16
/* Highest degree actually initialized at runtime.  */
extern int MAXPOL;

/* Taylor series coefficients for various functions
 */
double patan[N + 1] = {0.0, 1.0,        0.0, -1.0 / 3.0,  0.0, 1.0 / 5.0,
                       0.0, -1.0 / 7.0, 0.0, 1.0 / 9.0,   0.0, -1.0 / 11.0,
                       0.0, 1.0 / 13.0, 0.0, -1.0 / 15.0, 0.0};

double psin[N + 1] = {0.0, 1.0,
                      0.0, -1.0 / 6.0,
                      0.0, 1.0 / 120.0,
                      0.0, -1.0 / 5040.0,
                      0.0, 1.0 / 362880.0,
                      0.0, -1.0 / 39916800.0,
                      0.0, 1.0 / 6227020800.0,
                      0.0, -1.0 / 1.307674368e12,
                      0.0};

double pcos[N + 1] = {1.0,
                      0.0,
                      -1.0 / 2.0,
                      0.0,
                      1.0 / 24.0,
                      0.0,
                      -1.0 / 720.0,
                      0.0,
                      1.0 / 40320.0,
                      0.0,
                      -1.0 / 3628800.0,
                      0.0,
                      1.0 / 479001600.0,
                      0.0,
                      -1.0 / 8.7179291e10,
                      0.0,
                      1.0 / 2.0922789888e13};

double pasin[N + 1] = {0.0, 1.0,
                       0.0, 1.0 / 6.0,
                       0.0, 3.0 / 40.0,
                       0.0, 15.0 / 336.0,
                       0.0, 105.0 / 3456.0,
                       0.0, 945.0 / 42240.0,
                       0.0, 10395.0 / 599040.0,
                       0.0, 135135.0 / 9676800.0,
                       0.0};

/* Square root of 1 + x.  */
double psqrt[N + 1] = {1.0,
                       1. / 2.,
                       -1. / 8.,
                       1. / 16.,
                       -5. / 128.,
                       7. / 256.,
                       -21. / 1024.,
                       33. / 2048.,
                       -429. / 32768.,
                       715. / 65536.,
                       -2431. / 262144.,
                       4199. / 524288.,
                       -29393. / 4194304.,
                       52003. / 8388608.,
                       -185725. / 33554432.,
                       334305. / 67108864.,
                       -9694845. / 2147483648.};

/* Arctangent of the ratio num/den of two polynomials.
 */
HCEPHES_API void hcephes_polatn(double num[], double den[], double ans[], int nn) {
    double a, t;
    double *polq, *polu, *polt;
    int i;

    if (nn > N) {
        hcephes_mtherr("polatn", HCEPHES_OVERFLOW);
        return;
    }
    /* arctan( a + b ) = arctan(a) + arctan( b/(1 + ab + a**2) ) */
    t = num[0];
    a = den[0];
    if ((t == 0.0) && (a == 0.0)) {
        t = num[1];
        a = den[1];
    }
    t = atan2(t, a); /* arctan(num/den), the ANSI argument order */
    polq = (double *)malloc(((unsigned long)MAXPOL + 1) * sizeof(double));
    polu = (double *)malloc(((unsigned long)MAXPOL + 1) * sizeof(double));
    polt = (double *)malloc(((unsigned long)MAXPOL + 1) * sizeof(double));
    hcephes_polclr(polq, MAXPOL);
    i = hcephes_poldiv(den, nn, num, nn, polq);
    a = polq[0];                    /* a */
    polq[0] = 0.0;                  /* b */
    hcephes_polmov(polq, nn, polu); /* b */
    /* Form the polynomial
       1 + ab + a**2
       where a is a scalar.  */
    for (i = 0; i <= nn; i++)
        polu[i] *= a;
    polu[0] += 1.0 + a * a;
    hcephes_poldiv(polu, nn, polq, nn, polt);  /* divide into b */
    hcephes_polsbt(polt, nn, patan, nn, polu); /* arctan(b)  */
    polu[0] += t;                              /* plus arctan(a) */
    hcephes_polmov(polu, nn, ans);
    free(polt);
    free(polu);
    free(polq);
}

/* Square root of a polynomial.
 * Assumes the lowest degree nonzero term is dominant
 * and of even degree.  An error message is given
 * if the Newton iteration does not converge.
 */
HCEPHES_API void hcephes_polsqt(double pol[], double ans[], int nn) {
    double t;
    double *x, *y;
    int i, n;

    if (nn > N) {
        hcephes_mtherr("polatn", HCEPHES_OVERFLOW);
        return;
    }
    x = (double *)malloc(((unsigned long)MAXPOL + 1) * sizeof(double));
    y = (double *)malloc(((unsigned long)MAXPOL + 1) * sizeof(double));
    hcephes_polmov(pol, nn, x);
    hcephes_polclr(y, MAXPOL);

    /* Find lowest degree nonzero term.  */
    t = 0.0;
    for (n = 0; n < nn; n++) {
        if (x[n] != 0.0)
            goto nzero;
    }
    hcephes_polmov(y, nn, ans);
    return;

nzero:

    if (n > 0) {
        if (n & 1) {
            printf("error, sqrt of odd polynomial\n");
            return;
        }
        /* Divide by x^n.  */
        y[n] = x[n];
        hcephes_poldiv(y, nn, pol, N, x);
    }

    t = x[0];
    for (i = 1; i <= nn; i++)
        x[i] /= t;
    x[0] = 0.0;
    /* series development sqrt(1+x) = 1  +  x / 2  -  x**2 / 8  +  x**3 / 16
       hopes that first (constant) term is greater than what follows   */
    hcephes_polsbt(x, nn, psqrt, nn, y);
    t = sqrt(t);
    for (i = 0; i <= nn; i++)
        y[i] *= t;

    /* If first nonzero coefficient was at degree n > 0, multiply by
       x^(n/2).  */
    if (n > 0) {
        hcephes_polclr(x, MAXPOL);
        x[n / 2] = 1.0;
        hcephes_polmul(x, nn, y, nn, y);
    }

    hcephes_polmov(y, nn, ans);
    free(y);
    free(x);
}

/* Sine of a polynomial.
 * The computation uses
 *     sin(a+b) = sin(a) cos(b) + cos(a) sin(b)
 * where a is the constant term of the polynomial and
 * b is the sum of the rest of the terms.
 * Since sin(b) and cos(b) are computed by series expansions,
 * the value of b should be small.
 */
HCEPHES_API void hcephes_polsin(double x[], double y[], int nn) {
    double a, sc;
    double *w, *c;
    int i;

    if (nn > N) {
        hcephes_mtherr("polatn", HCEPHES_OVERFLOW);
        return;
    }
    w = (double *)malloc(((unsigned long)MAXPOL + 1) * sizeof(double));
    c = (double *)malloc(((unsigned long)MAXPOL + 1) * sizeof(double));
    hcephes_polmov(x, nn, w);
    hcephes_polclr(c, MAXPOL);
    hcephes_polclr(y, nn);
    /* a, in the description, is x[0].  b is the polynomial x - x[0].  */
    a = w[0];
    /* c = cos (b) */
    w[0] = 0.0;
    hcephes_polsbt(w, nn, pcos, nn, c);
    sc = sin(a);
    /* sin(a) cos (b) */
    for (i = 0; i <= nn; i++)
        c[i] *= sc;
    /* y = sin (b)  */
    hcephes_polsbt(w, nn, psin, nn, y);
    sc = cos(a);
    /* cos(a) sin(b) */
    for (i = 0; i <= nn; i++)
        y[i] *= sc;
    hcephes_poladd(c, nn, y, nn, y);
    free(c);
    free(w);
}

/* Cosine of a polynomial.
 * The computation uses
 *     cos(a+b) = cos(a) cos(b) - sin(a) sin(b)
 * where a is the constant term of the polynomial and
 * b is the sum of the rest of the terms.
 * Since sin(b) and cos(b) are computed by series expansions,
 * the value of b should be small.
 */
HCEPHES_API void hcephes_polcos(double x[], double y[], int nn) {
    double a, sc;
    double *w, *c;
    int i;

    if (nn > N) {
        hcephes_mtherr("polatn", HCEPHES_OVERFLOW);
        return;
    }
    w = (double *)malloc(((unsigned long)MAXPOL + 1) * sizeof(double));
    c = (double *)malloc(((unsigned long)MAXPOL + 1) * sizeof(double));
    hcephes_polmov(x, nn, w);
    hcephes_polclr(c, MAXPOL);
    hcephes_polclr(y, nn);
    a = w[0];
    w[0] = 0.0;
    /* c = cos(b)  */
    hcephes_polsbt(w, nn, pcos, nn, c);
    sc = cos(a);
    /* cos(a) cos(b)  */
    for (i = 0; i <= nn; i++)
        c[i] *= sc;
    /* y = sin(b) */
    hcephes_polsbt(w, nn, psin, nn, y);
    sc = sin(a);
    /* sin(a) sin(b) */
    for (i = 0; i <= nn; i++)
        y[i] *= sc;
    hcephes_polsub(y, nn, c, nn, y);
    free(c);
    free(w);
}
