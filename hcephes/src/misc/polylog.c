#include "hcephes.h"

/* polylog(4, 1-x) = zeta(4) - x zeta(3) + x^2 A4(x)/B4(x)
   0 <= x <= 0.125
   Theoretical peak absolute error 4.5e-18  */
static double A4[13] = {
    3.056144922089490701751E-2,  3.243086484162581557457E-1,
    2.877847281461875922565E-1,  7.091267785886180663385E-2,
    6.466460072456621248630E-3,  2.450233019296542883275E-4,
    4.031655364627704957049E-6,  2.884169163909467997099E-8,
    8.680067002466594858347E-11, 1.025983405866370985438E-13,
    4.233468313538272640380E-17, 4.959422035066206902317E-21,
    1.059365867585275714599E-25,
};
static double B4[12] = {
    /* 1.000000000000000000000E0, */
    2.821262403600310974875E0,   1.780221124881327022033E0,
    3.778888211867875721773E-1,  3.193887040074337940323E-2,
    1.161252418498096498304E-3,  1.867362374829870620091E-5,
    1.319022779715294371091E-7,  3.942755256555603046095E-10,
    4.644326968986396928092E-13, 1.913336021014307074861E-16,
    2.240041814626069927477E-20, 4.784036597230791011855E-25,
};

HCEPHES_API double hcephes_polylog(int n, double x) {
    double h, k, p, s, t, u, xc, z;
    int i, j;

    /*  This recurrence provides formulas for n < 2.

        d                 1
        --   Li (x)  =   ---  Li   (x)  .
        dx     n          x     n-1

    */

    if (n == -1) {
        p = 1.0 - x;
        u = x / p;
        s = u * u + u;
        return s;
    }

    if (n == 0) {
        s = x / (1.0 - x);
        return s;
    }

    /* Not implemented for n < -1.
       Not defined for x > 1.  Use cpolylog if you need that.  */
    if (x > 1.0 || n < -1) {
        hcephes_mtherr("polylog", HCEPHES_DOMAIN);
        return 0.0;
    }

    if (n == 1) {
        s = -log(1.0 - x);
        return s;
    }

    /* Argument +1 */
    if (x == 1.0 && n > 1) {
        s = hcephes_zetac((double)n) + 1.0;
        return s;
    }

    /* Argument -1.
                          1-n
       Li (-z)  = - (1 - 2   ) Li (z)
         n                       n
     */
    if (x == -1.0 && n > 1) {
        /* Li_n(1) = zeta(n) */
        s = hcephes_zetac((double)n) + 1.0;
        s = s * (hcephes_powi(2.0, 1 - n) - 1.0);
        return s;
    }

    /*  Inversion formula:
     *                                                   [n/2]   n-2r
     *                n                  1     n           -  log    (z)
     *  Li (-z) + (-1)  Li (-1/z)  =  - --- log (z)  +  2  >  ----------- Li
     * (-1) n               n              n!                -   (n - 2r)!    2r
     *                                                    r=1
     */
    if (x < -1.0 && n > 1) {
        double q, w;
        int r;

        w = log(-x);
        s = 0.0;
        for (r = 1; r <= n / 2; r++) {
            j = 2 * r;
            p = hcephes_polylog(j, -1.0);
            j = n - j;
            if (j == 0) {
                s = s + p;
                break;
            }
            q = (double)j;
            q = pow(w, q) * p / hcephes_fac(j);
            s = s + q;
        }
        s = 2.0 * s;
        q = hcephes_polylog(n, 1.0 / x);
        if (n & 1)
            q = -q;
        s = s - q;
        s = s - pow(w, (double)n) / hcephes_fac(n);
        return s;
    }

    if (n == 2) {
        if (x < 0.0 || x > 1.0)
            return (hcephes_spence(1.0 - x));
    }

    /*  The power series converges slowly when x is near 1.  For n = 3, this
        identity helps:

        Li (-x/(1-x)) + Li (1-x) + Li (x)
          3               3          3
                       2                               2                 3
         = Li (1) + (pi /6) log(1-x) - (1/2) log(x) log (1-x) + (1/6) log (1-x)
             3
    */

    if (n == 3) {
        p = x * x * x;
        if (x > 0.8) {
            /* Thanks to Oscar van Vlijmen for detecting an error here.  */
            u = log(x);
            s = u * u * u / 6.0;
            xc = 1.0 - x;
            s = s - 0.5 * u * u * log(xc);
            s = s + HCEPHES_PI * HCEPHES_PI * u / 6.0;
            s = s - hcephes_polylog(3, -xc / x);
            s = s - hcephes_polylog(3, xc);
            s = s + hcephes_zetac(3.0);
            s = s + 1.0;
            return s;
        }
        /* Power series  */
        t = p / 27.0;
        t = t + .125 * x * x;
        t = t + x;

        s = 0.0;
        k = 4.0;
        do {
            p = p * x;
            h = p / (k * k * k);
            s = s + h;
            k += 1.0;
        } while (fabs(h / s) > 1.1e-16);
        return (s + t);
    }

    if (n == 4) {
        if (x >= 0.875) {
            u = 1.0 - x;
            s = hcephes_polevl(u, A4, 12) / hcephes_p1evl(u, B4, 12);
            s = s * u * u - 1.202056903159594285400 * u;
            s += 1.0823232337111381915160;
            return s;
        }
        goto pseries;
    }

    if (x < 0.75)
        goto pseries;

    /*  This expansion in powers of log(x) is especially useful when
        x is near 1.

        See also the pari gp calculator.

                          inf                  j
                           -    z(n-j) (log(x))
        polylog(n,x)  =    >   -----------------
                           -           j!
                          j=0

          where

          z(j) = Riemann zeta function (j), j != 1

                                  n-1
                                   -
          z(1) =  -log(-log(x)) +  >  1/k
                                   -
                                  k=1
      */

    z = log(x);
    h = -log(-z);
    for (i = 1; i < n; i++)
        h = h + 1.0 / i;
    p = 1.0;
    s = hcephes_zetac((double)n) + 1.0;
    for (j = 1; j <= n + 1; j++) {
        p = p * z / j;
        if (j == n - 1)
            s = s + h * p;
        else
            s = s + (hcephes_zetac((double)(n - j)) + 1.0) * p;
    }
    j = n + 3;
    z = z * z;
    for (;;) {
        p = p * z / ((j - 1) * j);
        h = (hcephes_zetac((double)(n - j)) + 1.0);
        h = h * p;
        s = s + h;
        if (fabs(h / s) < HCEPHES_MACHEP)
            break;
        j += 2;
    }
    return s;

pseries:

    p = x * x * x;
    k = 3.0;
    s = 0.0;
    do {
        p = p * x;
        k += 1.0;
        h = p / hcephes_powi(k, n);
        s = s + h;
    } while (fabs(h / s) > HCEPHES_MACHEP);
    s += x * x * x / hcephes_powi(3.0, n);
    s += x * x / hcephes_powi(2.0, n);
    s += x;
    return s;
}
