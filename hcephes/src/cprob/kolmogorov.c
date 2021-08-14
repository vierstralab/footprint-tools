#include "hcephes.h"

/* Exact Smirnov statistic, for one-sided test.  */
HCEPHES_API double hcephes_smirnov(int n, double e) {
    int v, nn;
    double evn, omevn, p, t, c, hcephes_lgamnp1;

    if (n <= 0 || e < 0.0 || e > 1.0)
        return (-1.0);
    nn = (int)floor((double)n * (1.0 - e));
    p = 0.0;
    if (n < 1013) {
        c = 1.0;
        for (v = 0; v <= nn; v++) {
            evn = e + ((double)v) / n;
            p += c * pow(evn, (double)(v - 1)) * pow(1.0 - evn, (double)(n - v));
            /* Next combinatorial term; worst case error = 4e-15.  */
            c *= ((double)(n - v)) / (v + 1);
        }
    } else {
        hcephes_lgamnp1 = hcephes_lgam((double)(n + 1));
        for (v = 0; v <= nn; v++) {
            evn = e + ((double)v) / n;
            omevn = 1.0 - evn;
            if (fabs(omevn) > 0.0) {
                t = hcephes_lgamnp1 - hcephes_lgam((double)(v + 1)) -
                    hcephes_lgam((double)(n - v + 1)) + (v - 1) * log(evn) +
                    (n - v) * log(omevn);
                if (t > -HCEPHES_MAXLOG)
                    p += exp(t);
            }
        }
    }
    return (p * e);
}

/* Kolmogorov's limiting distribution of two-sided test, returns
   probability that sqrt(n) * max deviation > y,
   or that max deviation > y/sqrt(n).
   The approximation is useful for the tail of the distribution
   when n is large.  */
HCEPHES_API double hcephes_kolmogorov(double y) {
    double p, t, r, sign, x;

    x = -2.0 * y * y;
    sign = 1.0;
    p = 0.0;
    r = 1.0;
    do {
        t = exp(x * r * r);
        p += sign * t;
        if (t == 0.0)
            break;
        r += 1.0;
        sign = -sign;
    } while ((t / p) > 1.1e-16);
    return (p + p);
}

/* Functional inverse of Smirnov distribution
   finds e such that smirnov(n,e) = p.  */
HCEPHES_API double hcephes_smirnovi(int n, double p) {
    double e, t, dpde;

    if (p <= 0.0 || p > 1.0) {
        hcephes_mtherr("smirnovi", HCEPHES_DOMAIN);
        return 0.0;
    }
    /* Start with approximation p = exp(-2 n e^2).  */
    e = sqrt(-log(p) / (2.0 * n));
    do {
        /* Use approximate derivative in Newton iteration. */
        t = -2.0 * n * e;
        dpde = 2.0 * t * exp(t * e);
        if (fabs(dpde) > 0.0)
            t = (p - hcephes_smirnov(n, e)) / dpde;
        else {
            hcephes_mtherr("smirnovi", HCEPHES_UNDERFLOW);
            return 0.0;
        }
        e = e + t;
        if (e >= 1.0 || e <= 0.0) {
            hcephes_mtherr("smirnovi", HCEPHES_OVERFLOW);
            return 0.0;
        }
    } while (fabs(t / e) > 1e-10);
    return (e);
}

/* Functional inverse of Kolmogorov statistic for two-sided test.
   Finds y such that kolmogorov(y) = p.
   If e = hcephes_smirnovi (n,p), then hcephes_kolmogi(2 * p) / sqrt(n) should
   be close to e.  */
HCEPHES_API double hcephes_kolmogi(double p) {
    double y, t, dpdy;

    if (p <= 0.0 || p > 1.0) {
        hcephes_mtherr("kolmogi", HCEPHES_DOMAIN);
        return 0.0;
    }
    /* Start with approximation p = 2 exp(-2 y^2).  */
    y = sqrt(-0.5 * log(0.5 * p));
    do {
        /* Use approximate derivative in Newton iteration. */
        t = -2.0 * y;
        dpdy = 4.0 * t * exp(t * y);
        if (fabs(dpdy) > 0.0)
            t = (p - hcephes_kolmogorov(y)) / dpdy;
        else {
            hcephes_mtherr("kolmogi", HCEPHES_UNDERFLOW);
            return 0.0;
        }
        y = y + t;
    } while (fabs(t / y) > 1e-10);
    return (y);
}
