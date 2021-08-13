#include "hcephes.h"

HCEPHES_API double hcephes_stdtr(int k, double t) {
    double x, rk, z, f, tz, p, xsqk;
    int j;

    if (k <= 0) {
        hcephes_mtherr("stdtr", HCEPHES_DOMAIN);
        return (0.0);
    }

    if (t == 0)
        return (0.5);

    if (t < -2.0) {
        rk = k;
        z = rk / (rk + t * t);
        p = 0.5 * hcephes_incbet(0.5 * rk, 0.5, z);
        return (p);
    }

    /*	compute integral from -t to + t */

    if (t < 0)
        x = -t;
    else
        x = t;

    rk = k; /* degrees of freedom */
    z = 1.0 + (x * x) / rk;

    /* test if k is odd or even */
    if ((k & 1) != 0) {

        /*	computation for odd k	*/

        xsqk = x / sqrt(rk);
        p = atan(xsqk);
        if (k > 1) {
            f = 1.0;
            tz = 1.0;
            j = 3;
            while ((j <= (k - 2)) && ((tz / f) > HCEPHES_MACHEP)) {
                tz *= (j - 1) / (z * j);
                f += tz;
                j += 2;
            }
            p += f * xsqk / z;
        }
        p *= 2.0 / HCEPHES_PI;
    }

    else {

        /*	computation for even k	*/

        f = 1.0;
        tz = 1.0;
        j = 2;

        while ((j <= (k - 2)) && ((tz / f) > HCEPHES_MACHEP)) {
            tz *= (j - 1) / (z * j);
            f += tz;
            j += 2;
        }
        p = f * x / sqrt(z * rk);
    }

    /*	common exit	*/

    if (t < 0)
        p = -p; /* note destruction of relative accuracy */

    p = 0.5 + 0.5 * p;
    return (p);
}

HCEPHES_API double hcephes_stdtri(int k, double p) {
    double t, rk, z;
    int rflg;

    if (k <= 0 || p <= 0.0 || p >= 1.0) {
        hcephes_mtherr("stdtri", HCEPHES_DOMAIN);
        return (0.0);
    }

    rk = k;

    if (p > 0.25 && p < 0.75) {
        if (p == 0.5)
            return (0.0);
        z = 1.0 - 2.0 * p;
        z = hcephes_incbi(0.5, 0.5 * rk, fabs(z));
        t = sqrt(rk * z / (1.0 - z));
        if (p < 0.5)
            t = -t;
        return (t);
    }
    rflg = -1;
    if (p >= 0.5) {
        p = 1.0 - p;
        rflg = 1;
    }
    z = hcephes_incbi(0.5 * rk, 0.5, 2.0 * p);

    if (HUGE_VAL * z < rk)
        return (rflg * HUGE_VAL);
    t = sqrt(rk / z - rk);
    return (rflg * t);
}
