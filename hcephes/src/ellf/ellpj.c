#include "hcephes.h"

HCEPHES_API int hcephes_ellpj(double u, double m, double *sn, double *cn, double *dn, double *ph) {
    double ai, b, phi, t, twon;
    double a[9], c[9];
    int i;

    /* Check for special cases */

    if (m < 0.0 || m > 1.0) {
        hcephes_mtherr("ellpj", HCEPHES_DOMAIN);
        *sn = 0.0;
        *cn = 0.0;
        *ph = 0.0;
        *dn = 0.0;
        return (-1);
    }
    if (m < 1.0e-9) {
        t = sin(u);
        b = cos(u);
        ai = 0.25 * m * (u - t * b);
        *sn = t - ai * b;
        *cn = b + ai * t;
        *ph = u - ai;
        *dn = 1.0 - 0.5 * m * t * t;
        return (0);
    }

    if (m >= 0.9999999999) {
        ai = 0.25 * (1.0 - m);
        b = cosh(u);
        t = tanh(u);
        phi = 1.0 / b;
        twon = b * sinh(u);
        *sn = t + ai * (twon - u) / (b * b);
        *ph = 2.0 * atan(exp(u)) - HCEPHES_PIO2 + ai * (twon - u) / b;
        ai *= t * phi;
        *cn = phi - ai * (twon - u);
        *dn = phi + ai * (twon + u);
        return (0);
    }

    /*	A. G. M. scale		*/
    a[0] = 1.0;
    b = sqrt(1.0 - m);
    c[0] = sqrt(m);
    twon = 1.0;
    i = 0;

    while (fabs(c[i] / a[i]) > HCEPHES_MACHEP) {
        if (i > 7) {
            hcephes_mtherr("ellpj", HCEPHES_OVERFLOW);
            goto done;
        }
        ai = a[i];
        ++i;
        c[i] = (ai - b) / 2.0;
        t = sqrt(ai * b);
        a[i] = (ai + b) / 2.0;
        b = t;
        twon *= 2.0;
    }

done:

    /* backward recurrence */
    phi = twon * a[i] * u;
    do {
        t = c[i] * sin(phi) / a[i];
        b = phi;
        phi = (asin(t) + phi) / 2.0;
    } while (--i);

    t = sin(phi);
    *sn = t;
    *cn = cos(phi);
    /* Thanks to Hartmut Henkel for reporting a bug here:  */
    *dn = sqrt(1.0 - m * t * t);
    *ph = phi;
    return (0);
}
