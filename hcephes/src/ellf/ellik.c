#include "hcephes.h"

HCEPHES_API double hcephes_ellik(double phi, double m) {
    double a, b, c, e, temp, t, K;
    int d, mod, sign, npio2;

    if (m == 0.0)
        return (phi);
    a = 1.0 - m;
    if (a == 0.0) {
        if (fabs(phi) >= HCEPHES_PIO2) {
            hcephes_mtherr("ellik", HCEPHES_SING);
            return (HUGE_VAL);
        }
        return (log(tan((HCEPHES_PIO2 + phi) / 2.0)));
    }
    npio2 = (int)floor(phi / HCEPHES_PIO2);
    if (npio2 & 1)
        npio2 += 1;
    if (npio2) {
        K = hcephes_ellpk(m);
        phi = phi - npio2 * HCEPHES_PIO2;
    } else
        K = 0.0;
    if (phi < 0.0) {
        phi = -phi;
        sign = -1;
    } else
        sign = 0;
    b = sqrt(a);
    t = tan(phi);
    if (fabs(t) > 10.0) {
        /* Transform the amplitude */
        e = 1.0 / (b * t);
        /* ... but avoid multiple recursions.  */
        if (fabs(e) < 10.0) {
            e = atan(e);
            if (npio2 == 0)
                K = hcephes_ellpk(m);
            temp = K - hcephes_ellik(e, m);
            goto done;
        }
    }
    a = 1.0;
    c = sqrt(m);
    d = 1;
    mod = 0;

    while (fabs(c / a) > HCEPHES_MACHEP) {
        temp = b / a;
        phi = phi + atan(t * temp) + mod * HCEPHES_PI;
        mod = (int)((phi + HCEPHES_PIO2) / HCEPHES_PI);
        t = t * (1.0 + temp) / (1.0 - temp * t * t);
        c = (a - b) / 2.0;
        temp = sqrt(a * b);
        a = (a + b) / 2.0;
        b = temp;
        d += d;
    }

    temp = (atan(t) + mod * HCEPHES_PI) / (d * a);

done:
    if (sign < 0)
        temp = -temp;
    temp += npio2 * K;
    return (temp);
}
