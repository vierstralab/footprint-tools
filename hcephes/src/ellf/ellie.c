#include "hcephes.h"

HCEPHES_API double hcephes_ellie(double phi, double m) {
    double a, b, c, e, temp;
    double lphi, t, E;
    int d, mod, npio2, sign;

    if (m == 0.0)
        return (phi);
    lphi = phi;
    npio2 = (int)floor(lphi / HCEPHES_PIO2);
    if (npio2 & 1)
        npio2 += 1;
    lphi = lphi - npio2 * HCEPHES_PIO2;
    if (lphi < 0.0) {
        lphi = -lphi;
        sign = -1;
    } else {
        sign = 1;
    }
    a = 1.0 - m;
    E = hcephes_ellpe(m);
    if (a == 0.0) {
        temp = sin(lphi);
        goto done;
    }
    t = tan(lphi);
    b = sqrt(a);
    /* Thanks to Brian Fitzgerald <fitzgb@mml0.meche.rpi.edu>
       for pointing out an instability near odd multiples of pi/2.  */
    if (fabs(t) > 10.0) {
        /* Transform the amplitude */
        e = 1.0 / (b * t);
        /* ... but avoid multiple recursions.  */
        if (fabs(e) < 10.0) {
            e = atan(e);
            temp = E + m * sin(lphi) * sin(e) - hcephes_ellie(e, m);
            goto done;
        }
    }
    c = sqrt(m);
    a = 1.0;
    d = 1;
    e = 0.0;
    mod = 0;

    while (fabs(c / a) > HCEPHES_MACHEP) {
        temp = b / a;
        lphi = lphi + atan(t * temp) + mod * HCEPHES_PI;
        mod = (int)((lphi + HCEPHES_PIO2) / HCEPHES_PI);
        t = t * (1.0 + temp) / (1.0 - temp * t * t);
        c = (a - b) / 2.0;
        temp = sqrt(a * b);
        a = (a + b) / 2.0;
        b = temp;
        d += d;
        e += c * sin(lphi);
    }

    temp = E / hcephes_ellpk(m);
    temp *= (atan(t) + mod * HCEPHES_PI) / (d * a);
    temp += e;

done:

    if (sign < 0)
        temp = -temp;
    temp += npio2 * E;
    return (temp);
}
