#include "hcephes.h"

HCEPHES_API double hcephes_yn(int n, double x) {
    double an, anm1, anm2, r;
    int k, sign;

    if (n < 0) {
        n = -n;
        if ((n & 1) == 0) /* -1**n */
            sign = 1;
        else
            sign = -1;
    } else
        sign = 1;

    if (n == 0)
        return (sign * hcephes_y0(x));
    if (n == 1)
        return (sign * hcephes_y1(x));

    /* test for overflow */
    if (x <= 0.0) {
        hcephes_mtherr("yn", HCEPHES_SING);
        return (-HUGE_VAL);
    }

    /* forward recurrence on n */

    anm2 = hcephes_y0(x);
    anm1 = hcephes_y1(x);
    k = 1;
    r = 2 * k;
    do {
        an = r * anm1 / x - anm2;
        anm2 = anm1;
        anm1 = an;
        r += 2.0;
        ++k;
    } while (k < n);

    return (sign * an);
}
