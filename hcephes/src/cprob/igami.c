#include "hcephes.h"

HCEPHES_API double hcephes_igami(double a, double y0) {
    double x0, x1, x, yl, yh, y, d, lgm, dithresh;
    int i, dir;

    /* bound the solution */
    x0 = HUGE_VAL;
    yl = 0;
    x1 = 0;
    yh = 1.0;
    dithresh = 5.0 * HCEPHES_MACHEP;

    /* approximation to inverse function */
    d = 1.0 / (9.0 * a);
    y = (1.0 - d - hcephes_ndtri(y0) * sqrt(d));
    x = a * y * y * y;

    lgm = hcephes_lgam(a);

    for (i = 0; i < 10; i++) {
        if (x > x0 || x < x1)
            goto ihalve;
        y = hcephes_igamc(a, x);
        if (y < yl || y > yh)
            goto ihalve;
        if (y < y0) {
            x0 = x;
            yl = y;
        } else {
            x1 = x;
            yh = y;
        }
        /* compute the derivative of the function at this point */
        d = (a - 1.0) * log(x) - x - lgm;
        if (d < -HCEPHES_MAXLOG)
            goto ihalve;
        d = -exp(d);
        /* compute the step to the next approximation of x */
        d = (y - y0) / d;
        if (fabs(d / x) < HCEPHES_MACHEP)
            goto done;
        x = x - d;
    }

/* Resort to interval halving if Newton iteration did not converge. */
ihalve:

    d = 0.0625;
    if (x0 == HUGE_VAL) {
        if (x <= 0.0)
            x = 1.0;
        while (x0 == HUGE_VAL) {
            x = (1.0 + d) * x;
            y = hcephes_igamc(a, x);
            if (y < y0) {
                x0 = x;
                yl = y;
                break;
            }
            d = d + d;
        }
    }
    d = 0.5;
    dir = 0;

    for (i = 0; i < 400; i++) {
        x = x1 + d * (x0 - x1);
        y = hcephes_igamc(a, x);
        lgm = (x0 - x1) / (x1 + x0);
        if (fabs(lgm) < dithresh)
            break;
        lgm = (y - y0) / y0;
        if (fabs(lgm) < dithresh)
            break;
        if (x <= 0.0)
            break;
        if (y >= y0) {
            x1 = x;
            yh = y;
            if (dir < 0) {
                dir = 0;
                d = 0.5;
            } else if (dir > 1)
                d = 0.5 * d + 0.5;
            else
                d = (y0 - yl) / (yh - yl);
            dir += 1;
        } else {
            x0 = x;
            yl = y;
            if (dir > 0) {
                dir = 0;
                d = 0.5;
            } else if (dir < -1)
                d = 0.5 * d;
            else
                d = (y0 - yl) / (yh - yl);
            dir -= 1;
        }
    }
    if (x == 0.0)
        hcephes_mtherr("igami", HCEPHES_UNDERFLOW);

done:
    return (x);
}
