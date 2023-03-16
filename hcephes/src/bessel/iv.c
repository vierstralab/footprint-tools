#include "hcephes.h"

HCEPHES_API double hcephes_iv(double v, double x) {
    int sign;
    double t, ax;

    /* If v is a negative integer, invoke symmetry */
    t = floor(v);
    if (v < 0.0) {
        if (t == v) {
            v = -v; /* symmetry */
            t = -t;
        }
    }
    /* If x is negative, require v to be an integer */
    sign = 1;
    if (x < 0.0) {
        if (t != v) {
            hcephes_mtherr("iv", HCEPHES_DOMAIN);
            return (0.0);
        }
        if (v != 2.0 * floor(v / 2.0))
            sign = -1;
    }

    /* Avoid logarithm singularity */
    if (x == 0.0) {
        if (v == 0.0)
            return (1.0);
        if (v < 0.0) {
            hcephes_mtherr("iv", HCEPHES_OVERFLOW);
            return (HUGE_VAL);
        } else
            return (0.0);
    }

    ax = fabs(x);
    t = v * log(0.5 * ax) - x;
    t = sign * exp(t) / hcephes_gamma(v + 1.0);
    ax = v + 0.5;
    return (t * hcephes_hyperg(ax, 2.0 * ax, 2.0 * x));
}
