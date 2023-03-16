#include "hcephes.h"

HCEPHES_API double hcephes_fdtrc(int ia, int ib, double x) {
    double a, b, w;

    if ((ia < 1) || (ib < 1) || (x < 0.0)) {
        hcephes_mtherr("fdtrc", HCEPHES_DOMAIN);
        return (0.0);
    }
    a = ia;
    b = ib;
    w = b / (b + a * x);
    return (hcephes_incbet(0.5 * b, 0.5 * a, w));
}

HCEPHES_API double hcephes_fdtr(int ia, int ib, double x) {
    double a, b, w;

    if ((ia < 1) || (ib < 1) || (x < 0.0)) {
        hcephes_mtherr("fdtr", HCEPHES_DOMAIN);
        return (0.0);
    }
    a = ia;
    b = ib;
    w = a * x;
    w = w / (b + w);
    return (hcephes_incbet(0.5 * a, 0.5 * b, w));
}

HCEPHES_API double hcephes_fdtri(int ia, int ib, double y) {
    double a, b, w, x;

    // added by danilo
    y = 1.0 - y;

    if ((ia < 1) || (ib < 1) || (y <= 0.0) || (y > 1.0)) {
        hcephes_mtherr("fdtri", HCEPHES_DOMAIN);
        return (0.0);
    }
    a = ia;
    b = ib;
    /* Compute probability for x = 0.5.  */
    w = hcephes_incbet(0.5 * b, 0.5 * a, 0.5);
    /* If that is greater than y, then the solution w < .5.
       Otherwise, solve at 1-y to remove cancellation in (b - b*w).  */
    if (w > y || y < 0.001) {
        w = hcephes_incbi(0.5 * b, 0.5 * a, y);
        x = (b - b * w) / (a * w);
    } else {
        w = hcephes_incbi(0.5 * a, 0.5 * b, 1.0 - y);
        x = b * w / (a * (1.0 - w));
    }
    return (x);
}
