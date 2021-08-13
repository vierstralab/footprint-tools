#include "hcephes.h"

HCEPHES_API double hcephes_chdtrc(double df, double x) {

    if ((x < 0.0) || (df < 1.0)) {
        hcephes_mtherr("chdtrc", HCEPHES_DOMAIN);
        return (0.0);
    }
    return (hcephes_igamc(df / 2.0, x / 2.0));
}

HCEPHES_API double hcephes_chdtr(double df, double x) {

    if ((x < 0.0) || (df < 1.0)) {
        hcephes_mtherr("chdtr", HCEPHES_DOMAIN);
        return (0.0);
    }
    return (hcephes_igam(df / 2.0, x / 2.0));
}

HCEPHES_API double hcephes_chdtri(double df, double y) {
    double x;

    if ((y < 0.0) || (y > 1.0) || (df < 1.0)) {
        hcephes_mtherr("chdtri", HCEPHES_DOMAIN);
        return (0.0);
    }

    x = hcephes_igami(0.5 * df, y);
    return (2.0 * x);
}
