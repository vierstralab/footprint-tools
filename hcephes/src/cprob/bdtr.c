#include "hcephes.h"

HCEPHES_API double hcephes_bdtrc(int k, int n, double p) {
    double dk, dn;

    if ((p < 0.0) || (p > 1.0))
        goto domerr;
    if (k < 0)
        return (1.0);

    if (n < k) {
    domerr:
        hcephes_mtherr("bdtrc", HCEPHES_DOMAIN);
        return (0.0);
    }

    if (k == n)
        return (0.0);
    dn = n - k;
    if (k == 0) {
        if (p < .01)
            dk = -hcephes_expm1(dn * hcephes_log1p(-p));
        else
            dk = 1.0 - pow(1.0 - p, dn);
    } else {
        dk = k + 1;
        dk = hcephes_incbet(dk, dn, p);
    }
    return (dk);
}

HCEPHES_API double hcephes_bdtr(int k, int n, double p) {
    double dk, dn;

    if ((p < 0.0) || (p > 1.0))
        goto domerr;
    if ((k < 0) || (n < k)) {
    domerr:
        hcephes_mtherr("bdtr", HCEPHES_DOMAIN);
        return (0.0);
    }

    if (k == n)
        return (1.0);

    dn = n - k;
    if (k == 0) {
        dk = pow(1.0 - p, dn);
    } else {
        dk = k + 1;
        dk = hcephes_incbet(dn, dk, 1.0 - p);
    }
    return (dk);
}

HCEPHES_API double hcephes_bdtri(int k, int n, double y) {
    double dk, dn, p;

    if ((y < 0.0) || (y > 1.0))
        goto domerr;
    if ((k < 0) || (n <= k)) {
    domerr:
        hcephes_mtherr("bdtri", HCEPHES_DOMAIN);
        return (0.0);
    }

    dn = n - k;
    if (k == 0) {
        if (y > 0.8)
            p = -hcephes_expm1(hcephes_log1p(y - 1.0) / dn);
        else
            p = 1.0 - pow(y, 1.0 / dn);
    } else {
        dk = k + 1;
        p = hcephes_incbet(dn, dk, 0.5);
        if (p > 0.5)
            p = hcephes_incbi(dk, dn, 1.0 - y);
        else
            p = 1.0 - hcephes_incbi(dn, dk, y);
    }
    return (p);
}
