#include "hcephes.h"

HCEPHES_API double hcephes_nbdtrc(int k, int n, double p) {
    double dk, dn;

    if ((p < 0.0) || (p > 1.0))
        goto domerr;
    if (k < 0) {
    domerr:
        hcephes_mtherr("nbdtr", HCEPHES_DOMAIN);
        return (0.0);
    }

    dk = k + 1;
    dn = n;
    return (hcephes_incbet(dk, dn, 1.0 - p));
}

HCEPHES_API double hcephes_nbdtr(int k, int n, double p) {
    double dk, dn;

    if ((p < 0.0) || (p > 1.0))
        goto domerr;
    if (k < 0) {
    domerr:
        hcephes_mtherr("nbdtr", HCEPHES_DOMAIN);
        return (0.0);
    }
    dk = k + 1;
    dn = n;
    return (hcephes_incbet(dn, dk, p));
}

HCEPHES_API double hcephes_nbdtri(int k, int n, double p) {
    double dk, dn, w;

    if ((p < 0.0) || (p > 1.0))
        goto domerr;
    if (k < 0) {
    domerr:
        hcephes_mtherr("nbdtri", HCEPHES_DOMAIN);
        return (0.0);
    }
    dk = k + 1;
    dn = n;
    w = hcephes_incbi(dn, dk, p);
    return (w);
}
