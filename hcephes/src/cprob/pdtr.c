#include "hcephes.h"

HCEPHES_API double hcephes_pdtrc(int k, double m) {
    double v;

    if ((k < 0) || (m <= 0.0)) {
        hcephes_mtherr("pdtrc", HCEPHES_DOMAIN);
        return (0.0);
    }
    v = k + 1;
    return (hcephes_igam(v, m));
}

HCEPHES_API double hcephes_pdtr(int k, double m) {
    double v;

    if ((k < 0) || (m <= 0.0)) {
        hcephes_mtherr("pdtr", HCEPHES_DOMAIN);
        return (0.0);
    }
    v = k + 1;
    return (hcephes_igamc(v, m));
}

HCEPHES_API double hcephes_pdtri(int k, double y) {
    double v;

    if ((k < 0) || (y < 0.0) || (y >= 1.0)) {
        hcephes_mtherr("pdtri", HCEPHES_DOMAIN);
        return (0.0);
    }
    v = k + 1;
    v = hcephes_igami(v, y);
    return (v);
}
