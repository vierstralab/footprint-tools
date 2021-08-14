#include "hcephes.h"

HCEPHES_API double hcephes_gdtr(double a, double b, double x) {

    if (x < 0.0) {
        hcephes_mtherr("gdtr", HCEPHES_DOMAIN);
        return (0.0);
    }
    return (hcephes_igam(b, a * x));
}

HCEPHES_API double hcephes_gdtrc(double a, double b, double x) {

    if (x < 0.0) {
        hcephes_mtherr("gdtrc", HCEPHES_DOMAIN);
        return (0.0);
    }
    return (hcephes_igamc(b, a * x));
}
