#include "hcephes.h"

/*  NIST value (1999): 2 pi h c^2 = 3.741 7749(22) �� 10-16 W m2  */
double planck_c1 = 3.7417749e-16;
/*  NIST value (1999):  h c / k  = 0.014 387 69 m K */
double planck_c2 = 0.01438769;

HCEPHES_API double hcephes_plancki(double w, double T) {
    double b, h, y, bw;

    b = T / planck_c2;
    bw = b * w;

    if (bw > 0.59375) {
        y = b * b;
        h = y * y;
        /* Right tail.  */
        y = hcephes_planckc(w, T);
        /* pi^4 / 15  */
        y = 6.493939402266829149096 * planck_c1 * h - y;
        return y;
    }

    h = exp(-planck_c2 / (w * T));
    y = 6. * hcephes_polylog(4, h) * bw;
    y = (y + 6. * hcephes_polylog(3, h)) * bw;
    y = (y + 3. * hcephes_polylog(2, h)) * bw;
    y = (y - hcephes_log1p(-h)) * bw;
    h = w * w;
    h = h * h;
    y = y * (planck_c1 / h);
    return y;
}

HCEPHES_API double hcephes_planckc(double w, double T) {
    double b, d, p, u, y;

    b = T / planck_c2;
    d = b * w;
    if (d <= 0.59375) {
        y = 6.493939402266829149096 * planck_c1 * b * b * b * b;
        return (y - hcephes_plancki(w, T));
    }
    u = 1.0 / d;
    p = u * u;

    y = -236364091. * p / 45733251691757079075225600000.;
    y = (y + 77683. / 352527500984795136000000.) * p;
    y = (y - 174611. / 18465726242060697600000.) * p;
    y = (y + 43867. / 107290978560589824000.) * p;
    y = ((y - 3617. / 202741834014720000.) * p + 1. / 1270312243200.) * p;
    y = ((y - 691. / 19615115520000.) * p + 1. / 622702080.) * p;
    y = ((((y - 1. / 13305600.) * p + 1. / 272160.) * p - 1. / 5040.) * p + 1. / 60.) *
        p;
    y = y - 0.125 * u + 1. / 3.;
    y = y * planck_c1 * b / (w * w * w);
    return y;
}

HCEPHES_API double hcephes_planckd(double w, double T) {
    return (planck_c2 / ((w * w * w * w * w) * (exp(planck_c2 / (w * T)) - 1.0)));
}

/* Wavelength, w, of maximum radiation at given temperature T.
   c2/wT = constant
   Wein displacement law.
  */
HCEPHES_API double hcephes_planckw(double T) {
    return (planck_c2 / (4.96511423174427630 * T));
}
