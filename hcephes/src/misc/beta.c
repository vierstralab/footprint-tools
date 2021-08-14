#include "hcephes.h"

#define MAXGAM 34.84425627277176174

HCEPHES_API double hcephes_beta(double a, double b) {
    double y;
    int sign;

    sign = 1;

    if (a <= 0.0) {
        if (a == floor(a))
            goto over;
    }
    if (b <= 0.0) {
        if (b == floor(b))
            goto over;
    }

    y = a + b;
    if (fabs(y) > MAXGAM) {
        int sgngam;
        y = hcephes_lgam_sgn(y, &sgngam);
        sign *= sgngam; /* keep track of the sign */
        y = hcephes_lgam_sgn(b, &sgngam) - y;
        sign *= sgngam;
        y = hcephes_lgam_sgn(a, &sgngam) + y;
        sign *= sgngam;
        if (y > HCEPHES_MAXLOG) {
        over:
            hcephes_mtherr("beta", HCEPHES_OVERFLOW);
            return (sign * HUGE_VAL);
        }
        return (sign * exp(y));
    }

    y = hcephes_gamma(y);
    if (y == 0.0)
        goto over;

    if (a > b) {
        y = hcephes_gamma(a) / y;
        y *= hcephes_gamma(b);
    } else {
        y = hcephes_gamma(b) / y;
        y *= hcephes_gamma(a);
    }

    return (y);
}

/* Natural log of |beta|.  Return the sign of beta in sgngam.  */

HCEPHES_API double hcephes_lbeta(double a, double b) {
    double y;
    int sign;

    sign = 1;

    if (a <= 0.0) {
        if (a == floor(a))
            goto over;
    }
    if (b <= 0.0) {
        if (b == floor(b))
            goto over;
    }

    y = a + b;
    if (fabs(y) > MAXGAM) {
        int sgngam;
        y = hcephes_lgam_sgn(y, &sgngam);
        sign *= sgngam; /* keep track of the sign */
        y = hcephes_lgam_sgn(b, &sgngam) - y;
        sign *= sgngam;
        y = hcephes_lgam_sgn(a, &sgngam) + y;
        sign *= sgngam;
        sgngam = sign;
        return (y);
    }

    y = hcephes_gamma(y);
    if (y == 0.0) {
    over:
        hcephes_mtherr("lbeta", HCEPHES_OVERFLOW);
        return (sign * HUGE_VAL);
    }

    if (a > b) {
        y = hcephes_gamma(a) / y;
        y *= hcephes_gamma(b);
    } else {
        y = hcephes_gamma(b) / y;
        y *= hcephes_gamma(a);
    }

    if (y < 0)
        y = -y;

    return (log(y));
}
