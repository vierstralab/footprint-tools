#include "hcephes.h"

/* Chebyshev coefficients for reciprocal cephes_gamma function
 * in interval 0 to 1.  Function is 1/(x cephes_gamma(x)) - 1
 */

static double R[] = {
    3.13173458231230000000E-17, -6.70718606477908000000E-16, 2.20039078172259550000E-15,
    2.47691630348254132600E-13, -6.60074100411295197440E-12, 5.13850186324226978840E-11,
    1.08965386454418662084E-9,  -3.33964630686836942556E-8,  2.68975996440595483619E-7,
    2.96001177518801696639E-6,  -8.04814124978471142852E-5,  4.16609138709688864714E-4,
    5.06579864028608725080E-3,  -6.41925436109158228810E-2,  -4.98558728684003594785E-3,
    1.27546015610523951063E-1};

static char name[] = "rcephes_gamma";

HCEPHES_API double rcephes_gamma(double x) {
    double w, y, z;
    int sign;

    if (x > 34.84425627277176174) {
        hcephes_mtherr(name, HCEPHES_UNDERFLOW);
        return (1.0 / HUGE_VAL);
    }
    if (x < -34.034) {
        w = -x;
        z = sin(HCEPHES_PI * w);
        if (z == 0.0)
            return (0.0);
        if (z < 0.0) {
            sign = 1;
            z = -z;
        } else
            sign = -1;

        y = log(w * z) - log(HCEPHES_PI) + hcephes_lgam(w);
        if (y < -HCEPHES_MAXLOG) {
            hcephes_mtherr(name, HCEPHES_UNDERFLOW);
            return (sign * 1.0 / HUGE_VAL);
        }
        if (y > HCEPHES_MAXLOG) {
            hcephes_mtherr(name, HCEPHES_OVERFLOW);
            return (sign * HUGE_VAL);
        }
        return (sign * exp(y));
    }
    z = 1.0;
    w = x;

    while (w > 1.0) /* Downward recurrence */
    {
        w -= 1.0;
        z *= w;
    }
    while (w < 0.0) /* Upward recurrence */
    {
        z /= w;
        w += 1.0;
    }
    if (w == 0.0) /* Nonpositive integer */
        return (0.0);
    if (w == 1.0) /* Other integer */
        return (1.0 / z);

    y = w * (1.0 + hcephes_chbevl(4.0 * w - 2.0, R, 16)) / z;
    return (y);
}
