#include "hcephes.h"

#define M 128.0
#define MINV .0078125

HCEPHES_API double hcephes_expx2(double x, int sign) {
    double u, u1, m, f;

    x = fabs(x);
    if (sign < 0)
        x = -x;

    /* Represent x as an exact multiple of M plus a residual.
       M is a power of 2 chosen so that exp(m * m) does not overflow
       or underflow and so that |x - m| is small.  */
    m = MINV * floor(M * x + 0.5);
    f = x - m;

    /* x^2 = m^2 + 2mf + f^2 */
    u = m * m;
    u1 = 2 * m * f + f * f;

    if (sign < 0) {
        u = -u;
        u1 = -u1;
    }

    if ((u + u1) > HCEPHES_MAXLOG)
        return (HUGE_VAL);

    /* u is exact, u1 is small.  */
    u = exp(u) * exp(u1);
    return (u);
}
