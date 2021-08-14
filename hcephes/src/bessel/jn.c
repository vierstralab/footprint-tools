#include "hcephes.h"

HCEPHES_API double hcephes_jn(int n, double x) {
    double pkm2, pkm1, pk, xk, r, ans;
    int k, sign;

    if (n < 0) {
        n = -n;
        if ((n & 1) == 0) /* -1**n */
            sign = 1;
        else
            sign = -1;
    } else
        sign = 1;

    if (x < 0.0) {
        if (n & 1)
            sign = -sign;
        x = -x;
    }

    if (n == 0)
        return (sign * hcephes_j0(x));
    if (n == 1)
        return (sign * hcephes_j1(x));
    if (n == 2)
        return (sign * (2.0 * hcephes_j1(x) / x - hcephes_j0(x)));

    if (x < HCEPHES_MACHEP)
        return (0.0);

    k = 53;

    pk = 2 * (n + k);
    ans = pk;
    xk = x * x;

    do {
        pk -= 2.0;
        ans = pk - (xk / ans);
    } while (--k > 0);
    ans = x / ans;

    /* backward recurrence */

    pk = 1.0;
    pkm1 = 1.0 / ans;
    k = n - 1;
    r = 2 * k;

    do {
        pkm2 = (pkm1 * r - pk * x) / x;
        pk = pkm1;
        pkm1 = pkm2;
        r -= 2.0;
    } while (--k > 0);

    if (fabs(pk) > fabs(pkm1))
        ans = hcephes_j1(x) / pk;
    else
        ans = hcephes_j0(x) / pkm1;
    return (sign * ans);
}
