#include "hcephes.h"

static double A[] = {8.33333333333333333333E-2, -2.10927960927960927961E-2,
                     7.57575757575757575758E-3, -4.16666666666666666667E-3,
                     3.96825396825396825397E-3, -8.33333333333333333333E-3,
                     8.33333333333333333333E-2};

#define EUL 0.57721566490153286061

HCEPHES_API double hcephes_psi(double x) {
    double p, q, nz, s, w, y, z;
    int i, n, negative;

    negative = 0;
    nz = 0.0;

    if (x <= 0.0) {
        negative = 1;
        q = x;
        p = floor(q);
        if (p == q) {
            hcephes_mtherr("psi", HCEPHES_SING);
            return (HUGE_VAL);
        }
        /* Remove the zeros of tan(HCEPHES_PI x)
         * by subtracting the nearest integer from x
         */
        nz = q - p;
        if (nz != 0.5) {
            if (nz > 0.5) {
                p += 1.0;
                nz = q - p;
            }
            nz = HCEPHES_PI / tan(HCEPHES_PI * nz);
        } else {
            nz = 0.0;
        }
        x = 1.0 - x;
    }

    /* check for positive integer up to 10 */
    if ((x <= 10.0) && (x == floor(x))) {
        y = 0.0;
        n = (int)x;
        for (i = 1; i < n; i++) {
            w = i;
            y += 1.0 / w;
        }
        y -= EUL;
        goto done;
    }

    s = x;
    w = 0.0;
    while (s < 10.0) {
        w += 1.0 / s;
        s += 1.0;
    }

    if (s < 1.0e17) {
        z = 1.0 / (s * s);
        y = z * hcephes_polevl(z, A, 6);
    } else
        y = 0.0;

    y = log(s) - (0.5 / s) - y - w;

done:

    if (negative) {
        y -= nz;
    }

    return (y);
}
