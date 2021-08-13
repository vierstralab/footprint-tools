#include "hcephes.h"

#define BIG (1.0 / HCEPHES_MACHEP)

/* Add fractions. */

HCEPHES_API void hcephes_radd(fract *f1, fract *f2, fract *f3) {
    double gcd, d1, d2, gcn, n1, n2;

    n1 = f1->n;
    d1 = f1->d;
    n2 = f2->n;
    d2 = f2->d;
    if (n1 == 0.0) {
        f3->n = n2;
        f3->d = d2;
        return;
    }
    if (n2 == 0.0) {
        f3->n = n1;
        f3->d = d1;
        return;
    }

    gcd = hcephes_euclid(&d1, &d2); /* common divisors of denominators */
    gcn = hcephes_euclid(&n1, &n2); /* common divisors of numerators */
    /* Note, factoring the numerators
     * makes overflow slightly less likely.
     */
    f3->n = (n1 * d2 + n2 * d1) * gcn;
    f3->d = d1 * d2 * gcd;
    hcephes_euclid(&f3->n, &f3->d);
}

/* Subtract fractions. */

HCEPHES_API void hcephes_rsub(fract *f1, fract *f2, fract *f3) {
    double gcd, d1, d2, gcn, n1, n2;

    n1 = f1->n;
    d1 = f1->d;
    n2 = f2->n;
    d2 = f2->d;
    if (n1 == 0.0) {
        f3->n = n2;
        f3->d = d2;
        return;
    }
    if (n2 == 0.0) {
        f3->n = -n1;
        f3->d = d1;
        return;
    }

    gcd = hcephes_euclid(&d1, &d2);
    gcn = hcephes_euclid(&n1, &n2);
    f3->n = (n2 * d1 - n1 * d2) * gcn;
    f3->d = d1 * d2 * gcd;
    hcephes_euclid(&f3->n, &f3->d);
}

/* Multiply fractions. */

HCEPHES_API void hcephes_rmul(fract *ff1, fract *ff2, fract *ff3) {
    double d1, d2, n1, n2;

    n1 = ff1->n;
    d1 = ff1->d;
    n2 = ff2->n;
    d2 = ff2->d;

    if ((n1 == 0.0) || (n2 == 0.0)) {
        ff3->n = 0.0;
        ff3->d = 1.0;
        return;
    }
    hcephes_euclid(&n1, &d2); /* cross cancel common divisors */
    hcephes_euclid(&n2, &d1);
    ff3->n = n1 * n2;
    ff3->d = d1 * d2;
    /* Report overflow. */
    if ((fabs(ff3->n) >= BIG) || (fabs(ff3->d) >= BIG)) {
        hcephes_mtherr("rmul", HCEPHES_OVERFLOW);
        return;
    }
    /* hcephes_euclid( &ff3->n, &ff3->d );*/
}

/* Divide fractions. */

HCEPHES_API void hcephes_rdiv(fract *ff1, fract *ff2, fract *ff3) {
    double d1, d2, n1, n2;

    n1 = ff1->d; /* Invert ff1, then multiply */
    d1 = ff1->n;
    if (d1 < 0.0) { /* keep denominator positive */
        n1 = -n1;
        d1 = -d1;
    }
    n2 = ff2->n;
    d2 = ff2->d;
    if ((n1 == 0.0) || (n2 == 0.0)) {
        ff3->n = 0.0;
        ff3->d = 1.0;
        return;
    }

    hcephes_euclid(&n1, &d2); /* cross cancel any common divisors */
    hcephes_euclid(&n2, &d1);
    ff3->n = n1 * n2;
    ff3->d = d1 * d2;
    /* Report overflow. */
    if ((fabs(ff3->n) >= BIG) || (fabs(ff3->d) >= BIG)) {
        hcephes_mtherr("rdiv", HCEPHES_OVERFLOW);
        return;
    }
    /* hcephes_euclid( &ff3->n, &ff3->d );*/
}

/* Euclidean algorithm
 *   reduces fraction to lowest terms,
 *   returns greatest common divisor.
 */

HCEPHES_API double hcephes_euclid(double *num, double *den) {
    double n, d, q, r;

    n = *num; /* Numerator. */
    d = *den; /* Denominator. */

    /* Make numbers positive, locally. */
    if (n < 0.0)
        n = -n;
    if (d < 0.0)
        d = -d;

    /* Abort if numbers are too big for integer arithmetic. */
    if ((n >= BIG) || (d >= BIG)) {
        hcephes_mtherr("euclid", HCEPHES_OVERFLOW);
        return (1.0);
    }

    /* Divide by zero, gcd = 1. */
    if (d == 0.0)
        return (1.0);

    /* Zero. Return 0/1, gcd = denominator. */
    if (n == 0.0) {
        /*
                if( *den < 0.0 )
                        *den = -1.0;
                else
                        *den = 1.0;
        */
        *den = 1.0;
        return (d);
    }

    while (d > 0.5) {
        /* Find integer part of n divided by d. */
        q = floor(n / d);
        /* Find remainder after dividing n by d. */
        r = n - d * q;
        /* The next fraction is d/r. */
        n = d;
        d = r;
    }

    if (n < 0.0)
        hcephes_mtherr("euclid", HCEPHES_UNDERFLOW);

    *num /= n;
    *den /= n;
    return (n);
}
