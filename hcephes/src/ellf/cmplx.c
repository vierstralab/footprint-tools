#include "hcephes.h"

#ifdef _MSC_VER
#if (_MSC_VER <= 1500)
#define isnan(x) _isnan(x)
#endif
#endif

cmplx czero = {0.0, 0.0};
extern cmplx czero;
cmplx cone = {1.0, 0.0};
extern cmplx cone;

/*	c = b + a	*/

HCEPHES_API void hcephes_cadd(register cmplx *a, register cmplx *b, cmplx *c) {

    c->r = b->r + a->r;
    c->i = b->i + a->i;
}

/*	c = b - a	*/

HCEPHES_API void hcephes_csub(register cmplx *a, register cmplx *b, cmplx *c) {

    c->r = b->r - a->r;
    c->i = b->i - a->i;
}

/*	c = b * a */

HCEPHES_API void hcephes_cmul(register cmplx *a, register cmplx *b, cmplx *c) {
    double y;

    y = b->r * a->r - b->i * a->i;
    c->i = b->r * a->i + b->i * a->r;
    c->r = y;
}

/*	c = b / a */

HCEPHES_API void hcephes_cdiv(register cmplx *a, register cmplx *b, cmplx *c) {
    double y, p, q, w;

    y = a->r * a->r + a->i * a->i;
    p = b->r * a->r + b->i * a->i;
    q = b->i * a->r - b->r * a->i;

    if (y < 1.0) {
        w = HUGE_VAL * y;
        if ((fabs(p) > w) || (fabs(q) > w) || (y == 0.0)) {
            c->r = HUGE_VAL;
            c->i = HUGE_VAL;
            hcephes_mtherr("cdiv", HCEPHES_OVERFLOW);
            return;
        }
    }
    c->r = p / y;
    c->i = q / y;
}

/*	b = a
   Caution, a `short' is assumed to be 16 bits wide.  */

HCEPHES_API void hcephes_cmov(void *a, void *b) {
    register short *pa, *pb;
    int i;

    pa = (short *)a;
    pb = (short *)b;
    i = 8;
    do
        *pb++ = *pa++;
    while (--i);
}

HCEPHES_API void hcephes_cneg(register cmplx *a) {

    a->r = -a->r;
    a->i = -a->i;
}

#define PREC 27
#define MAXEXP 1024
#define MINEXP -1077

HCEPHES_API double hcephes_cabs(register cmplx *z) {
    double x, y, b, re, im;
    int ex, ey, e;

    /* Note, hcephes_cabs(HUGE_VAL,HCEPHES_NAN) = HUGE_VAL. */
    if (z->r == HUGE_VAL || z->i == HUGE_VAL || z->r == -HUGE_VAL || z->i == -HUGE_VAL)
        return (HUGE_VAL);

    if (isnan(z->r))
        return (z->r);
    if (isnan(z->i))
        return (z->i);

    re = fabs(z->r);
    im = fabs(z->i);

    if (re == 0.0)
        return (im);
    if (im == 0.0)
        return (re);

    /* Get the exponents of the numbers */
    x = frexp(re, &ex);
    y = frexp(im, &ey);

    /* Check if one number is tiny compared to the other */
    e = ex - ey;
    if (e > PREC)
        return (re);
    if (e < -PREC)
        return (im);

    /* Find approximate exponent e of the geometric mean. */
    e = (ex + ey) >> 1;

    /* Rescale so mean is about 1 */
    x = ldexp(re, -e);
    y = ldexp(im, -e);

    /* Hypotenuse of the right triangle */
    b = sqrt(x * x + y * y);

    /* Compute the exponent of the answer. */
    y = frexp(b, &ey);
    ey = e + ey;

    /* Check it for overflow and underflow. */
    if (ey > MAXEXP) {
        hcephes_mtherr("cabs", HCEPHES_OVERFLOW);
        return (HUGE_VAL);
    }
    if (ey < MINEXP)
        return (0.0);

    /* Undo the scaling */
    b = ldexp(b, e);
    return (b);
}

HCEPHES_API void hcephes_csqrt(cmplx *z, cmplx *w) {
    cmplx q, s;
    double x, y, r, t;

    x = z->r;
    y = z->i;

    if (y == 0.0) {
        if (x < 0.0) {
            w->r = 0.0;
            w->i = sqrt(-x);
            return;
        } else {
            w->r = sqrt(x);
            w->i = 0.0;
            return;
        }
    }

    if (x == 0.0) {
        r = fabs(y);
        r = sqrt(0.5 * r);
        if (y > 0)
            w->r = r;
        else
            w->r = -r;
        w->i = r;
        return;
    }

    /* Approximate  sqrt(x^2+y^2) - x  =  y^2/2x - y^4/24x^3 + ... .
     * The relative error in the first term is approximately y^2/12x^2 .
     */
    if ((fabs(y) < 2.e-4 * fabs(x)) && (x > 0)) {
        t = 0.25 * y * (y / x);
    } else {
        r = hcephes_cabs(z);
        t = 0.5 * (r - x);
    }

    r = sqrt(t);
    q.i = r;
    q.r = y / (2.0 * r);
    /* Heron iteration in complex arithmetic */
    hcephes_cdiv(&q, z, &s);
    hcephes_cadd(&q, &s, w);
    w->r *= 0.5;
    w->i *= 0.5;
}

HCEPHES_API double hcephes_hypot(double x, double y) {
    cmplx z;

    z.r = x;
    z.i = y;
    return (hcephes_cabs(&z));
}
