#include "hcephes.h"

static double stop = 1.37e-17;

HCEPHES_API double hcephes_onef2(double a, double b, double c, double x, double *err) {
    double n, a0, sum, t;
    double an, bn, cn, max, z;

    an = a;
    bn = b;
    cn = c;
    a0 = 1.0;
    sum = 1.0;
    n = 1.0;
    t = 1.0;
    max = 0.0;

    do {
        if (an == 0)
            goto done;
        if (bn == 0)
            goto error;
        if (cn == 0)
            goto error;
        if ((a0 > 1.0e34) || (n > 200))
            goto error;
        a0 *= (an * x) / (bn * cn * n);
        sum += a0;
        an += 1.0;
        bn += 1.0;
        cn += 1.0;
        n += 1.0;
        z = fabs(a0);
        if (z > max)
            max = z;
        if (sum != 0)
            t = fabs(a0 / sum);
        else
            t = z;
    } while (t > stop);

done:

    *err = fabs(HCEPHES_MACHEP * max / sum);

    goto xit;

error:
    *err = 1.0e38;

xit:

    return (sum);
}

HCEPHES_API double hcephes_threef0(double a, double b, double c, double x, double *err) {
    double n, a0, sum, t, conv, conv1;
    double an, bn, cn, max, z;

    an = a;
    bn = b;
    cn = c;
    a0 = 1.0;
    sum = 1.0;
    n = 1.0;
    t = 1.0;
    max = 0.0;
    conv = 1.0e38;
    conv1 = conv;

    do {
        if (an == 0.0)
            goto done;
        if (bn == 0.0)
            goto done;
        if (cn == 0.0)
            goto done;
        if ((a0 > 1.0e34) || (n > 200))
            goto error;
        a0 *= (an * bn * cn * x) / n;
        an += 1.0;
        bn += 1.0;
        cn += 1.0;
        n += 1.0;
        z = fabs(a0);
        if (z > max)
            max = z;
        if (z >= conv) {
            if ((z < max) && (z > conv1))
                goto done;
        }
        conv1 = conv;
        conv = z;
        sum += a0;
        if (sum != 0)
            t = fabs(a0 / sum);
        else
            t = z;
    } while (t > stop);

done:

    t = fabs(HCEPHES_MACHEP * max / sum);

    max = fabs(conv / sum);
    if (max > t)
        t = max;

    goto xit;

error:
    t = 1.0e38;

xit:

    *err = t;
    return (sum);
}

HCEPHES_API double hcephes_struve(double v, double x) {
    double y, ya, f, g, h, t;
    double onef2err, threef0err;

    f = floor(v);
    if ((v < 0) && (v - f == 0.5)) {
        y = hcephes_jv(-v, x);
        f = 1.0 - f;
        g = 2.0 * floor(f / 2.0);
        if (g != f)
            y = -y;
        return (y);
    }
    t = 0.25 * x * x;
    f = fabs(x);
    g = 1.5 * fabs(v);
    if ((f > 30.0) && (f > g)) {
        onef2err = 1.0e38;
        y = 0.0;
    } else {
        y = hcephes_onef2(1.0, 1.5, 1.5 + v, -t, &onef2err);
    }

    if ((f < 18.0) || (x < 0.0)) {
        threef0err = 1.0e38;
        ya = 0.0;
    } else {
        ya = hcephes_threef0(1.0, 0.5, 0.5 - v, -1.0 / t, &threef0err);
    }

    f = sqrt(HCEPHES_PI);
    h = pow(0.5 * x, v - 1.0);

    if (onef2err <= threef0err) {
        g = hcephes_gamma(v + 1.5);
        y = y * h * t / (0.5 * f * g);
        return (y);
    } else {
        g = hcephes_gamma(v + 0.5);
        ya = ya * h / (f * g);
        ya = ya + hcephes_yv(v, x);
        return (ya);
    }
}

/* Bessel function of noninteger order
 */

HCEPHES_API double hcephes_yv(double v, double x) {
    double y, t;
    int n;

    y = floor(v);
    if (y == v) {
        n = (int)v;
        y = hcephes_yn(n, x);
        return (y);
    }
    t = HCEPHES_PI * v;
    y = (cos(t) * hcephes_jv(v, x) - hcephes_jv(-v, x)) / sin(t);
    return (y);
}
