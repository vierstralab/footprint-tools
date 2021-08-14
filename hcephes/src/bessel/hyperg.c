#include "hcephes.h"

static double hcephes_hy1f1p(double a, double b, double x, double *err);
static double hcephes_hy1f1a(double a, double b, double x, double *err);

HCEPHES_API double hcephes_hyperg(double a, double b, double x) {
    double asum, psum, acanc, pcanc, temp;
    pcanc = 0.0;
    acanc = 0.0;

    /* See if a Kummer transformation will help */
    temp = b - a;
    if (fabs(temp) < 0.001 * fabs(a))
        return (exp(x) * hcephes_hyperg(temp, b, -x));

    psum = hcephes_hy1f1p(a, b, x, &pcanc);
    if (pcanc < 1.0e-15)
        goto done;

    /* try asymptotic series */

    asum = hcephes_hy1f1a(a, b, x, &acanc);

    /* Pick the result with less estimated error */

    if (acanc < pcanc) {
        pcanc = acanc;
        psum = asum;
    }

done:
    if (pcanc > 1.0e-12)
        hcephes_mtherr("hyperg", HCEPHES_PLOSS);

    return (psum);
}

/* Power series summation for confluent hypergeometric function		*/

static double hcephes_hy1f1p(double a, double b, double x, double *err) {
    double n, a0, sum, t, u, temp;
    double an, bn, maxt, pcanc;

    /* set up for power series summation */
    an = a;
    bn = b;
    a0 = 1.0;
    sum = 1.0;
    n = 1.0;
    t = 1.0;
    maxt = 0.0;

    while (t > HCEPHES_MACHEP) {
        if (bn == 0) /* check bn first since if both	*/
        {
            hcephes_mtherr("hyperg", HCEPHES_SING);
            return (HUGE_VAL); /* an and bn are zero it is	*/
        }
        if (an == 0) /* a singularity		*/
            return (sum);
        if (n > 200)
            goto pdone;
        u = x * (an / (bn * n));

        /* check for blowup */
        temp = fabs(u);
        if ((temp > 1.0) && (maxt > (HUGE_VAL / temp))) {
            pcanc = 1.0; /* estimate 100% error */
            goto blowup;
        }

        a0 *= u;
        sum += a0;
        t = fabs(a0);
        if (t > maxt)
            maxt = t;
        /*
                if( (maxt/fabs(sum)) > 1.0e17 )
                        {
                        pcanc = 1.0;
                        goto blowup;
                        }
        */
        an += 1.0;
        bn += 1.0;
        n += 1.0;
    }

pdone:

    /* estimate error due to roundoff and cancellation */
    if (sum != 0.0)
        maxt /= fabs(sum);
    maxt *= HCEPHES_MACHEP; /* this way avoids multiply overflow */
    pcanc = fabs(HCEPHES_MACHEP * n + maxt);

blowup:

    *err = pcanc;

    return (sum);
}

/*							hcephes_hy1f1a()
 */
/* asymptotic formula for hypergeometric function:
 *
 *        (    -a
 *  --    ( |z|
 * |  (b) ( -------- 2f0( a, 1+a-b, -1/x )
 *        (  --
 *        ( |  (b-a)
 *
 *
 *                                x    a-b                     )
 *                               e  |x|                        )
 *                             + -------- 2f0( b-a, 1-a, 1/x ) )
 *                                --                           )
 *                               |  (a)                        )
 */

static double hcephes_hy1f1a(double a, double b, double x, double *err) {
    double h1, h2, t, u, temp, acanc, asum, err1, err2;

    if (x == 0) {
        acanc = 1.0;
        asum = HUGE_VAL;
        goto adone;
    }
    temp = log(fabs(x));
    t = x + temp * (a - b);
    u = -temp * a;

    if (b > 0) {
        temp = hcephes_lgam(b);
        t += temp;
        u += temp;
    }

    h1 = hcephes_hyp2f0(a, a - b + 1, -1.0 / x, 1, &err1);

    temp = exp(u) / hcephes_gamma(b - a);
    h1 *= temp;
    err1 *= temp;

    h2 = hcephes_hyp2f0(b - a, 1.0 - a, 1.0 / x, 2, &err2);

    if (a < 0)
        temp = exp(t) / hcephes_gamma(a);
    else
        temp = exp(t - hcephes_lgam(a));

    h2 *= temp;
    err2 *= temp;

    if (x < 0.0)
        asum = h1;
    else
        asum = h2;

    acanc = fabs(err1) + fabs(err2);

    if (b < 0) {
        temp = hcephes_gamma(b);
        asum *= temp;
        acanc *= fabs(temp);
    }

    if (asum != 0.0)
        acanc /= fabs(asum);

    acanc *= 30.0; /* fudge factor, since error of asymptotic formula
                    * often seems this much larger than advertised */

adone:

    *err = acanc;
    return (asum);
}

/*							hcephes_hyp2f0()
 */

HCEPHES_API double hcephes_hyp2f0(double a, double b, double x, int type, double *err) {
    double a0, alast, t, tlast, maxt;
    double n, an, bn, u, sum, temp;

    an = a;
    bn = b;
    a0 = 1.0e0;
    alast = 1.0e0;
    sum = 0.0;
    n = 1.0e0;
    t = 1.0e0;
    tlast = 1.0e9;
    maxt = 0.0;

    do {
        if (an == 0)
            goto pdone;
        if (bn == 0)
            goto pdone;

        u = an * (bn * x / n);

        /* check for blowup */
        temp = fabs(u);
        if ((temp > 1.0) && (maxt > (HUGE_VAL / temp)))
            goto error;

        a0 *= u;
        t = fabs(a0);

        /* terminating condition for asymptotic series */
        if (t > tlast)
            goto ndone;

        tlast = t;
        sum += alast; /* the sum is one term behind */
        alast = a0;

        if (n > 200)
            goto ndone;

        an += 1.0e0;
        bn += 1.0e0;
        n += 1.0e0;
        if (t > maxt)
            maxt = t;
    } while (t > HCEPHES_MACHEP);

pdone: /* series converged! */

    /* estimate error due to roundoff and cancellation */
    *err = fabs(HCEPHES_MACHEP * (n + maxt));

    alast = a0;
    goto done;

ndone: /* series did not converge */

    /* The following "Converging factors" are supposed to improve accuracy,
     * but do not actually seem to accomplish very much. */

    n -= 1.0;
    x = 1.0 / x;

    switch (type) /* "type" given as subroutine argument */
    {
    case 1:
        alast *= (0.5 + (0.125 + 0.25 * b - 0.5 * a + 0.25 * x - 0.25 * n) / x);
        break;

    case 2:
        alast *= 2.0 / 3.0 - b + 2.0 * a + x - n;
        break;

    default:;
    }

    /* estimate error due to roundoff, cancellation, and nonconvergence */
    *err = HCEPHES_MACHEP * (n + maxt) + fabs(a0);

done:
    sum += alast;
    return (sum);

/* series blew up: */
error:
    *err = HUGE_VAL;
    hcephes_mtherr("hyperg", HCEPHES_TLOSS);
    return (sum);
}
