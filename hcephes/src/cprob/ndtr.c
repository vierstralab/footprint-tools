#include "hcephes.h"

/* Define this macro to suppress error propagation in exp(x^2)
   by using the hcephes_expx2 function.  The tradeoff is that doing so
   generates two calls to the exponential function instead of one.  */
#define USE_EXPXSQ 1

static double P[] = {
    2.46196981473530512524E-10, 5.64189564831068821977E-1, 7.46321056442269912687E0,
    4.86371970985681366614E1,   1.96520832956077098242E2,  5.26445194995477358631E2,
    9.34528527171957607540E2,   1.02755188689515710272E3,  5.57535335369399327526E2};
static double Q[] = {
    /* 1.00000000000000000000E0,*/
    1.32281951154744992508E1, 8.67072140885989742329E1, 3.54937778887819891062E2,
    9.75708501743205489753E2, 1.82390916687909736289E3, 2.24633760818710981792E3,
    1.65666309194161350182E3, 5.57535340817727675546E2};
static double R[] = {5.64189583547755073984E-1, 1.27536670759978104416E0,
                     5.01905042251180477414E0,  6.16021097993053585195E0,
                     7.40974269950448939160E0,  2.97886665372100240670E0};
static double S[] = {
    /* 1.00000000000000000000E0,*/
    2.26052863220117276590E0, 9.39603524938001434673E0, 1.20489539808096656605E1,
    1.70814450747565897222E1, 9.60896809063285878198E0, 3.36907645100081516050E0};
static double T[] = {9.60497373987051638749E0, 9.00260197203842689217E1,
                     2.23200534594684319226E3, 7.00332514112805075473E3,
                     5.55923013010394962768E4};
static double U[] = {
    /* 1.00000000000000000000E0,*/
    3.35617141647503099647E1, 5.21357949780152679795E2, 4.59432382970980127987E3,
    2.26290000613890934246E4, 4.92673942608635921086E4};

#define UTHRESH 37.519379347

HCEPHES_API double hcephes_ndtr(double a) {
    double x, y, z;

    x = a * HCEPHES_SQRTH;
    z = fabs(x);

    /* if( z < HCEPHES_SQRTH ) */
    if (z < 1.0)
        y = 0.5 + 0.5 * hcephes_erf(x);

    else {
#ifdef USE_EXPXSQ
        /* See below for erfce. */
        y = 0.5 * hcephes_erfce(z);
        /* Multiply by exp(-x^2 / 2)  */
        z = hcephes_expx2(a, -1);
        y = y * sqrt(z);
#else
        y = 0.5 * hcephes_erfc(z);
#endif
        if (x > 0)
            y = 1.0 - y;
    }

    return (y);
}

/* Exponentially scaled erfc function
   exp(x^2) erfc(x)
   valid for x > 1.
   Use with ndtr and hcephes_expx2.  */
HCEPHES_API double hcephes_erfce(double x) {
    double p, q;

    if (x < 8.0) {
        p = hcephes_polevl(x, P, 8);
        q = hcephes_p1evl(x, Q, 8);
    } else {
        p = hcephes_polevl(x, R, 5);
        q = hcephes_p1evl(x, S, 6);
    }
    return (p / q);
}

HCEPHES_API double hcephes_erf(double x) {
    double y, z;

    if (fabs(x) > 1.0)
        return (1.0 - hcephes_erfc(x));
    z = x * x;
    y = x * hcephes_polevl(z, T, 4) / hcephes_p1evl(z, U, 5);
    return (y);
}

HCEPHES_API double hcephes_erfc(double a) {
    double p, q, x, y, z;

    if (a < 0.0)
        x = -a;
    else
        x = a;

    if (x < 1.0)
        return (1.0 - hcephes_erf(a));

    z = -a * a;

    if (z < -HCEPHES_MAXLOG) {
    under:
        hcephes_mtherr("erfc", HCEPHES_UNDERFLOW);
        if (a < 0)
            return (2.0);
        else
            return (0.0);
    }

#ifdef USE_EXPXSQ
    /* Compute z = exp(z).  */
    z = hcephes_expx2(a, -1);
#else
    z = exp(z);
#endif
    if (x < 8.0) {
        p = hcephes_polevl(x, P, 8);
        q = hcephes_p1evl(x, Q, 8);
    } else {
        p = hcephes_polevl(x, R, 5);
        q = hcephes_p1evl(x, S, 6);
    }
    y = (z * p) / q;

    if (a < 0)
        y = 2.0 - y;

    if (y == 0.0)
        goto under;

    return (y);
}
