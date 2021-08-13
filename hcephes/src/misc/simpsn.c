#include "hcephes.h"

#define NCOTE 8

/* 8th order formula */
static double simcon[] = {
    3.488536155202821869E-2,  2.076895943562610229E-1,   -3.27336860670194003527E-2,
    3.7022927689594356261E-1, -1.6014109347442680776E-1,
};

HCEPHES_API double hcephes_simpsn(double f[],  /* tabulated function */
                                  double delta /* spacing of arguments */
) {
    extern double simcon[];
    double ans;
    int i;

    ans = simcon[NCOTE / 2] * f[NCOTE / 2];
    for (i = 0; i < NCOTE / 2; i++)
        ans += simcon[i] * (f[i] + f[NCOTE - i]);

    return (ans * delta * NCOTE);
}
