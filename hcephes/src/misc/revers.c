#include "hcephes.h"
#include <stdlib.h>

extern int MAXPOL; /* initialized by polini() */

HCEPHES_API void hcephes_revers(double y[], double x[], int n) {
    double *yn, *yp, *ysum;
    int j;

    if (y[1] == 0.0)
        hcephes_mtherr("revers", HCEPHES_DOMAIN);

    j = (int)((unsigned long)(MAXPOL + 1) * sizeof(double));
    yn = (double *)malloc((size_t)j);
    yp = (double *)malloc((size_t)j);
    ysum = (double *)malloc((size_t)j);

    hcephes_polmov(y, n, yn);
    hcephes_polclr(ysum, n);
    x[0] = 0.0;
    x[1] = 1.0 / y[1];
    for (j = 2; j <= n; j++) {
        /* A_(j-1) times the expansion of y^(j-1)  */
        hcephes_polmul(&x[j - 1], 0, yn, n, yp);
        /* The expansion of the sum of A_k y^k up to k=j-1 */
        hcephes_poladd(yp, n, ysum, n, ysum);
        /* The expansion of y^j */
        hcephes_polmul(yn, n, y, n, yn);
        /* The coefficient A_j to make the sum up to k=j equal to zero */
        x[j] = -ysum[j] / yn[j];
    }
    free(yn);
    free(yp);
    free(ysum);
}
