#include "hcephes.h"

HCEPHES_API double hcephes_chbevl(double x, double array[], int n) {
    double b0, b1, b2, *p;
    int i;

    p = array;
    b0 = *p++;
    b1 = 0.0;
    i = n - 1;

    do {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + *p++;
    } while (--i);

    return (0.5 * (b0 - b2));
}
