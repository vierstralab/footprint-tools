#include "hcephes.h"

HCEPHES_API double hcephes_btdtr(double a, double b, double x) {
    return (hcephes_incbet(a, b, x));
}
