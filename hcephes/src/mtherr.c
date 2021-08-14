#include "hcephes.h"

#include <stdio.h>

int merror = 0;

static char *ermsg[7] = {
    "unknown",     /* error code 0 */
    "domain",      /* error code 1 */
    "singularity", /* et seq.      */
    "overflow",    "underflow", "total loss of precision", "partial loss of precision"};

HCEPHES_API int hcephes_mtherr(char *name, int code) {

    //printf("\n%s ", name);

    merror = code;

    if ((code <= 0) || (code >= 7))
        code = 0;
    //printf("%s error\n", ermsg[code]);

    return 0;
}
