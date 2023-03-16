// Fast windowing functions implemented in native C
#ifndef __WINDOWING_H__
#define __WINDOWING_H__

#include <math.h>
#include "hcephes.h"

typedef double (*func_t)(const double * const, int);
typedef double (*weighted_func_t)(const double* const, const double* const, int);

double fast_sum(const double* const x, int k)
{
    int i;
    double s = 0.0;

    for (i = 0; i < k; i++)
    {
        s += x[i];
    }

    return s;
}

double fast_product(const double* const x, int k)
{
    int i;
    double p = 1.0;

    for (i = 0; i < k; i++)
    {
        p *= x[i];
    }

    return p;
}


double fast_fishers_combined(const double* const x, int k)
{
    int i;
    double p, chi = 0.0;

    for (i = 0; i < k; i++)
    {
        chi += log(x[i]);
    }
    chi *= -2.0;
    
    p = hcephes_chdtrc( (double) 2.0 * k, chi);
    return p;
}

double fast_stouffers_z(const double* const x, int k)
{
    int i;
    double p, z;
    double s = 0.0;
    
    for (i = 0; i < k; ++i)
    {
        s += hcephes_ndtri(1.0 - x[i]);
    }
    z = s / sqrt((double)k);
    p = hcephes_ndtr(-z);

    return p;
}

double* fast_windowing_func(const double* const x, int n, int hw, func_t func_ptr)
{
    int i, k = 2 * hw + 1;
    double* tmp  = (double*) malloc(k * sizeof(double));
    double* res = (double*) calloc(n, sizeof(double));

    for (i = hw; i < n-hw; ++i)
    {
        memcpy(&tmp[0], &x[i-hw], k * sizeof(double));
        res[i] = func_ptr(tmp, k);
    }

    free(tmp);

    return res;
}

double fast_weighted_stouffers_z(const double* const x, const double* const w, int k)
{
    int i;
    double p, z;
    double s = 0.0;
    double sw = 0.0;

    for (i = 0; i < k; ++i)
    {
        s += w[i] * hcephes_ndtri(1.0 - x[i]);
        sw += w[i] * w[i];
    }
    z = s / sqrt(sw);
    p = hcephes_ndtr(-z);

    return p;	
}

double* fast_weighted_windowing_func(const double* const x, const double* const w, int n, int hw, weighted_func_t func_ptr)
{
    int i, k = 2 * hw + 1;
    double* xtmp = (double*) malloc(k * sizeof(double));
    double* wtmp = (double*) malloc(k * sizeof(double));
    double* res = (double*) calloc(n, sizeof(double));

    for (i = hw; i < n-hw; ++i)
    {
        memcpy(&xtmp[0], &x[i-hw], k * sizeof(double));
        memcpy(&wtmp[0], &w[i-hw], k * sizeof(double));
        res[i] = func_ptr(xtmp, wtmp, k);
    }

    free(xtmp);
    free(wtmp);

    return res;
}

#endif