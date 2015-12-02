#include "cephes.h"

extern double c_log( double );
extern double c_sqrt( double );
extern double c_chdtrc( double, double );
extern double c_ndtri( double );
extern double c_ndtr( double );

typedef double (*func_t)(const double * const, int);

double fast_fishers_combined(const double* const x, int n)
{
	int i;
	double p, chi = 0.0;

	for (i = 0; i < n; i++)
	{
		chi += c_log(x[i]);
	}
	chi *= -2.0;
	p = c_chdtrc( (double) 2.0 * n, chi);
	return p;
}

double fast_stouffers_z(const double* const x, int n)
{
	int i;
	double p, z, s = 0.0;
	for (i = 0; i < n; ++i)
	{
		s += c_ndtri(1.0 - x[i]);
	}
	z = s / c_sqrt((double)n);
	p = c_ndtr(-z);

	return p;
}

double* fast_windowing_func(const double* const x, int n, int w, func_t func_ptr)
{
	int i, fw = 2*w + 1;
	double* tmp  = (double*) malloc(fw * sizeof(double));
	double* res = (double*) calloc(n, sizeof(double));

	for (i = w; i < n-w; ++i)
	{
		memcpy(&tmp[0], &x[i-w], fw * sizeof(double));
		res[i] = func_ptr(tmp, fw);
	}

	free(tmp);

	return res;
}