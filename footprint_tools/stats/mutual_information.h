// Copyright 2016 Jeff Vierstra

#ifndef __MUTUAL_INFORMATION_H__
#define __MUTUAL_INFORMATION_H__

#include <stdio.h>

#include "cephes.h"

extern double c_log( double );

double* joint_probability(int* x, int* y, int n)
{
	int i, k;
	double p[4] = {0.0};

	for(i = 0; i < n; ++i) {
		k = x[i];
		k ^= (-y[i] ^ k) & (1 << 1);
		p[k] += 1.0;
	}

	for(k = 0; k < 4; ++k)
	{
		p[k] /= n;
	}

	return &p[0];
}

double entropy(double p)
{
	double h = 0.0;

	if(p > 0.0 && p < 1.0)
	{
		h += p * c_log(p);
		h += (1.0-p) * c_log(1.0-p);
	}

	return -h;
}

double joint_entropy(double* p)
{
	int k;
	double h = 0.0;

	for(k = 0; k < 4; ++k)
	{
		if(p[k] > 0.0 && p[k] < 1.0)
		{
			h += p[k] * c_log(p[k]);
		}
	}

	return -h;
}


#endif