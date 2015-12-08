// Copyright 2015 Jeff Vierstra

#ifndef __PREDICT_H_
#define __PREDICT_H_

#include "cephes.h"
#include "smoothing.h"

typedef struct result
{
	double *exp;
	double *win;
} result_t;

void free_result_t(result_t* res)
{
	free(res->exp);
	free(res->win);
	free(res);
}

result_t* fast_predict(double const *obs, double const *probs, int l, int half_window_width, int smoothing_half_window_width, double smoothing_clip)
{
	//probs = l
	int i, j;
	//int w = (half_window_width * 2);

	// Create memory space for windows
	double *win_counts	= (double*) calloc(l, sizeof(double));
	double *win_probs 	= (double*) calloc(l, sizeof(double));
	double *exp_counts	= (double*) calloc(l, sizeof(double));

	// Make windowed counts and probabilities
	for (i = half_window_width; i < l-half_window_width; i++)
	{
		for (j = -half_window_width; j < half_window_width; j++)
		{
			win_counts[i] += obs[i+j];
			win_probs[i] += probs[i+j]; 
		}
	}

	if (smoothing_half_window_width > 0)
	{
		double *smoothed_win_counts = windowed_trimmed_mean(&win_counts[0], l, smoothing_half_window_width, smoothing_clip);
		// free
		free(win_counts);
		// set
		win_counts = smoothed_win_counts;
	}

	// Compute expected
	for (i = half_window_width; i < l-half_window_width; i++)
	{
		exp_counts[i] = c_round( (probs[i] / win_probs[i]) * win_counts[i] );
	}

	result_t *res = (result_t*) malloc(sizeof(result_t));
	res->exp = exp_counts;
	res->win = win_counts;

	// Clean up
	free(win_probs);

	// Return
	return res;
}

#endif