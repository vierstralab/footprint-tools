
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

result_t* fast_predict(double const *obs, double const *probs, int l, int half_window_width)
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

	// Smooth windows
	double *smoothed_win_counts = windowed_trimmed_mean(&win_counts[0], l, 50, 0.05);
	memcpy(&win_counts[0], &smoothed_win_counts[0], l*sizeof(double));
	free(smoothed_win_counts);

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
