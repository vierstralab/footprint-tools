#ifndef __MUTUAL_INFORMATION_H__
#define __MUTUAL_INFORMATION_H__

#include <stdio.h>

#include "cephes.h"

extern double c_log( double );

int pairwise_mutual_information(int* mat, int nx, int ny, double *result)
{
	//
	int i, j, k;

	//histograms
	size_t n_bins = 4;

	double *hx, *hy, *hxy;
	hx = (double*) calloc(n_bins, sizeof(double));
	hy = (double*) calloc(n_bins, sizeof(double));
	hxy = (double*) calloc(n_bins * n_bins, sizeof(double));

	assert(hx == NULL);
	assert(hy == NULL);
	assert(hxy == NULL);

	int hi, hj;

	double pxy, px, py;
	
	//for(i = 0; i < nx*ny; i++)
	//{
	//	fprintf(stderr, "%d\t%d\n", i, mat[i]);
	//}
	//return 0;

	for(i = 0; i < ny; i++)
	{

		for(j = 0; j < ny; j++)
		{
			//fprintf(stderr, "Cols: %i\t%i\n", i, j );

			//clear histogram
			for(hi = 0; hi < n_bins*n_bins; hi++)
			{
				if(hi < n_bins)
				{
					hx[hi] = hy[hi] = 0.0;
				}
				hxy[hi] = 0.0;
			}

			for(k = 0; k < nx; k++)
			{
				hi = mat[(k * ny) + i];
				//if(hi < 0) { fprintf(stderr, "hi too small"); }
				//if(hi >= n_bins) { fprintf(stderr, "hi too big"); }

				hj = mat[(k * ny) + j];
				//if(hj < 0) { fprintf(stderr, "hj too small"); }
				//if(hj >= n_bins) { fprintf(stderr, "hj too big"); }

				//fprintf(stderr, "%d\t%d\n", hi, hj );

				hx[hi] += 1.0 / (double)nx;
				hy[hj] += 1.0 / (double)nx;
				hxy[(hi * n_bins) + hj] += 1.0 / (double)nx;
			}

			for(hi = 0; hi < n_bins; hi++)
			{
				for(hj = 0; hj < n_bins; hj++)
				{
					pxy = hxy[(hi * n_bins) + hj];
					px = hx[hi];
					py = hy[hj];

					//fprintf(stderr, "\t%0.2f\t%0.2f\t%0.2f\n", px, py, pxy);

					if(pxy > 0 && px > 0 && py > 0)
					{
						result[(i * ny) + j] += pxy * c_log((pxy)/(px*py)) ;
					}
				}
				//fprintf(stderr, "\n");
			
			}
		}
	}
	
	free(hx);
	free(hy);
	free(hxy);

	return 0;
}

double pairwise_mutual_information_0(int* x, int* y, int n)
{

	//histograms
	size_t n_bins = 4;
	double *hx, *hy, *hxy;
	
	hx = (double*) calloc(n_bins, sizeof(double));
	hy = (double*) calloc(n_bins, sizeof(double));
	hxy = (double*) calloc(n_bins * n_bins, sizeof(double));

	int i, j, xi, yi;

	//clear histogram
	for(i = 0; i < n_bins*n_bins; i++)
	{
		if(i < n_bins)
		{
			hx[i] = hy[i] = 0.0;
		}
		hxy[i] = 0.0;
	}

	// for "row", assign to respective histogram bins

	//fprintf(stderr, "%d\n", nx);

	for(j = 0; j < n; j++)
	{
		xi = x[j];
		yi = y[j];

		hx[xi] += 1.0;
		hy[yi] += 1.0;
		hxy[(xi * n_bins) + yi] += 1.0;
	}

	double px, py, pxy, pmi, mi = 0.0;

	for(i = 0; i < n_bins; i++)
	{
		for(j = 0; j < n_bins; j++)
		{
			px = hx[i] / (double)n;
			py = hy[j] / (double)n;
			pxy = hxy[(j * n_bins) + i] / (double)n;

			fprintf(stderr, "%d\t%d", xi, yi);
			fprintf(stderr, "\t%0.4f\t%0.4f\t%0.4f", px, py, pxy);

			if(px > 0 && py > 0 && pxy > 0)
			{
				pmi = pxy * c_log(pxy/(px*py));
				fprintf(stderr, "\t%0.4f\n", pmi);

				mi += pmi;
			}
			fprintf(stderr, "\n");

			//fprintf(stderr, "%0.2f\t", hxy[(xi * n_bins) + yi]);
		}
		//fprintf(stderr, "\n");
	}

	
	free(hx);
	free(hy);
	free(hxy);

	return mi;
}


#endif