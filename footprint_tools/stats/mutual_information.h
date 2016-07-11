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

#endif