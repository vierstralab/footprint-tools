#ifndef __MUTUAL_INFORMATION_H__
#define __MUTUAL_INFORMATION_H__

#include <stdio.h>
#include <stdlib.h>

#include "cephes.h"

#include "kdtree.h"

extern double c_log( double );



/*
a, b, c are offset values to switch from a 2d to a 3d array

a*i + b*i + c

2D: a = N2, b = 1, c = 0 -> N2*i + i + 0
3D: a = N3*N2, b = N2, c = n3 -> N3*N2*i + N3*j + n3 

*/
int pairwise_mutual_information(double* mat, int nx, int ny, int n_bins, double *result, int a, int b, int c)
{
	//
	int i, j, k, offset;

	//histograms
	//size_t n_bins = 4;

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

		for(j = i+1; j < ny; j++)
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
						offset = (a*i) + (b*j) + c;
						result[offset] += pxy * c_log((pxy)/(px*py));
						//result[(i * ny) + j] += pxy * c_log((pxy)/(px*py)) ;
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

static int rand_int(int n)
{
	int limit = RAND_MAX - RAND_MAX % n;
	int rnd;

	do
	{
		rnd = rand();
	}
	while(rnd >= limit);

	return rnd % n;
}

void shuffle_columns_fast(double* mat, int nx, int ny)
{

	int i, j, s, ii, jj, tmp;

	for(j=0; j<ny; j++)
	{
		for(i=nx-1; i>0; i--)
		{
			s = rand_int(i+1);

			ii = (i*ny) + j;
			jj = (s*ny) + j;
			tmp = mat[jj];
			mat[jj] = mat[ii];
			mat[ii] = tmp;
		}
	}

}


unsigned int get_msec(void)
{
	static struct timeval timeval, first_timeval;

	gettimeofday(&timeval, 0);

	if(first_timeval.tv_sec == 0) {
		first_timeval = timeval;
		return 0;
	}
	return (timeval.tv_sec - first_timeval.tv_sec) * 1000 + (timeval.tv_usec - first_timeval.tv_usec) / 1000;
}

int pairwise_mutual_information_kraskov(double* mat, int m, int n, double *result, int a, int b, int c)
{

	int i, j, k;

	kdtree_t **kdarr = (void*) malloc(sizeof(kdtree_t *) * n);

	double pt[2];

	for(i = 0; i < n; i++)
	{
		kdarr[i] = kd_create(1);

		for(k = 0; k < m; k++)
		{
			pt[0] = mat[(k*n) + i];
			kd_insert(kdarr[i], &pt[0]);
		}
	}

	// clean up
	for(i = 0; i < n; i++)
	{
		free(kdarr[i]);
	}
	free(kdarr);

	return 0;

}


#endif