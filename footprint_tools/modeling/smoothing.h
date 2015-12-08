// Copyright 2015 Jeff Vierstra

#ifndef __SMOOTHING_H_
#define __SMOOTHING_H_

/*
from Numerical Recipes
*/
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

double quickselect(double* arr, int n, int k)
{
    unsigned long i, ir, j, l, mid;
    double a, temp;

    l=0;
    ir=n-1;
    for(;;) {
        if (ir <= l+1) { 
            if (ir == l+1 && arr[ir] < arr[l]) {
                SWAP(arr[l],arr[ir]);
            }
            return arr[k];
        }
        else {
            mid=(l+ir) >> 1; 
            SWAP(arr[mid],arr[l+1]);
            if (arr[l] > arr[ir]) {
                SWAP(arr[l],arr[ir]);
            }
            if (arr[l+1] > arr[ir]) {
                SWAP(arr[l+1],arr[ir]);
            }
            if (arr[l] > arr[l+1]) {
                SWAP(arr[l],arr[l+1]);
            }
            i=l+1; 
            j=ir;
            a=arr[l+1]; 
            for (;;)
            { 
                do i++; while (arr[i] < a); 
                do j--; while (arr[j] > a); 
                if (j < i) break; 
                SWAP(arr[i],arr[j]);
            } 
            arr[l+1]=arr[j]; 
            arr[j]=a;
            if (j >= k) ir=j-1; 
            if (j <= k) l=i;
        }
    }
}


/*
from "Fast Computation of Trimmed Means" (Gleb Beliakov)
*/
double weighted(double x, double t1, double t2, double w1, double w2)
{
    if(x < t2 && x > t1)
        return x;
    if(x < t1)
        return 0;
    if(x > t2)
        return 0;
    if(x == t1) 
        return w1*x;
    return w2*x;
}

double trimmed_sum(double x[], int n, int k)
{
    int i;
    double w1, w2, OS1, OS2;
    double t = 0;

    OS1 = quickselect(x, n, k);
    OS2 = quickselect(x, n, n-k-1);

    // compute weights
    double a, b=0, c, d=0, dm=0, bm=0, r;

    for(i=0; i<n; i++)
    {
        r = x[i];
        if(r < OS1) bm += 1;
        else if(r == OS1) b += 1;
        if(r < OS2) dm += 1;
        else if(r == OS2) d += 1;
    }
    a = b + bm - k;
    c = n - k - dm;
    w1 = a/b;
    w2 = c/d;
    for(i=0; i<n; i++)
        t += weighted(x[i], OS1, OS2, w1, w2);
    return t;
}

double trimmed_mean(double x[], int n, int k)
{
    return trimmed_sum(x, n, k) / (n-2*k);
}


double* windowed_trimmed_mean(double* win, int l, int half_window_width, double clip)
{
    int i, j, k;
    int w = (half_window_width * 2) + 1;

    k = (int)( (double)w * clip );

    double *tmp = malloc(w * sizeof(double));
    double *res = calloc(l, sizeof(double));
    
    for (i = half_window_width; i < l-half_window_width; i++)
    {
        // index to start
        j = i - half_window_width;

        // copy data into new memory location
        memcpy(&tmp[0], &win[j], w * sizeof(double));
        
        // trimmed mean
        res[i] = trimmed_mean(tmp, w, k);

    }

    free(tmp);

    return res;
}

#endif