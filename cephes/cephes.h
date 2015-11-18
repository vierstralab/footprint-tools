#ifndef __CEPHES_H__
#define __CEPHES_H__

/*
Declaration for functions called (in)directly by Cython libraries

2015-11-18: Initial commit and function definitions
*/

// calls to some of cephes library overlap math.h
// 
#define gamma c_gamma

//special functions
extern double polevl ( double, void *, int );
extern double p1evl ( double, void *, int );

extern double beta(double a, double b);
extern double lbeta(double a, double b);

extern double gamma ( double );
extern double igamc ( double, double );
extern double igam ( double, double );
extern double igami ( double, double );

extern double incbet( double, double, double );

// chi-sq   
extern double chdtr( double, double );

// normal
extern double ndtr( double );
extern double ndtri( double );

#endif
