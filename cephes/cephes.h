#ifndef __CEPHES_H__
#define __CEPHES_H__

/*
Declaration for functions called (in)directly by Cython libraries

2015-11-18: Initial commit and function definitions
*/

// calls to some of cephes library overlap math.h and 
// cause erratic behavoir when used with cython and numpy

#define exp		c_exp
#define log		c_log
#define log1p	c_log1p

#define round	c_round

#define sqrt 	c_sqrt

#define gamma	c_gamma
#define lgam	c_lgamma

#define incbet	c_incbet
#define beta	c_beta
#define lbeta	c_lbeta

#define chdtr	c_chdtr
#define chdtrc	c_chdtrc

#define ndtr 	c_ndtr
#define ndtri	c_ndtri

#define psi		c_psi

//const
extern double PI;

//math
extern double exp ( double );
extern double log ( double );
extern double log1p ( double ); // unity.c
extern double round ( double );
extern double sqrt ( double );

//special functions
extern double polevl ( double, void *, int );
extern double p1evl ( double, void *, int );

extern double beta( double, double);
extern double lbeta( double, double);

extern double gamma ( double );
extern double lgam( double );
extern double igamc ( double, double );
extern double igam ( double, double );
extern double igami ( double, double );

extern double incbet( double, double, double );

// chi-sq   
extern double chdtr( double, double );
extern double chdtrc( double, double ); //complemented

// normal
extern double ndtr( double );
extern double ndtri( double );

// psi
extern double psi( double ); // psi.c

#endif
