#ifndef __CEPHES_H__
#define __CEPHES_H__

/*
Declaration for functions called (in)directly by Cython libraries

2015-11-18: Initial commit and function definitions
*/

// calls to some of cephes library overlap math.h and 
// cause erratic behavoir when used with cython and numpy

#define log		c_log
#define round	c_round
#define sqrt 	c_sqrt
#define gamma	c_gamma
#define incbet	c_incbet
#define chdtr	c_chdtr
#define chdtrc	c_chdtrc
#define ndtr 	c_ndtr
#define ndtri	c_ndtri

//math
extern double log ( double );
extern double round ( double );
extern double sqrt ( double );

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
extern double chdtrc( double, double); //complemented

// normal
extern double ndtr( double );
extern double ndtri( double );

#endif
