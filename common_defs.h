#ifndef COMMON_DEFS_H
#define COMMON_DEFS_H
#define _USE_MATH_DEFINES /* for M_PI constant */
#include <math.h>
// #define PI M_PI /* from math.h */

typedef struct
{
   long double r; /* real */
   long double i; /* imaginary */
} complex;

void bit_reverse_copy(complex* a, complex* a_rev_copy);

/*---------------------------------------------------------
** NAME: recursive_fft
**
** PURPOSE:
**    Implement the recursive FFT algorithm in CLRS for
**    evaluating polynomials at complex roots of unity.
**    Used as a helper function for multiplying polynomials.
**
** INPUTS:
**    a     Complex array of polynomial coefficients
**    n     Length of array (must be a power of 2)
**    inv   1 if performing inverse DFT, 0 otherwise
**
** OUTPUTS:
**    y     An array of n complex numbers representing the 
**          evaluation of the input polynomial at complex 
**          roots of unity. User must allocate this memory.  
**
** RETURNS: void
**-------------------------------------------------------*/
void recursive_fft(complex* a, complex* y, int n, int inv);

void print_poly(complex* a, int n);
void get_poly(complex* a, int n);
void poly_mul(complex* a, complex* b, int n);
void poly_add(complex*sum, complex* a, complex* b, int n);
void poly_sub(complex*sum, complex* a, complex* b, int n);
void poly_differentiate(complex*a, int n);
void poly_integrate(complex*a, int n);
void linear(double a0,double a1);
void quadratic(double t,double r,double s);
void roots(double*a, int n);

complex complex_mul(complex a, complex b);
complex complex_add(complex a, complex b);
complex complex_sub(complex a, complex b);
#endif