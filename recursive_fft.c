#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "common_defs.h"

void recursive_fft(complex* a, complex* y, int n, int inv)
{
   complex w, wn, twiddle;
   complex* a0;
   complex* a1;
   complex* y0;
   complex* y1;

   /* Base Case */
   if (n == 1)
   {
      y[0] = a[0];
      return;
   }

   /* Calculating principal nth root of unity (i.e. exp(2*M_PI*i/n)) */
   if (inv)
   {
      wn.r = cos(-2*M_PI/(double)n);
      wn.i = sin(-2*M_PI/(double)n);
   }
   else
   {
      wn.r = cos(2*M_PI/(double)n);
      wn.i = sin(2*M_PI/(double)n); 
   }
   w.r = 1.0;
   w.i = 0.0;
   
   /* allocating memory for even/odd coefficients and corresponding FFTs */
   a0 = (complex*)malloc((n/2) * sizeof(complex));
   a1 = (complex*)malloc((n/2) * sizeof(complex));
   y0 = (complex*)malloc((n/2) * sizeof(complex));
   y1 = (complex*)malloc((n/2) * sizeof(complex));
   
   /* Extracting even and odd coefficients */
   for (int i = 0; i < (n/2); i++)
   {
      a0[i] = a[2*i];
      a1[i] = a[2*i+1];
   }

   /* Calculating 2 FFTs of size n/2 */
   recursive_fft(a0, y0, n/2, inv);
   recursive_fft(a1, y1, n/2, inv);
   
   /* Combining results from half-size FFTs */
   for (int k = 0; k < (n/2); k++)
   {
      twiddle  = complex_mul(w, y1[k]);
      y[k]     = complex_add(y0[k], twiddle);
      y[k+n/2] = complex_sub(y0[k], twiddle);
      
     
      w = complex_mul(w, wn);
   }
   
   free(a0);
   free(a1);
   free(y0);
   free(y1);
   
   return;
}
