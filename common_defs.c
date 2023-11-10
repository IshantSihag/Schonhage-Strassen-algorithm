#include <stdio.h>
#include <stdlib.h>
#include "common_defs.h"

complex complex_mul(complex a, complex b)
{
   complex ans;
   
   ans.r = a.r * b.r - a.i * b.i;
   ans.i = a.r * b.i + a.i * b.r;
   
   return ans;
}

complex complex_add(complex a, complex b)
{
   complex ans;
   
   ans.r = a.r + b.r;
   ans.i = a.i + b.i;
   
   return ans;
}

complex complex_sub(complex a, complex b)
{
   complex ans;
   
   ans.r = a.r - b.r;
   ans.i = a.i - b.i;
   
   return ans;
}


void print_poly(complex* a, int n)
{
   // printf("\nPrinting coefficients for x^k:\n");
   // for (int i = 0; i < (n); i++)
   // {
   //    printf("[%d] = %.0f\n", i, a[i].r);
   // }
   int flag=0;
   for(int i=n-1; i>=0; i--)
   {
      // printf("x^%d=%lf\n", i, a[i].r);
      // if(flag==0)
      //    continue;
      if(a[i].r<1e-10)
         continue;
      // flag=1;
      if(i)
      {
         if(flag==0)
            printf("%lgx^%d ", a[i].r, i);
         else
            printf("%+lgx^%d ", a[i].r, i);

      }
      else
         printf("%+lg\n", a[i]);
      flag=1;
   }
   printf("\n");
}

void get_poly(complex* a, int n)
{
   if(n==0)
      printf("Polynomial with 0 degree is not valid\n");
   double coeff;
   int next_power_of_2 = 1;
   while (next_power_of_2 < n)
      next_power_of_2 <<= 1;

   /* Read coefficients from stdin */
   printf("Enter the coefficient of :\n");
   for (int i = n-1; i >=0; i--)
   {
      printf("x^%d : ",i);
      scanf("%lf",&coeff);
      a[i].r = coeff;
      a[i].i = 0.0;
   }
   
   /* Pad the rest with zeros */
   for (int i = n; i < (2 * next_power_of_2); i++)
   {
      a[i].r = 0.0; a[i].i = 0.0;
   }
}

void poly_mul(complex* a, complex* b, int n)
{
   complex* ya;
   complex* yb;
   int j;
   
   /* Allocate storage for fft results */
   ya = (complex*)malloc(n * sizeof(complex));
   yb = (complex*)malloc(n * sizeof(complex));
   
   /* DFT of A and B */
   recursive_fft(a, ya, n, 0);
   recursive_fft(b, yb, n, 0);
   
   /* Pointwise Multiplication */
   for (j = 0; j < n; j++)
      ya[j] = complex_mul(ya[j], yb[j]);
      
   /* Inverse DFT (swapped input and output arrays) */
   recursive_fft(ya, a, n, 1);
   
   /* Divide real part by n */
   for (j = 0; j < (n-1); j++)
      a[j].r = a[j].r/n;
      
   free(ya);
   free(yb);
}


void poly_add(complex*sum, complex* a, complex* b, int n)
{
   /* Pointwise Multiplication */
   for (int j = 0; j < n; j++)
      sum[j] = complex_add(a[j], b[j]);
      
}

void poly_sub(complex*sum, complex* a, complex* b, int n)
{
   /* Pointwise Multiplication */
   for (int j = 0; j < n; j++)
      sum[j] = complex_sub(a[j], b[j]);
      
}

void poly_differentiate(complex*a, int n)
{
   for(int i=0; i<n-1; i++)
      a[i].r = a[i+1].r*(i+1);
   a[n-1].r=0;
}

void poly_integrate(complex*a, int n)
{
   for(int i=n; i>0; i--)
      a[i].r=(a[i-1].r)/i;
   // printf("%lf \n %lg", a[2].r, a[2].r);
}


void linear(double a0,double a1)
{
    double e=-a0/a1;
    printf("%.3lf, ",e);
}

void quadratic(double t,double r,double s)
{
    double det,x1,x2,x1r,x1i;
    det=(r*r)-(4*s*t);
    if(det>=0)
    {
        x1= (-r+sqrt(det))/2*t;
        x2= (-r-sqrt(det))/2*t;
        printf("%.3lf, %.3lf, ",x2,x1);
    }
    else
    {
        x1r=-r/2;
        x1i=sqrt(fabs(det))/2;
        printf("%.3lf + %.3lf j, %.3lf - %.3lf j, ",x1r,x1i,x1r,x1i);
    }
}

void roots(double*a, int n)
{
   int w,j,i;
   double t,det,p,q,x1=0.00,x2=0.00,x1r,x1i,r=0.1,s=0.1,ds,dr,b[100],c[100],d[100];
   b[4]=0,b[3]=0;
   
   w=n;
   if(w==1) // Linear Case
   {
      linear(a[0],a[1]);
      w--;
   }
   if(w==2) // Quadratic Case
   {
      quadratic(a[2],a[1],a[0]);
      w=w-2;
   }
   while(w>=3) // Bairstow's Case
   {
      for(j=1;j<=50;j++) {
         b[n]=a[n];
         b[n-1]=a[n-1]-r*b[n];
         for(i=n-2;i>=1;i--)
         {
               b[i]=a[i]-r*b[i+1]-s*b[i+2];
         }
         b[0]=a[0]-s*b[2];
         c[n]=b[n];
         c[n-1]=b[n-1]-r*c[n];
         for(i=n-2;i>=2;i--)
         {
               c[i]=b[i]-r*c[i+1]-s*c[i+2];
         }
         c[1]=-s*c[3];
         d[n]=b[n];
         d[n-1]=b[n-1]-r*d[n];
         for(i=n-2;i>=3;i--)
         {
               d[i]=b[i]-r*d[i+1]-s*d[i+2];
         }
         d[2]=b[2]-s*d[4];
         dr=(b[0]*d[3]-b[1]*d[2])/((d[2]-r*d[3])*d[2]+s*d[3]*d[3]);
         ds=(-b[1]*s*d[3]-b[0]*(d[2]-r*d[3]))/((d[2]-r*d[3])*d[2]+s*d[3]*d[3]);
         p=r-dr;
         q=s-ds;
         r=p;
         s=q;
      }
      t=1;
      quadratic(t,r,s);
      w=w-2;
      for(i=n;i>=0;i--){
         a[n-i]=b[n-i];
      }
      for(i=n;i>=0;i--) {
         a[n-i]=a[n-i+2];
      }
   }
   if(w==2) // Quadratic Case
   {
      quadratic(b[4],b[3],b[2]);
      w=w-2;
   }
   if(w==1) // Linear Case
   {
      linear(b[2],b[3]);
      w--;
   }
}
