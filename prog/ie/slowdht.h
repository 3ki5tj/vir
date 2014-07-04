#ifndef SLOWDHT_H__
#define SLOWDHT_H__



/* Enhanced version of GSL's discrete Hankel transform
 * Features:
 *  o long double and __float128
 *  o using disk for large memory transforms */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "xdouble.h"
#include <gsl/gsl_sf_bessel.h>



#ifdef SLOWDHT

#define xdht slowdht
#define XDHT(f) slowdht_ ## f

#else

#include <gsl/gsl_dht.h>
#define xdht gsl_dht
#define XDHT(f) gsl_dht_ ## f
#define gsl_dht_newx(size, nu, xmax, usedisk) \
        gsl_dht_new(size, nu, xmax)

#endif



typedef struct slowdht_struct {
  size_t    size;
  xdouble   nu;
  xdouble   inu;
  xdouble   xmax;
  xdouble   kmax;
  xdouble   *j;    /* j_{nu, s} = j[s] */
  xdouble   *J2;   /* J_{nu+1}^2(j_m) */
  xdouble   *Jjj;  /* J_nu(j_i j_m / j_M) */
  char      Jjjfn[FILENAME_MAX];
  FILE      *Jjjfp;
  xdouble   *Jjjarr;
} slowdht;



char *slowdht_tmpdir = NULL; /* ended with / */
int slowdht_quiet = 0;
unsigned slowdht_block = 32; /* block size for reading/writing the temporary file */



/* J(nu, x) */
__inline xdouble besselJnu(xdouble nu, xdouble x)
{
  int inu = (int)(nu + .5);
  if ( FABS(nu - inu) < 1e-7 )
    return JN(inu, x);
  else
    return gsl_sf_bessel_Jnu((double) nu, (double) x);
}



/* return dJ(nu, x) / dx */
__inline xdouble besseldJnu(xdouble nu, xdouble x)
{
  int inu = (int) (nu + .5);
  xdouble x1, x2;
  if ( FABS(nu - inu) < 1e-7 ) {
    if ( inu == 0 ) return -J1(x);
    x1 = JN(inu - 1, x);
    x2 = JN(inu + 1, x);
  } else if ( FABS(nu - 0.5) < 1e-7 ) {
    x1 = COS(x);
    x2 = SIN(x)/x/2;
    return (x1 - x2) * SQRT(2/PI/x);
  } else {
    x1 = gsl_sf_bessel_Jnu((double) nu - 1, (double) x);
    x2 = gsl_sf_bessel_Jnu((double) nu + 1, (double) x);
  }
  return (x1 - x2) / 2;
}



/* compute the mth zero of the Bessel function
 * for xdouble == long double or __float128
 * use the GSL function get an estimate,
 * then refine it by Newton-Raphson's method  */
__inline xdouble besselzeroJnu(xdouble nu, int m)
{
  xdouble x, x1, y, yabs, dy, y1, y1abs;
  int i;

  x = gsl_sf_bessel_zero_Jnu((double) nu, m);
  /* the Newton-Raphson's method fails at x = 0 */
  if ( FABS(x) < 1e-8 ) return x;
  y = besselJnu(nu, x);
  yabs = FABS(y);
  for ( i = 0; ; i++ ) {
    dy = besseldJnu(nu, x);
    x1 = x - y/dy;
    y1 = besselJnu(nu, x1);
    y1abs = FABS(y1);
    if ( y1abs >= yabs ) break;
    x = x1;
    y = y1;
    yabs = y1abs;
  }
  //printf("iter %d, x %.20" XDBLPRNF "f, y %g\n", i, x(double) y); getchar();
  return x;
}


enum {
  SLOWDHT_NODISK = 0,     /* do not use disk */
  SLOWDHT_USEDISK = 1,    /* use disk is necessary */
  SLOWDHT_FORCEDISK = 2,  /* always use the disk to save memory */
  SLOWDHT_NOJJJ = 3       /* do not allocate the Jjj matrix */
};



#define slowdht_new(size, nu, xmax) \
        slowdht_newx(size, nu, xmax, SLOWDHT_USEDISK)

__inline slowdht *slowdht_newx(size_t size, xdouble nu, xdouble xmax,
    int usedisk)
{
  slowdht *dht;
  size_t i, tabsize, m, id, wb;
  xdouble x;

  if ( (dht = calloc(1, sizeof(*dht))) == NULL ) {
    fprintf(stderr, "no memory for dht\n");
    exit(1);
  }
  dht->size = size;
  dht->nu = nu;
  dht->inu = (int) (nu + .5);
  if ( FABS(dht->nu - dht->inu) > 1e-3 )
    dht->inu = -1;
  dht->xmax = xmax;

  if ((dht->j = calloc(dht->size + 2, sizeof(xdouble))) == NULL) {
    fprintf(stderr, "no memory for dht->j %d\n", (int) dht->size);
    exit(1);
  }

  /* compute the zeros of the Bessel functions
   * Note, the 0th zero of nu = 0 does not exist */
  dht->j[0] = 0;
  for ( i = 1; i <= size + 1; i++ ) {
    dht->j[i] = besselzeroJnu(nu, i);
  }
  dht->kmax = dht->j[dht->size+1] / xmax;

  printf("XXX\n");
  if ((dht->J2 = calloc(dht->size + 1, sizeof(xdouble))) == NULL) {
    fprintf(stderr, "no memory for dht->J2 %d\n", (int) dht->size);
    exit(1);
  }
  for ( i = 0; i <= size; i++ ) {
    x = besselJnu(nu + 1, dht->j[i]);
    dht->J2[i] = x*x;
  }

#if HAVEF128
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat"
#pragma GCC diagnostic ignored "-Wformat-extra-args"
#endif

  sprintf(dht->Jjjfn, "DHTnu%" XDBLPRNF "gsz%ld%s.dat",
      dht->nu, (long) dht->size, STRPREC);

#if HAVEF128
#pragma GCC diagnostic pop
#endif

  dht->Jjj = NULL;
  dht->Jjjfp = NULL;
  dht->Jjjarr = NULL;
  if ( usedisk == SLOWDHT_NOJJJ )
    return dht;

  tabsize = size * (size + 1) / 2;
  if (  usedisk != SLOWDHT_FORCEDISK
    && (dht->Jjj = calloc(tabsize, sizeof(xdouble))) != NULL ) {
    /* compute the table and save it in the memory */
    fprintf(stderr, "slowdht: using memory to save the transform coefficients (%gG)...\n",
        (double) size*(size+1)/2*sizeof(xdouble)/(1024.*1024*1024));
    id = 0;
    for ( m = 0; m < size; m++ ) {
      if ( !slowdht_quiet )
        fprintf(stderr, "computing the coefficients... %ld/%ld = %7.4lf%% \r",
            (long) id, (long) tabsize, id*100./tabsize);
      for ( i = m;  i < size; i++ ) {
        x = dht->j[m+1] * dht->j[i+1] / dht->j[size + 1];
        x = besselJnu(dht->nu, x);
        dht->Jjj[ id++ ] = x;
      }
    }
  } else if ( usedisk == SLOWDHT_NODISK ) {
    fprintf(stderr, "no memory for dht->Jjj %ld\n", (long)tabsize);
  } else {
    /* compute the coefficients and save it to the disk */
    if ( (dht->Jjjarr = calloc(size*slowdht_block, sizeof(xdouble)))
          == NULL ) {
      fprintf(stderr, "no memory for buffer array\n");
      exit(1);
    }

    /* open the Jjj file */
    if ( (dht->Jjjfp = fopen(dht->Jjjfn, "rb")) != NULL ) {
      /* check if the file is good */
      long fpos;
      fseek(dht->Jjjfp, 0, SEEK_END);
      fpos = ftell(dht->Jjjfp);
      if ( (size_t) fpos == size * size * sizeof(xdouble) ) {
        goto END_JFILE;
      } else {
        fprintf(stderr, "slowdht: %s bad size %ld vs %ld\n",
            dht->Jjjfn, fpos, (long) (size*size*sizeof(xdouble)));
      }
    }

    /* if the file does not exist, open it in update mode */
    dht->Jjjfp = fopen(dht->Jjjfn, "wb+");
    if ( dht->Jjjfp == NULL ) {
      fprintf(stderr, "cannot open %s\n", dht->Jjjfn);
      exit(1);
    }

    fprintf(stderr, "slowdht: to save the transform coefficients (%gG) to %s...\n",
        (double) size*size*sizeof(xdouble)/(1024.*1024*1024), dht->Jjjfn);

    for ( m = 0; m < size; m += slowdht_block ) {
      size_t k;

      id = 0;
      for ( k = 0; k < slowdht_block; k++ ) {
        if ( m+k >= size ) break;
        for ( i = 0; i < size; i++ ) {
          x = dht->j[m+k+1] * dht->j[i+1] / dht->j[size + 1];
          x = besselJnu(dht->nu, x);
          dht->Jjjarr[id++] = x;
        }
      }

      if ( (wb = fwrite(dht->Jjjarr, sizeof(xdouble), size*k, dht->Jjjfp))
           != size*k ) {
        fprintf(stderr, "\nfwrite: failed at m %ld, written %ld. \n",
            (long) m, (long) wb);
        exit(1);
      }
      if ( !slowdht_quiet )
        fprintf(stderr, "%s: computing the coefficients... %ld/%ld = %7.4lf%% \r",
            dht->Jjjfn, (long) m, (long) size, m*100./size);
    } /* end of the loop of m */

END_JFILE:
    ; /* open the Jjj file */
  }
  return dht;
}



__inline xdouble slowdht_x_sample(slowdht *dht, int n)
{
  if ( (size_t) n <= dht->size )
    return dht->xmax * dht->j[n+1] / dht->j[dht->size+1];
  return 0;
}



__inline xdouble slowdht_k_sample(slowdht *dht, int n)
{
  if ( (size_t) n <= dht->size )
    return dht->j[n+1] / dht->xmax;
  return 0;
}



/* get the pair index from 0 to n*(n + 1)/2 - 1 */
__inline int getpridx(int i, int j, int n)
{
  if ( i < 0 || i >= n || j < 0 || j >= n ) {
    fprintf(stderr, "bad index error i %d, j %d, n %d\n", i, j, n);
    exit(1);
  }
  if ( i > j ) { int i1 = i; i = j; j = i1; }
  return n*i - (i + 1)*i/2 + j;
  /* i = 0, j = 0, 1, ..., n-1 [n]
   * i = 1, j = 1, 2, ..., n-1 [2*n-1]
   * i = 2, j = 2, 3, ..., n-1 [3*n-3]
   * ...
   * i = i, ...,             [(i+1)*n - i*(i+1)/2]
   * id = i*n-i*(i-1)/2 + (j-i).
   * */
}


__inline int slowdht_apply(const slowdht *dht, xdouble *inp, xdouble *out)
{
  size_t m, i, size = dht->size, wb, k, id, rows;
  xdouble x, y;

  if ( dht->Jjj == NULL && dht->Jjjarr == NULL ) {
    fprintf(stderr, "cannot apply discrete Hankel transform\n");
    return -1;
  }
  if ( dht->Jjjfp != NULL ) rewind(dht->Jjjfp);

  for ( m = 0; m < size; m += slowdht_block ) {
    /* read several row of the transform coefficients */
    if ( (rows = slowdht_block) > size - m ) rows = size - m;

    if ( dht->Jjjarr ) {
      if ( dht->Jjjfp != NULL ) {
        if ( (wb = fread(dht->Jjjarr, sizeof(xdouble), size*rows, dht->Jjjfp))
             != size*rows ) {
          fprintf(stderr, "fread: failed at m %ld, read %ld. \n",
              (long) m, (long) wb);
          exit(1);
        }
      }
      if ( !slowdht_quiet )
        fprintf(stderr, "DHT %ld/%ld = %7.4f%%. \r",
            (long) m, (long) size, (double) m * 100 / size);
    }

    id = 0;
    for ( k = 0; k < rows; k++ ) {
      y = 0;
      for ( i = 0; i < size; i++ ) {
        if ( dht->Jjjarr != NULL ) { /* from disk buffer */
          x = dht->Jjjarr[id++];
        } else if ( dht->Jjj != NULL ) { /* from memory */
          x = dht->Jjj[ getpridx(m+k, i, size) ];
        } else { /* compute directly */
          x = dht->j[m+k+1] * dht->j[i+1] / dht->j[size+1];
          x = besselJnu(dht->nu, x);
        }
        y += inp[i] * x / dht->J2[i+1];
      }
      x = dht->xmax / dht->j[size+1];
      out[m+k] = 2*(x*x)*y;
    }
  }
  return 0;
}



__inline void slowdht_free(slowdht *dht)
{
  if ( dht->Jjjfp ) fclose(dht->Jjjfp);
  if ( dht->Jjjarr ) free(dht->Jjjarr);
  if ( dht->Jjj ) free(dht->Jjj);
  free(dht->j);
  free(dht->J2);
  free(dht);
}



#endif /* SLOWDHT_H__ */

