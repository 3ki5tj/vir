#ifndef GR_H__
#define GR_H__



/* loading correlation functions */



/* load the reference g(r) that comes from MC/MD
 * the input file should be the output of ljrdf_save()
 * in the LJ module of zcom.h, with the name "rdfxxx.dat" */
__inline double *loadgrraw(const char *fn, int *npt, double **ri)
{
  double *gr = NULL, x, y, z, xmin, dx, l;
  int version, i, j, rows, d, nfr, num;
  unsigned flags;
  char s[1024], *p;
  FILE *fp;

  *ri = NULL;

  xfopen(fp, fn, "r", return NULL);

  /* read the number of points */
  if ( fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "bad file %s\n", fn);
    goto EXIT;
  }
  if ( (i = sscanf(s + 1, " %d %x | %d %d %lf %lf | ",
        &version, &flags, &rows, npt, &xmin, &dx))
      != 6 ) {
    fprintf(stderr, "cannot load information from %s, %d\n%s", fn, i, s);
    goto EXIT;
  }
  if ( (p = strstr(s, "RDF")) != NULL ) {
    if (4 == sscanf(p + 3, "%d %d %d %lf", &nfr, &d, &num, &l)) {
      rho = num / POW(l, d);
      printf("set density to %g\n", (double) rho);
    }
  }

  /* load the array */
  xnew(gr, *npt);
  xnew(*ri, *npt);
  for ( i = 0; i < *npt; i++ )
    (*ri)[i] = (i + .5) * dx;
  for ( i = 0; i < *npt; i++ ) {
    if ( fgets(s, sizeof s, fp) == NULL
      || sscanf(s, "%lf %lf %lf", &x, &y, &z) != 3 ) {
      fprintf(stderr, "%s stopped at i %d\n", fn, i);
      break;
    }
    j = (int) (x/dx);
    if ( j < *npt ) gr[j] = z;
  }
  /* remove the last entry for safe */
  *npt -= 1;
EXIT:
  fclose(fp);
  return gr;
}



#ifndef INTERP
#define INTERP
/* given the array (xi, yi), evaluate the value at x */
__inline double interp(double x, double *xi, double *yi,
    int imin, int imax, int cutimin, int cutimax)
{
  int i;
  double gam;

  for ( i = imin; i < imax; i++ ) if ( xi[i] > x ) break;
  if ( i == imax ) { /* extrapolate */
    if ( cutimax ) return 0;
    i--;
  } else if ( i == imin ) { /* extrapolate */
    if ( cutimin ) return 0;
    i++;
  }
  /* linear interpolation */
  gam = (xi[i] - x) / (xi[i] - xi[i-1]);
  return gam * yi[i-1] + (1 - gam) * yi[i];
}
#endif



/* load g(r) */
__inline int loadgr(const char *fn, int npt, xdouble *gr,
    const xdouble *ri, xdouble sdr)
{
  double *ri0 = NULL, *gr0 = NULL, r, r0max;
  int i, n0 = 0;

  /* load the raw data */
  gr0 = loadgrraw(fn, &n0, &ri0);
  if ( gr0 == NULL ) return -1;
  r0max = ri0[n0-1];

  if (sdr > 0) {
    double s = 0, sg = 0;
    for ( i = 0; i < n0; i++ )
      if ( ri0[i] > r0max - sdr ) {
        s += 1;
        sg += gr0[i];
      }
    sg /= s;
    printf("scale g(r) by %g\n", sg);
    for ( i = 0; i < n0; i++ )
      gr0[i] /= sg;
  }

  for ( i = 0; i < npt; i++ ) {
    r = (double) ri[i];
    if ( r > r0max ) {
      gr[i] = 1;
    } else {
      gr[i] = interp(r, ri0, gr0, 0, n0, 1, 1);
    }
  }
  free(gr0);
  free(ri0);
  return 0;
}



/* load the reference w(r) = lny(r) that comes from MC/MD
 * the input file should be the output of mc1.c
 * in the LJ module of zcom.h, with the name pmfxxx.dat */
__inline double *loadwrraw(const char *fn, int *npt, double **ri)
{
  double *wr = NULL, x, y, y3, y4, y5, dx, tp, den;
  int i, j, d, num;
  char s[1024];
  FILE *fp;

  *ri = NULL;

  xfopen(fp, fn, "r", return NULL);

  /* read the number of points */
  if ( fgets(s, sizeof s, fp) == NULL || s[0] != '#') {
    fprintf(stderr, "bad file %s\n", fn);
    goto EXIT;
  }
  if ( (i = sscanf(s + 1, " %d %lf %d %d %lf %lf",
        npt, &dx, &d, &num, &tp, &den))
      != 6 ) {
    fprintf(stderr, "cannot load information from %s, %d\n%s", fn, i, s);
    goto EXIT;
  }
  if ( FABS(den - rho) > 1e-6 ) {
    fprintf(stderr, "Warning: rho %g vs. %g (%s)\n", (double) rho, den, fn);
  }

  rho = den;
  T = tp;
  beta = 1/T;
  fprintf(stderr, "set density %g, T %g\n", (double) rho, (double) T);

  /* load the array */
  xnew(wr, *npt);
  xnew(*ri, *npt);
  for ( i = 0; i < *npt; i++ )
    (*ri)[i] = (i + .5) * dx;
  for ( i = 0; i < *npt; i++ ) {
    if ( fgets(s, sizeof s, fp) == NULL
      || sscanf(s, "%lf %lf %lf %lf %lf", &x, &y, &y3, &y4, &y5) != 5 ) {
      fprintf(stderr, "%s stopped at i %d\n", fn, i);
      break;
    }
    j = (int) (x/dx);
    if ( j < *npt )
      wr[j] = y/tp;
  }
  /* remove the last entry for safe */
  *npt -= 1;
EXIT:
  fclose(fp);
  return wr;
}



/* load w(r) */
__inline int loadwr(const char *fn, int npt, xdouble *wr, const xdouble *ri)
{
  double *ri0 = NULL, *wr0 = NULL, r, r0max;
  int i, n0 = 0;

  /* load the raw data */
  wr0 = loadwrraw(fn, &n0, &ri0);
  //printf("loaded wr %s n0 %d %p\n", fn, n0, wr0); getchar();
  if ( wr0 == NULL ) return -1;
  r0max = ri0[n0-1];

  for ( i = 0; i < npt; i++ ) {
    r = (double) ri[i];
    if ( r > r0max ) {
      wr[i] = 0;
    } else {
      wr[i] = interp(r, ri0, wr0, 0, n0, 0, 1);
    }
  }
  free(wr0);
  free(ri0);
  return 0;
}



#endif /* defined(GR_H__) */
