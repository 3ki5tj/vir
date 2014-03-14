#define ZCOM_PICK
#define ZCOM_RNG
#define ZCOM_ARGOPT
#define ZCOM_MDS
#include "zcom.h"



char *fninp = NULL;
int outdim = 3;
char *fnout = NULL;
int nsteps = 2000;
float amp = 1.0f;
float dt = 0.01f;



static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  argopt_add(ao, NULL, "!%s", &fninp, "input position file");
  argopt_add(ao, "-o", NULL, &fnout, "out position file");
  argopt_add(ao, "-d", "%d", &outdim, "output dimension");
  argopt_add(ao, "-1", "%d", &nsteps, "number of steps");
  argopt_add(ao, "-a", "%f", &amp, "move amplitude");
  argopt_add(ao, "-t", "%f", &dt, "amplitude along the force");
  argopt_parse(ao, argc, argv);
  if ( fninp == NULL ) argopt_help(ao);
  argopt_close(ao);
  if ( fnout == NULL ) {
    char *p;
    fnout = ssnew(256);
    sscpy(fnout, "mds");
    sscat(fnout, fninp);
    if ((p = strchr(fnout, '.')) == NULL)
      p = fnout + strlen(fnout);
    sprintf(p, "d%d.pos", outdim);
  }
}



static real *load(const char *fn, int *dim, int *n)
{
  FILE *fp;
  char s[4096], *p;
  int ln, i, j, next, frame = -1;
  real *x;

  xfopen(fp, fn, "r", return NULL);
  if ( fgets(s, sizeof s, fp) == NULL
    || (s[0] != '#')
    || sscanf(s + 1, "%d %d", dim, n) != 2 ) {
    fprintf(stderr, "%s: first line corrupted %s", fn, s);
    fclose(fp);
    return NULL;
  }
  xnew(x, (*n) * (*dim));
  for ( ln = 0; fgets(s, sizeof s, fp); ln++ ) {
    if ( (i = ln % (*n)) == 0 ) {
      if ( feof(fp) ) break;
    }
    for ( p = s, j = 0; j < *dim; j++, p += next ) {
      if ( sscanf(p, "%lf%n", &x[i*(*dim) + j], &next) != 1 ) {
        fprintf(stderr, "corrupted at ln %d, i %d, j %d\n%s", ln, i, j, s);
        goto END;
      }
    }
    sscanf(p, "%d", &frame);
  }
  printf("file dim %d, n %d, frame %d\n", *dim, *n, frame);
END:
  fclose(fp);
  return x;
}



/* compute the distance matrix */
static real *getdismat(real *x, int dim, int n)
{
  real *dismat, dx, dis;
  int i, j, k;

  xnew(dismat, n * n);
  for ( i = 0; i < n; i++ ) {
    dismat[i*n + i] = 0;
    for ( j = i + 1; j < n; j++ ) {
      for ( dis = 0, k = 0; k < dim; k++ ) {
        dx = x[i*dim + k] - x[j*dim + k];
        dis += dx * dx;
      }
      dis = (real) sqrt(dis);
      dismat[i*n + j] = dismat[j*n + j] = dis;
    }
  }
  return dismat;
}



static real domds(real *dismat, real *nx, int n, int ndim)
{
  real *y, *ym, *f;
  real ene, enemin;
  int i, j, k;

  xnew(y, n * ndim);
  xnew(f, n * ndim);
  xnew(ym, n * ndim);
  for ( i = 0; i < n * ndim; i++ )
    y[i] = nx[i];
  for ( enemin = 1e9, k = 0; k < nsteps; k++ ) {
    /* perturb the current configuration */
    mds_force(y, f, dismat, n, ndim);
    if (k > 0) {
      for ( i = 0; i < n; i++ )
        for ( j = 0; j < ndim; j++ )
          y[i*ndim + j] += (rnd0() - .5) * amp + f[i*ndim + j] * dt;
    }
    for ( i = 0; i < n * ndim; i++ )
      ym[i] = y[i];
    ene = mds_min0(ym, dismat, n, ndim, 2e-5);
    if ( ene < enemin ) {
      enemin = ene;
      for ( i = 0; i < n * ndim; i++ ) nx[i] = ym[i];
      printf("round %d: ene %g\n", k, enemin);
    }
    if (enemin < 1e-10) break;
  }
  free(y);
  free(ym);
  free(f);
  return enemin;
}



static int save(const char *fn, real *x, int n, int dim,
    int dim0, real *dismat, real *dismat0, real ene)
{
  int i, j;
  real *diff, dx, diffmax = 0, diffmin = 1e10;
  FILE *fp;

  xfopen(fp, fn, "w", return -1);

  xnew(diff, n);
  for ( i = 0; i < n; i++ ) diff[i] = 0;
  for ( i = 0; i < n - 1; i++ )
    for ( j = i + 1; j < n; j++ ) {
      dx = fabs(dismat[i*n + j] - dismat0[i*n + j]);
      diff[i] += dx;
      diff[j] += dx;
    }

  for ( i = 0; i < n; i++ )
    if ( diff[i] > diffmax )
      diffmax = diff[i];
    else if ( diff[i] < diffmin )
      diffmin = diff[i];
  fprintf(fp, "# %d %d %d %g %g %g\n", dim, n, dim0, ene, diffmax, diffmin);

  /* print out the coordinates */
  for ( i = 0; i < n; i++ ) {
    for ( j = 0; j < dim; j++ )
      fprintf(fp, "%.6f ", x[i*dim + j]);
    fprintf(fp, "%g\n", diff[i]);
  }
  /* print out distance matrix */
  for ( i = 0; i < n - 1; i++ )
    for ( j = i + 1; j < n; j++ )
      fprintf(fp, "# %d %d %g %g\n",
          i, j, dismat[i*n + j], dismat0[i*n + j]);
  fclose(fp);
  free(diff);
  return -1;
}



int main(int argc, char **argv)
{
  int dim, n, i, j;
  real *x, *dismat, *nx, *ndismat;
  real ene;

  doargs(argc, argv);
  x = load(fninp, &dim, &n);
  dismat = getdismat(x, dim, n);
  xnew(nx, n * outdim);
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < outdim; j++ )
      nx[i*outdim + j] = ( j < dim ) ? x[i*dim + j] : 0;
  ene = domds(dismat, nx, n, outdim);
  ndismat = getdismat(nx, outdim, n);
  save(fnout, nx, n, outdim, dim, ndismat, dismat, ene);
  printf("minimal energy %g, saved to %s\n", ene, fnout);
  free(x);
  free(nx);
  free(dismat);
  return 0;
}

