/* merge data from different number of edges */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <gmp.h>


int order = 12;
char *lspattern = "gaussfn%de*E*ldblmpf.dat";



static char **collectlist(int order, int *nfns)
{
  const char *cmdout = "output.lst";
  char cmd0[1024], cmd[1024], **fns = NULL, *p;
  FILE *fp;
  int i;

  sprintf(cmd0, lspattern, order);
  sprintf(cmd, "/bin/ls %s > %s", cmd0, cmdout);
  system(cmd);
  if ( (fp = fopen(cmdout, "r")) == NULL) return NULL;
  /* count the number of lines */
  for ( *nfns = 0; !feof(fp); )
    if (fgetc(fp) == '\n') (*nfns)++;
  printf("%d files\n", *nfns);
  if ((fns = calloc(*nfns, sizeof(*fns))) == NULL)
    exit(1);
  if ((fns[0] = calloc(*nfns*FILENAME_MAX, sizeof(char))) == NULL)
    exit(1);

  rewind(fp);
  for ( i = 0; i < *nfns; i++ ) {
    if ( i > 0 ) fns[i] = fns[0] + i*FILENAME_MAX;
    fgets(fns[i], FILENAME_MAX, fp);
    /* trim the trailing spaces */
    for ( p = fns[i] + strlen(fns[i]) - 1; p >= fns[i] && isspace(*p); )
      *p-- = '\0';
    printf("%4d: %s\n", i + 1, fns[i]);
  }
  fclose(fp);
  remove(cmdout);
  return fns;
}



static int merge(int order, int nfns, char **fns)
{
  int i, dim, d, d1;
  FILE *fp;
  char buf[512], fn[80];
  mpf_t *arr;

  if ( nfns == 0 ) return -1;

  /* get the maximal dimension */
  if ((fp = fopen(fns[0], "r")) == NULL) exit(1);
  fgets(buf, sizeof buf, fp);
  if ( buf[0] != '#' || sscanf(buf + 1, " D %d", &dim) != 1 ) {
    fprintf(stderr, "cannot get dimensions from %s\n%s", fns[0], buf);
    exit(1);
  }
  fclose(fp);

  mpf_set_default_prec(256);

  /* initialize the array */
  if ( (arr = calloc(dim + 1, sizeof(*arr))) == NULL )
    exit(1);
  for ( d = 0; d <= dim; d++ ) {
    mpf_init(arr[d]);
    mpf_set_si(arr[d], 0);
  }

  /* add up all parts of the virial coefficients */
  for ( i = 0; i < nfns; i++ ) {
    if ( (fp = fopen(fns[i], "r")) == NULL ) {
      fprintf(stderr, "cannot open %s\n", fns[i]);
      exit(1);
    }
    fgets(buf, sizeof buf, fp);
    for (d = 1; d <= dim; d++) {
      if (fgets(buf, sizeof buf, fp) == NULL) {
        fprintf(stderr, "%s corrupted on line %d\n", fns[i], d);
        break;
      }
      gmp_sscanf(buf, "%d%Ff", &d1, arr[0]);
      mpf_add(arr[d], arr[d], arr[0]);
    }
    fclose(fp);
  }

  /* write the output */
  sprintf(fn, "gaussfn%dldblmpf.dat", order);
  if ( (fp = fopen(fn, "w")) == NULL ) {
    fprintf(stderr, "cannot write %s\n", fn);
    exit(1);
  }
  fprintf(fp, "# D %d B%d  ldbl mpf 256\n", dim, order);
  for ( d = 1; d <= dim; d++ ) {
    gmp_fprintf(fp, "%4d %+.76Fe\n", d, arr[d]);
  }
  fclose(fp);
  fprintf(stderr, "saved results to %s\n", fn);

  for ( d = 0; d <= dim; d++ )
    mpf_clear(arr[d]);
  free(arr);
  return 0;
}



int main(int argc, char **argv)
{
  char **fns;
  int nfns;

  if (argc > 1) order = atoi(argv[1]);
  printf("merging all order %d files\n", order);
  fns = collectlist(order, &nfns);
  merge(order, nfns, fns);
  free(fns[0]);
  free(fns);
  return 0;
}

