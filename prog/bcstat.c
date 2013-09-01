/* compute the probability of obtaining biconnected diagrams
 * by the standard Monte Carlo integration */
#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"



/* compute the probability of obtaining a biconnected diagram */
static void bcstat(int nmax, double nsteps)
{
  int i, j, n;
  code_t vs;
  rvn_t *v, u;
  double t, *cnt;
  dg_t *g;

  g = dg_open(nmax);
  xnew(v, nmax + 1);
  xnew(cnt, nmax + 1);
  for (n = 0; n <= nmax; n++) cnt[n] = 0;

  rvn_zero(v[0]);
  for (t = 1; t <= nsteps; t += 1) {
    /* linear chain */
    for (n = 1; n < nmax; n++)
      rvn_add(v[n], v[n - 1], rvn_rndball0(u));

    /* construct the diagram */
    dg_empty(g);
    for (i = 0; i < nmax; i++)
      for (j = i + 1; j < nmax; j++)
        if (rvn_dist2(v[i], v[j]) < 1)
          dg_link(g, i, j);

    /* test biconnectivity of chain subdiagrams */
    for (n = 3; n <= nmax; n++) { /* over n-vertex graphs */
      for (i = 0; i <= nmax - n; i++) { /* offset */
        /* construct the vertex set */
        vs = (((code_t) 1u << n) - 1) << i;
        if (dg_biconnectedvs(g, vs))
          cnt[n] += 1.;
      }
    }
  }

  /* compute the ratio of biconnected diagrams */
  printf("%2d %-14d", D, 1);
  for (n = 3; n <= nmax; n++) {
    cnt[n] /= (nmax - n + 1) * nsteps;
    printf(" %.12f", cnt[n]);
  }
  printf("\n");
  free(v);
  free(cnt);
  dg_close(g);
}



int main(int argc, char *argv[])
{
  int nmax = 30;
  double nsteps = 1e7;

  if (argc > 1) nmax = atoi(argv[1]);
  if (argc > 2) nsteps = atof(argv[2]);
  bcstat(nmax, nsteps);
  mtsave(NULL);
  return 0;
}

