#include "dgaut.h"
#include "dgcsep.h"
#include "testutil.h"


static void testfoo(int n)
{
  dg_t *g1, *g2, *g1a, *g2a;
  g1 = dg_open(n);
  dg_full(g1);
  dg_unlink(g1, 2, 3); dg_unlink(g1, 3, 1); dg_unlink(g1, 1, 0);
  g2 = dg_open(n);
  dg_full(g2);
  dg_unlink(g2, 0, 3); dg_unlink(g2, 0, 2); dg_unlink(g2, 2, 1);
  g1a = dg_open(n);
  g2a = dg_open(n);

  printf("before:\n");
  dg_print(g1);
  dg_print(g2);
  dg_canlabel(g1a, g1);
  dg_canlabel(g2a, g2);
  printf("after:\n");
  dg_print(g1a);
  dg_print(g2a);
  dg_close(g1); dg_close(g1a);
  dg_close(g2); dg_close(g2a);
}



static void testequipart(void)
{
  dg_t *g;
  int pairs[][2] = {
    {0, 1}, {1, 2}, {0, 3}, {3, 6},
    {3, 4}, {4, 5}, {1, 4}, {4, 7},
    {6, 7}, {7, 8}, {2, 5}, {5, 8},
    {-1, -1}};
  dgpart_t part;

  g = dg_open(9);
  dg_linkpairs(g, pairs);
  printf("3x3 graph\n");
  dg_print(g);
  dgpart_unit(&part, g->n);
  dg_equipart(&part, g);
  dgpart_print(&part);

  printf("recursive\n");
  dgpart_unit(&part, g->n);
  while (part.nc < g->n) {
    int ip;
    dgword_t vs, b;

    dg_equipart(&part, g);
    if (part.nc == g->n) break;
    /* find the first unfixed cell */
    for (ip = 0; ip < part.nc; ip++)
      if (part.cnt[ip] > 1)
        break;
    /* artificially break the first cell */
    vs = part.cs[ip];
    b = vs & (-vs);
    part.cs[part.nc] = vs ^ b;
    part.cnt[part.nc] = part.cnt[ip] - 1;
    part.cs[ip] = b;
    part.cnt[ip] = 1;
    part.nc += 1;
    dgpart_print(&part);
  }
  dgpart_print(&part);
  dg_close(g);
}



static void testspeed(int n, int nsamp, int nedmax)
{
  int t, ned, eql = 1, nequil = 1000, good = 0, tot = 0, level;
  dg_t *g, *ng;
  clock_t t0;
  double sum = 0, tsum[] = {0, 0, 0, 0, 0}, rnp = 0.1;
  const char *names[] = {
    "copy                     ",
    "degree sequence          ",
    "equitable partition      ",
    "deep equitable partition ",
    "canonical label          "};
#define LISTSIZE 1000
  dg_t *gls[LISTSIZE] = {NULL};
  int lscnt = 0, kk;

  printf("speed test for n %d, nsamp %d, nedmax %d\n",
      n, nsamp, nedmax);
  g = dg_open(n);
  ng = dg_open(n);
  for (t = 0; t < LISTSIZE; t++)
    gls[t] = dg_open(n);
  dg_full(g);
  ned = dg_nedges(g);

  for (t = 1; ; t++) {
    dg_rndswitchedge(g, &ned, nedmax, rnp);
    if (eql && t >= nequil) {
      t = 0;
      eql = 0;
      lscnt = 0;
      continue;
    }
    if ( ned <= nedmax ) {
      die_if (lscnt >= LISTSIZE, "lscnt %d, t %d\n", lscnt, t);
      dg_copy(gls[lscnt], g);
      lscnt++;
    }
    adjustrnp(ned, nedmax, t, 1000000, &good, &tot, &rnp);
    if (t % LISTSIZE != 0) continue;

    for (level = 0; level <= 4; level++) {
      t0 = clock();
      for (kk = 0; kk < lscnt; kk++)
        dg_repiso(ng, gls[kk], level);
      tsum[level] += clock() - t0;
    }
    sum += lscnt;
    lscnt = 0;
    if (sum >= nsamp) break;
  }
  printf("canonical label, n %d, samples %d, steps %d, nedmax %d\n",
      n, nsamp, t, nedmax);
  for (level = 0; level <= 4; level++) {
    tsum[level] /= CLOCKS_PER_SEC;
    printf("level %d, %s  time used: %gs/%g = %gmcs\n",
        level, names[level], tsum[level], sum, tsum[level] / sum * 1e6);
  }
  dg_close(g);
  dg_close(ng);
  for (t = 0; t < LISTSIZE; t++)
    dg_close(gls[t]);
}



int main(int argc, char **argv)
{
  int n = 12, nsamp = 10000000, nedmax = 1000000;

#ifdef N
  n = N;
#else
  testfoo(6);
  testequipart();
  if (argc >= 2) n = atoi(argv[1]);
#endif

  if (argc >= 3) nsamp = atoi(argv[2]);
  if (argc >= 4) nedmax = atoi(argv[3]);

  /* Timing on T60 with the default setting
   *                          general n   -DN=12
   * copy                       0.02mcs   0.01mcs
   * degree sequence            0.8mcs    0.26mcs
   * equitable partition        1.7mcs    1.0mcs
   * deep equitable partition   1.6mcs    1.1mcs
   * canonical label            2.4mcs    2.4mcs
   * */
  testspeed(n, nsamp, nedmax);

  return 0;
}

