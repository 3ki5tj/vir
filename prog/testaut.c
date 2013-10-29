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



static void testspeed(int n, int nsamp, int nedmax)
{
  int t, ned, eql = 1, nequil = 1000, isamp, good = 0, tot = 0;
  dg_t *g, *ng;
  clock_t t0;
  double sum = 0, tsum = 0, rnp = 0.1;

  printf("speed test for n %d, nsamp %d, nedmax %d\n",
      n, nsamp, nedmax);
  g = dg_open(n);
  ng = dg_open(n);
  dg_full(g);
  ned = dg_nedges(g);

  for (t = 1; ; t++) {
    dg_rndswitchedge(g, &ned, nedmax, rnp);
    if (eql && t >= nequil) {
      t = 0;
      eql = 0;
      continue;
    }
    if (t % 10 != 0) continue;
    adjustrnp(ned, nedmax, t, 1000000, &good, &tot, &rnp);
    if ( ned > nedmax ) continue;
    if ( dg_cliquesep(g) ) continue;

    t0 = clock();
    dg_canlabel(ng, g);
    tsum += clock() - t0;
    sum += dg_nedges(ng);
    if (++isamp >= nsamp) break;
  }
  tsum /= CLOCKS_PER_SEC;
  printf("canonical label, n %d, samples %d, steps %d, nedmax %d; "
      "time used: %gs/%d = %gmcs\n",
      n, nsamp, t, nedmax, tsum, nsamp, tsum / nsamp * 1e6);
  dg_close(g);
  dg_close(ng);
}



int main(int argc, char **argv)
{
  int n = 12, nsamp = 1000000, nedmax = 1000000;
  testfoo(6);
  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) nsamp = atoi(argv[2]);
  if (argc >= 4) nedmax = atoi(argv[3]);
  testspeed(n, nsamp, nedmax);
  return 0;
}

