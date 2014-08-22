#include "dgcryr.h"
#include "testutil.h"



static void test4(void)
{
  dg_t *g;
  int i, j, n = 4;

  g = dg_open(n);
  dg_link(g, 0, 1);
  dg_link(g, 0, 2);
  dg_link(g, 0, 3);
  dg_link(g, 1, 2);
  dg_link(g, 1, 3);
  //dg_link(g, 2, 3);
  dg_print(g);
  for ( i = 0; i < n; i++ )
    for ( j = i + 1; j < n; j++ ) {
      double fb1 = dgsc_yriter(g, i, j, 0);
      double fb1d = dgsc_yriter(g, i, j, 1);
      double fb1w = dgsc_yriter(g, i, j, -1);
      double fb2 = dgrjw_yrfb(g, i, j);
      printf("i %d, j %d, %g, %g | %g, %g\n", i, j, fb1, fb2, fb1d, fb1w);
    }
  dg_close(g);
}



static void test5(void)
{
  dg_t *g;
  int i, j, n = 5;

  g = dg_open(n);
  dg_link(g, 0, 1);
  dg_link(g, 1, 2);
  dg_link(g, 2, 3);
  dg_link(g, 3, 4);
  dg_link(g, 4, 0);
  dg_link(g, 2, 4);
  dg_print(g);
  for ( i = 0; i < n; i++ )
    for ( j = i + 1; j < n; j++ ) {
      double fb1 = dgsc_yriter(g, i, j, 0);
      double fb1d = dgsc_yriter(g, i, j, 1);
      double fb1w = dgsc_yriter(g, i, j, -1);
      double fb2 = dgrjw_yrfb(g, i, j);
      printf("i %d, j %d, %g, %g | %g, %g\n", i, j, fb1, fb2, fb1d, fb1w);
    }
  dg_close(g);
}



/* check if the Wheatley's method yields the same result with the direct method */
static void testconsist(int n, int nsamp, int nedmin)
{
  int i, j, t, ned, eql = 0, nequil = 1000, isamp = 0;
  dg_t *g;
  double fb1, fb2;
  double rnm = 0.5; /* rate for moves that increase edges */

  g = dg_open(n);
  dg_full(g);
  ned = dg_nedges(g);

  for (t = 1; ; t++) {
    dg_rndswitchedge0(g, &ned, nedmin, rnm, 1000000, 1);
    if (eql && t >= nequil) {
      t = 0;
      eql = 0; /* stop equilibration */
      continue;
    }

    //adjustrnm(ned, nedmin, t, 1000000, &good, &tot; &rnm);
    if (t % 100 != 0) continue;

    for ( i = 0; i < n; i++ ) {
      for ( j = i + 1; j < n; j++ ) {
        fb1 = dgsc_yriter(g, i, j, 0);
        fb2 = dgrjw_yrfb(g, i, j);
        if ( fabs(fb1 - fb2) > 1e-4 ) {
          dg_print(g);
          fprintf(stderr, "i %d, j %d %g, %g\n", i, j, fb1, fb2);
          getchar();
        }
      }
    }
    if (++isamp >= nsamp) break;
  }
  dg_close(g);
  printf("consistency check for n = %d, finished for %d samples\n", n, nsamp);
}



int main(void)
{
  test4();
  test5();
  testconsist(6, 1000, 0);
  return 0;
}
