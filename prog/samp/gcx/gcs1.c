/* radius-based cheap reaction coordinate */
#define ZCOM_PICK
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"

int nmin = 3;
int nmax = DG_NMAX - 1;



/* step in a pure state */
static int bcstep(dg_t *g, dg_t *ng, rvn_t *x, rvn_t xi, real amp)
{
  int i, j, n = g->n;

  i = (int) (n * rnd0());
  rvn_rnddisp(xi, x[i], amp);
  ng->n = n;
  dg_copy(ng, g);
  for (j = 0; j < n; j++) {
    if (i == j) continue;
    if (rvn_dist2(xi, x[j]) < 1)
      dg_link(ng, i, j);
    else
      dg_unlink(ng, i, j);
  }
  if ( dg_biconnected(ng) ) {
    dg_copy(g, ng);
    rvn_copy(x[i], xi);
    return 1;
  }
  return 0;
}



/* step in a restricted state r(0, n - 1) < r0 */
static int bcrstep(dg_t *g, dg_t *ng, rvn_t *x, rvn_t xi,
    real amp, real r0)
{
  int i, j, n = g->n;

  i = (int) (n * rnd0());
  rvn_rnddisp(xi, x[i], amp);
  if (i == 0) { /* check if the (0, n - 1) bond is broken */
    if (rvn_dist2(xi, x[n - 1]) > r0 * r0)
      return 0;
  } else if (i == n - 1) {
    if (rvn_dist2(xi, x[0]) > r0 * r0)
      return 0;
  }
  ng->n = n;
  dg_copy(ng, g);
  for (j = 0; j < n; j++) {
    if (i == j) continue;
    if (rvn_dist2(xi, x[j]) < 1)
      dg_link(ng, i, j);
    else
      dg_unlink(ng, i, j);
  }
  if ( dg_biconnected(ng) ) {
    dg_copy(g, ng);
    rvn_copy(x[i], xi);
    return 1;
  }
  return 0;
}




/* restricted state n + 1 to n, remove the last particle */
static int nmove_remove(dg_t *g) 
{
  if ( dg_biconnectedvs(g, (1u << (g->n - 1)) - 1) ) {
    return 1;
  }
  return 0;
}



/* the n + 1 simulation */
static double foo_n1(int nsteps, int n, real rc)
{
  dg_t *g, *ng;
  rvn_t x[DG_NMAX] = {{0}}, xi;
  int t, nacc = 0, ntot = 0, cacc = 0, ctot = 0;
  real amp = 1.5/D;

  g = dg_open(DG_NMAX - 1);
  ng = dg_open(DG_NMAX - 1);
  g->n = n + 1;
  dg_full(g);

  for (t = 1; t <= nsteps; t++) {
    if (rnd0() < 0.1) {
      ntot++;
      nacc += nmove_remove(g);
    } else {
      ctot++;
      cacc += bcrstep(g, ng, x, xi, amp, rc);
      die_if (rvn_dist(x[0], x[g->n - 1]) > rc, "bad t %d\n", t);
    }
  }
  printf("%d --> %d, nacc %g, cacc %g\n",
      n + 1, n, 1.*nacc/ntot, 1.*cacc/ctot);
  dg_close(g);
  dg_close(ng);
  return 1.*nacc/ntot;
}



/* add a particle */
static int nmove_add(dg_t *g, dg_t *ng, rvn_t *x, real *xi,
    int conn[], real rc)
{
  int deg, j, n = g->n;

  rvn_rndball(xi, rc);
  rvn_inc(xi, x[0]);
  for (deg = 0, j = 0; j < n; j++) {
    if (rvn_dist2(xi, x[j]) < 1) {
      conn[j] = 1;
      deg++;
    } else {
      conn[j] = 0;
    }
  }
  if (deg >= 2) { /* accept */
    ng->n = n;
    dg_copy(ng, g);
    ng->n = n + 1;
    for (j = 0; j < n; j++) {
      if (conn[j]) dg_link(ng, j, n);
      else dg_unlink(ng, j, n);
    }
    return 1;
  }
  return 0;
}




/* n to n + 1*/
static double foo_n(int nsteps, int n, real rc)
{
  dg_t *g, *ng;
  rvn_t x[DG_NMAX] = {{0}}, xi;
  int t, nacc = 0, ntot = 0, cacc = 0, ctot = 0;
  real amp = 1.5/D;
  int conn[DG_NMAX];

  g = dg_open(DG_NMAX - 1);
  ng = dg_open(DG_NMAX - 1);
  g->n = n;
  dg_full(g);

  for (t = 1; t <= nsteps; t++) {
    if (rnd0() < 0.1) {
      ntot++;
      nacc += nmove_add(g, ng, x, xi, conn, rc);
    } else {
      ctot++;
      cacc += bcstep(g, ng, x, xi, amp);
    }
  }
  printf("%d --> %d, nacc %g, cacc %g\n",
      n, n + 1, 1.*nacc/ntot, 1.*cacc/ctot);
  dg_close(g);
  dg_close(ng);
  return 1.*nacc/ntot;
}



/* scale the distance between x[n-1] and x[0] */
static int nmove_scale(dg_t *g, dg_t *ng, rvn_t *x, real *xi, real sig)
{
  int j, n1 = g->n - 1;

  rvn_diff(xi, x[n1], x[0]);
  rvn_smul(xi, sig);
  rvn_inc(xi, x[0]);
  ng->n = g->n;
  dg_copy(ng, g);
  for (j = 0; j < n1; j++) {
    if ( rvn_dist2(x[j], xi) < 1 )
      dg_link(ng, j, n1);
    else
      dg_unlink(ng, j, n1);
  }
  if ( dg_biconnected(ng) ) {
    return 1;
  }
  return 0;
}



/* restricted simulation */
static double foo_restricted(int nsteps, int n, real rc, real sig)
{
  dg_t *g, *ng;
  rvn_t x[DG_NMAX] = {{0}}, xi;
  int t, nacc = 0, ntot = 0, cacc = 0, ctot = 0;
  real amp = 1.5/D;

  g = dg_open(DG_NMAX - 1);
  ng = dg_open(DG_NMAX - 1);
  g->n = n + 1;
  dg_full(g);

  for (t = 1; t <= nsteps; t++) {
    if (rnd0() < 0.1) {
      ntot++;
      nacc += nmove_scale(g, ng, x, xi, sig);
    } else {
      ctot++;
      cacc += bcrstep(g, ng, x, xi, amp, rc);
      die_if (rvn_dist(x[0], x[g->n - 1]) > rc, "bad t %d\n", t);
    }
  }
  printf("%d restricted rc %g --> pure, nacc %.7f, cacc %.7f\n",
      n, rc, 1.*nacc/ntot, 1.*cacc/ctot);
  dg_close(g);
  dg_close(ng);
  return 1.*nacc/ntot;
}



/* pure-state simulation */
static double foo_pure(int nsteps, int n, real rc, real sig)
{
  dg_t *g, *ng;
  rvn_t x[DG_NMAX] = {{0}}, xi;
  int t, nacc = 0, ntot = 0, cacc = 0, ctot = 0;
  real amp = 1.5/D, rcs2 = rc * rc * sig * sig;
  int conn[DG_NMAX];

  g = dg_open(DG_NMAX - 1);
  ng = dg_open(DG_NMAX - 1);
  g->n = n;
  dg_full(g);

  for (t = 1; t <= nsteps; t++) {
    if (rnd0() < 0.1) {
      ntot++;
      if (rvn_dist2(x[0], x[g->n - 1]) < rcs2 ) {
        nacc += nmove_scale(g, ng, x, xi, 1./sig);
      }
    } else {
      ctot++;
      cacc += bcstep(g, ng, x, xi, amp);
    }
  }
  printf("%d pure --> restricted rc %g, nacc %.7f, cacc %.7f\n",
      n, rc, 1.*nacc/ntot, 1.*cacc/ctot);
  dg_close(g);
  dg_close(ng);
  return 1.*nacc/ntot;
}



int main(int argc, char **argv)
{
  real rc = 1, sig = 1;
  double r1, r2, r3, r4;
  int n = 31, nsteps = 1000000;

  if (argc > 1) n = atoi(argv[1]);
  if (argc > 2) nsteps = atoi(argv[2]);
  if (argc > 3) rc = (real) atof(argv[3]);

  //rc = 2 * pow(n, 1./D);
  printf("n %d, rc %g\n", n, rc);
  r2 = foo_n1(nsteps, n - 1, rc);
  r1 = foo_n(nsteps, n - 1, rc);

  if (argc > 4) sig = (real) atof(argv[4]);
  else sig = (real) pow(0.66 * n * (r2/r1), 1./D) / rc;
  printf("sig %g\n", sig);
  r3 = foo_restricted(nsteps, n, rc, sig);
  r4 = foo_pure(nsteps, n, rc, sig);
  printf("r1 %.7f, r2 %.7f, r12 %.7f\n"
      "r3 %.7f, r4 %.7f, r34 %.7f\nr %.7f\n",
      r1, r2, r1/r2, r3, r4, r3/r4, r1/r2*r3/r4*pow(rc*sig, D));
  mtsave(NULL);
  return 0;
}
