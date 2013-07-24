/* compute the ratio of two successive virial coefficients
 * importance sampler coupled to several uniform samplers
 * uniform samplers are distinguished by the number of clique
 * separators in the clique-separator decomposition
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 mcratn.c -lm */



#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"
#include "dgrjw.h"
#include "mcutil.h"


int n = 3; /* order */
double nequil = 100000; /* number of equilibration steps */
double nsteps = 10000000, nsteps2 = -1;
int usegrand = 0; /* normally distributed */
int nstfb = 0; /* interval of evaluting the weight */
int nstnmv = 100; /* frequency of the n move */
int nstrep = 100000000; /* interval of reporting */
int lookup = -1; /* if to use the lookup table */
int cachesize = 10; /* number of recently visited diagrams */

int nuniform; /* number of uniform samplers */
int nsystem; /* number of systems */
real mcamp[DG_NMAX];
double Z[DG_NMAX] = {1}; /* inverse weights (etimated partition function) */



/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;
  int i;
  double amp0 = 1.5, amp1 = 0.9, amp2 = 0.9;
  double Z1 = -1, Z2 = -1;

  ao = argopt_open(0);
  argopt_add(ao, "-n", "%d", &n, "order n");
  argopt_add(ao, "-0", "%lf", &nequil, "number of equilibration steps");
  argopt_add(ao, "-1", "%lf", &nsteps, "number of simulation steps");
  argopt_add(ao, "-2", "%lf", &nsteps2, "number of simulation steps, for n - 1");
  argopt_add(ao, "-c", "%d", &nuniform, "number of uniform sampler");
  argopt_add(ao, "-a", "%lf", &amp0, "MC amplitude for system unrestricted biconnected diagrams");
  argopt_add(ao, "-A", "%lf", &amp1, "MC amplitude for nonzero biconnected diagrams");
  argopt_add(ao, "--A2", "%lf", &amp2, "MC amplitude for the importance sampling");
  argopt_add(ao, "-Z", "%lf", &Z1, "relative partition function of system 1");
  argopt_add(ao, "--Z2", "%lf", &Z2, "relative partition function of system 2");
  argopt_add(ao, "-L", "%d", &lookup, "lookup table (-1: automatic)");
  argopt_add(ao, "-w", "%d", &nstfb, "interval of evaluating the weight");
  argopt_add(ao, "-y", "%d", &nstnmv, "interval of the n-move");
  argopt_add(ao, "-q", "%d", &nstrep, "interval of reporting");
  argopt_add(ao, "-G", "%b", &usegrand, "normally-distributed displacement");
  argopt_add(ao, "-C", "%d", &cachesize, "cache size of recently visited diagrams");
  argopt_parse(ao, argc, argv);

  /* decide whether to use lookup table or not */
  if (lookup < 0)
    lookup = (n <= DGMAP_NMAX);

  /* interval of computing fb */
  if (nstfb <= 0) {
    if (lookup)
      nstfb = 1;
    else /* allow more sampling for larger n */
      for (nstfb = 10, i = 8; i < n; i++) nstfb *= 3;
  }

  /* number of sampling systems
   * the maximal number of clique separators is n - 3, hence n - 2
   * different values */
  if (nuniform <= 0)
    nuniform = (n <= 4) ? 2 : (n - 2);
  nsystem = nuniform + 1;
  die_if (nuniform == 1, "nuniform = %d not supported\n", nuniform);

  /* Monte Carlo move amplitude */
  mcamp[0] = (real) amp0;
  if (nuniform > 1) {
    for (i = 1; i < nuniform - 1; i++)
      mcamp[i] = (real) (amp0 * pow(amp1/amp0, (i - nuniform + n - 2.)/(n - 3.)));
    mcamp[nuniform - 1] = (real) amp1;
  }
  mcamp[nsystem - 1] = (real) amp2;
  for (i = 0; i < nsystem; i++) mcamp[i] = mcamp[i]/D;

  /* default weights */
  if (Z1 < 0) {
    static double Z1arr[] = {1, 1, 1, 1, 0.29, 0.084, 0.044, 0.032, 0.026, 0.022, 0.019};
    Z1 = (n > 10) ? 0.2/n : Z1arr[n];
  }
  if (Z2 < 0) {
    static double Z2arr[] = {1, 1, 1, 1, 0.53, 0.3, 0.2, 0.17, 0.16};
    Z2 = (n > 8) ? 0.15 : Z2arr[n];
  }

  /* inverse weights / partition functions */
  Z[0] = 1;
  if (nuniform > 1) {
    for (i = 1; i < nuniform - 1; i++)
      Z[i] = pow(Z1, pow((i - nuniform + n - 2.)/(n - 3.0), 2.0));
    Z[nuniform - 1] = Z1;
  }
  Z[nsystem - 1] = Z2;

  /* the smaller system converges faster, so we spend less time on it */
  if (nsteps2 <= 0) nsteps2 = nsteps / 4;

  argopt_dump(ao);

  /* print out the MC amplitude and partition function */
  for (i = 0; i < nsystem; i++)
    printf("%2d: amp %.6f Z %.7f\n", i, mcamp[i], Z[i]);
  argopt_close(ao);
}



/* randomly displace a particle, and update the diagram */
INLINE int disp(rvn_t *x, int n, rvn_t xi, double amp,
    const dg_t *g, dg_t * RESTRICT ng, int usegrand)
{
  int i = (int) (rnd0() * n), j;

  if (usegrand)
    rvn_granddisp(xi, x[i], amp);
  else
    rvn_rnddisp(xi, x[i], amp);
  dg_copy(ng, g);
  for (j = 0; j < n; j++) {
    if (j == i) continue;
    if (rvn_dist2(xi, x[j]) >= 1)
      dg_unlink(ng, i, j);
    else
      dg_link(ng, i, j);
  }
  return i;
}



/* compute the probability of biconnectivity after adding and
 * removing a particle */
INLINE int nmove(rvn_t *x, int n, rvn_t xi,
    const dg_t *g, dg_t *sg, double *nrat)
{
  int i = -1;

  if (rnd0() < 0.5) {
    /* the rate of biconnectivity after adding a random vertex */
    nrat[0] += 1;
    nrat[1] += rvn_voladd(x, n, xi);
  } else {
    /* the rate of biconnectivity after removing a random vertex */
    nrat[2] += 1;
    nrat[3] += dgmc_nremove(g, sg, n, &i);
  }
  return i;
}



#define WTSYS2(wt, fb) { if ((wt = abs(fb)) == 0) wt = 1; }

/* change `sys' between `nsys' systems of biconnected diagrams
 * 1. all
 * 2. those with at least nsys - 2 - i clique separators
 * 3. those with |fb| as the weight
 * assuming fb is known */
INLINE int movesys_heatbath(int sys, int nsys, int fb, int ncsep,
    double Z[], double tacc[])
{
  static double proba[DG_NMAX + 1];
  double r;
  int i, wt;

  /* compute the relative weights */
  proba[0] = 0;
  proba[1] = 1./Z[0]; /* system 0 accepts all diagrams */
  for (i = 1; i < nsys - 1; i++)
    proba[i + 1] = (ncsep <= nsys - 2 - i) ? 1./Z[i] : 0;
  WTSYS2(wt, fb);
  proba[nsys] = (ncsep == 0) ? wt/Z[nsys - 1] : 0;
  /* compute the cumulative probabilities */
  for (i = 0; i < nsys; i++)
    proba[i + 1] += proba[i];

  r = rnd0() * proba[nsys];
  for (i = nsys - 1; i >= 0; i--)
    if (r > proba[i]) break;

  /* compute the acceptance probabilities */
  if (sys == nsys - 2) {
    tacc[0] += 1;
    if (i == nsys - 1) tacc[1] += 1;
  } else if (sys == nsys - 1) {
    tacc[2] += 1;
    if (i == nsys - 2) tacc[3] += 1;
  }
  return i;
}



/* change `sys' between `nsys' systems of biconnected diagrams
 * 1. all
 * 2. those with at least nsys - 2 - i clique separators
 * 3. those with |fb| as the weight
 * Metropolis version */
INLINE int movesys_metro(int sys, int nsys, int fb, int ncsep, int ifb,
    double Z[], double tacc[])
{
  int k, sys1, nwt, wt;
  double r;

  k = (int) (rnd0() * (nsys - 1));
  if (k >= sys) k++;
  die_if (k == sys || k >= nsys, "bad transition k %d, nsys %d\n", k, nsys);
  sys1 = sys;

  if (sys == nsys - 1) { /* importance sampler */
    WTSYS2(wt, fb);
    nwt = 1;
    ifb = 1; /* importance sampler always has fb */
  } else if (k == nsys - 1) { /* sys != nsys - 1 */
    wt = 1;
    WTSYS2(nwt, fb);
    if (ncsep != 0) ifb = 0;
  } else { /* k, sys != nsys - 1 */
    wt = nwt = 1;
    ifb = (k == 0) ? 1 : (ncsep <= nsys - 2 - k);
  }

  if ( ifb )
    if ((r = nwt/Z[k] * Z[sys]/wt) >= 1 || rnd0() < r)
      sys1 = k;

  /* compute the acceptance probabilities */
  if (sys == nsys - 2 && k == nsys - 1) {
    tacc[0] += 1;
    if (k == sys1) tacc[1] += 1;
  } else if (sys == nsys - 1 && k == nsys - 2) {
    tacc[2] += 1;
    if (k == sys1) tacc[3] += 1;
  }
  return sys1;
}



/* compute the sign of the virial coefficient,
 * assuming a lookup table */
INLINE void mcrat_lookup(int n, double nequil, double nsteps,
    double amp[], int usegrand,
    int nstnmv, double *nrat, double *tacc, double *fbav)
{
  rvn_t *x, xi;
  int i, it, fb = 0, nfb = 0, wt, nwt, ncl, nncl, nbc, eql;
  dg_t *g, *ng, *sg;
  double t, nzsm = 0, *fbsm, *hist, *cacc, *ctot;
  int sys, sys0, acc1;
  unqid_t gmapid;
  code_t ncode;

  die_if (n > DGMAP_NMAX, "no diagram map for n %d\n", n);
  xnew(x, n + 1);
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.5 / sqrt(D)), (real) (0.5 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  sg = dg_open(n - 1);
  mkgraph(g, x);
  die_if (!dg_biconnected_lookup(g), "initial diagram not biconnected D %d\n", D);
  fb = dg_hsfb_lookup(g);
  ncl = dg_ncsep_lookup(g);
  sys = 0;
  eql = 1;
  nrat[0] = nrat[2] = 1e-6; nrat[1] = nrat[3] = 0;
  tacc[0] = tacc[2] = 1e-6; tacc[1] = tacc[3] = 0;
  xnew(fbsm, nsystem);
  xnew(hist, nsystem);
  xnew(cacc, nsystem);
  xnew(ctot, nsystem);
  for (i = 0; i < nsystem; i++) {
    fbsm[i] = 0;
    hist[i] = 1e-6;
    cacc[i] = 0;
    ctot[i] = 1e-6;
  }

  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    /* randomly displace a particle */
    i = disp(x, n, xi, amp[sys], g, ng, usegrand);

    gmapid = dg_getmapidx(ng, &ncode);
    nbc = dg_biconnected_lookuplow(n, gmapid);
    nncl = dg_ncsep_lookuplow(ng, ncode);
    nfb = dg_hsfb_lookuplow(n, gmapid);

#define GETACC() { \
    acc1 = 0; \
    if (nbc) { \
      if (sys == 0) { \
        acc1 = 1; \
      } else if (sys < nsystem - 1) { \
        acc1 = (nncl <= nsystem - 2 - sys); \
      } else { /* sys == nsystem - 1 */ \
        if (nncl == nsystem - 1) { \
          WTSYS2(wt, fb); \
          WTSYS2(nwt, nfb); \
          if (nwt >= wt) acc1 = 1; \
          else acc1 = (rnd0() < 1.*nwt/wt); \
        } \
      } \
    } }
    GETACC();

    if ( acc1 ) { /* accept the move */
      rvn_copy(x[i], xi);
      dg_copy(g, ng);
      fb = nfb;
      ncl = nncl;
    }

    /* change `sys' */
    sys0 = sys;
    sys = movesys_heatbath(sys0, nsystem, fb, ncl, Z, tacc);
    //sys = movesys_metro(sys0, nsystem, fb, ncl, 1, Z, tacc);
/* // simplest transition probability for debugging
    Z[nuniform - 1] = 1;
    if (sys0 == 0 && ncl == 0) sys = nuniform - 1;
    else if (sys0 == nuniform - 1) sys = 0; */

#ifdef CHECK
    /* verify if fb has been correctly computed */
    {
      int fb1 = dg_hsfb(g), nfb1 = dg_hsfb(ng);
      int ncl1 = dg_ncsep(g), nncl1 = dg_ncsep(ng);
      if (ncl1 != ncl || fb != fb1
       || (ncl1 != 0 && fb1 != 0)
       || (sys < nsystem - 1 && ncl1 + sys + 2 > nsystem) ) {
        printf("sys %d->%d, ncl %d (%d), nncl %d(%d), fb %d (%d), acc %d, nfb %d (%d)\n",
            sys0, sys, ncl, ncl1, nncl, nncl1, fb, fb1, acc1, nfb, nfb1);
        dg_print(g);
        exit(1);
      }
      if (n <= 5 && sys == nsystem - 2 && (fb1 == 0 || ncl1 != 0)) {
        printf("abnormal last system %d, ncl %d, fb1 %d\n", sys, ncl, fb);
      }
    }
#endif

    if (eql) {
      if (t >= nequil) {
        printf("equilibrated at t %g, nstfb %d\n", t, nstfb);
        t = 0;
        it = 0;
        eql = 0;
      }
      continue;
    }

    /* accumulate data after equilibration */
    cacc[sys0] += acc1;
    ctot[sys0] += 1;
    if (sys == 0) nzsm += (ncl == 0);
    hist[sys] += 1;
    if (sys < nsystem - 1) {
      fbsm[sys] += fb;
    } else {
      WTSYS2(wt, fb);
      fbsm[sys] += fb/wt;
    }

    if (it % nstnmv == 0) { /* compute the acceptance rates */
      rvn_rmcom(x, n); /* remove the origin to the center of mass */
      if (sys == 0) nmove(x, n, xi, g, sg, nrat);
    }

#define REPORT() { \
      int k; \
      double ta1 = tacc[1]/tacc[0], ta2 = tacc[3]/tacc[2], fbsm, hsm; \
      fbav[0] = fbsm[nuniform]/hist[nuniform] * (ta1/(ta2 + 1e-12)) \
              * (Z[nuniform]/Z[nuniform - 1]) * (nzsm/hist[0]); \
      /* multiple histogram: fb = sum_k fbsm[k] / sum_k Z0/Zk hist[k] */ \
      /* and Zk ~ Z[k] hist[k], where Z[k] is runtime weight */ \
      for (fbsm = hsm = 0, k = 0; k < nuniform; k++) \
      { fbsm += fbsm[k]; hsm += hist[0]*Z[0]/Z[k]; } \
      fbav[1] = fbsm/hsm; /* from uniform sampling */ \
      for (k = 0; k < nuniform; k++) { \
        printf("%d ncsep %2d, fb%+12.8f, acc %7.4f, hist %9.6f\n", \
            k, nsystem - 2 - k, fbsm[k]/hist[0] * Z[k]/Z[0], \
            cacc[k]/ctot[k], hist[k]/hist[0]); } \
      printf("D %d, n %d, t%11g, fb %+.8f %+.8f%+7d, " \
          "nz %.7f; tacc %.7f %.7f; n+: %.7f, n-: %.7f\n", \
          D, n, t, fbav[0], fbav[1], fb, nzsm/hist[0], \
          ta1, ta2, nrat[1]/nrat[0], nrat[3]/nrat[2]); \
    }

    if (it % nstrep == 0) {
      it = 0;
      REPORT();
    }
  }
  if (it % nstrep != 1) REPORT();
  free(x);
  dg_close(g);
  dg_close(ng);
  dg_close(sg);
}



#include "dgque.h"

/* compute ncsep with a lookup table */
INLINE int dg_ncsep_direct(const dg_t *g, int nbc)
{
  if (nbc) {
    static dgque_t *q[DG_NMAX + 1];
    int n = g->n, k, ncl = 0;
    code_t c[4];

    die_if(n > 16, "bad n %d\n", n);
    if (!q[n]) {
      k = (n * (n - 1)/2 + 31) / 32;
      q[n] = dgque_open(cachesize, k);
    }
    dg_encode(g, c);
    if (dgque_find(q[n], c, &ncl) < 0) {
      ncl = dg_ncsep(g);
      dgque_add(q[n], c, ncl);
    }
    return ncl;
  } else return g->n - 2;
}



/* compute fb with a lookup table */
INLINE int dg_hsfb_direct(const dg_t *g, int ncl)
{
  if (ncl == 0) {
    static dgque_t *q[DG_NMAX + 1];
    int n = g->n, k, fb = 0;
    code_t c[4];

    die_if(n > 16, "bad n %d\n", n);
    if (!q[n]) {
      k = (n * (n - 1)/2 + 31) / 32;
      q[n] = dgque_open(cachesize, k);
    }
    dg_encode(g, c);
    if (dgque_find(q[n], c, &fb) < 0) {
      fb = dg_hsfbrjw(g);
      dgque_add(q[n], c, fb);
    }
    return fb;
    //return dg_hsfbrjw(g);
  } else return 0;
}



/* compute the sign of the virial coefficient by importance sampling
 * direct version without lookup table */
INLINE void mcrat_direct(int n, double nequil, double nsteps,
    double amp[], int usegrand, int nstfb,
    int nstnmv, double *nrat, double *tacc, double *fbav)
{
  rvn_t *x, xi;
  int i, it, fb = 0, nfb = 0, wt, nwt, ncl, nncl, nbc, eql;
  dg_t *g, *ng, *sg;
  double t, nzsm = 0, *fbsm, *hist, *cacc, *ctot;
  int sys, sys0, acc1;
  int needfb = 0; /* need to compute fb */
  int hasfb = 1; /* if fb has been computed for the step */
  int hasnfb = 1; /* if fb has been computed for the MC trial */
  int ifb = 1; /* if fb is randomly computed */

  xnew(x, n + 1);
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.5 / sqrt(D)), (real) (0.5 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  sg = dg_open(n - 1);
  mkgraph(g, x);
  die_if (!dg_biconnected(g), "initial diagram not biconnected D %d\n", D);

  ncl = dg_ncsep(g);
  fb = dg_hsfb_direct(g, 0);
  hasfb = 1;
  sys = 0;
  eql = 1;
  nrat[0] = nrat[2] = 1e-6; nrat[1] = nrat[3] = 0;
  tacc[0] = tacc[2] = 1e-6; tacc[1] = tacc[3] = 0;
  xnew(fbsm, nsystem);
  xnew(hist, nsystem);
  xnew(cacc, nsystem);
  xnew(ctot, nsystem);
  for (i = 0; i < nsystem; i++) {
    fbsm[i] = 0;
    hist[i] = 1e-6;
    cacc[i] = 0;
    ctot[i] = 1e-6;
  }

  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    /* randomly displace a particle */
    i = disp(x, n, xi, amp[sys], g, ng, usegrand);

    nbc = dg_biconnected(ng);
    /* since computing fb is expensive, we try to find the number of
     * nodes in the clique-separator decomposition
     * if this is zero, fb is very likely nonzero
     * we don't care the value if it is not biconnected  */
    //nncl = nbc ? dg_ncsep(ng) : (n - 2);
    /* the lookup table for ncsep is useful only for n <= 8
     * for larger n, the runtime of dg_ncsep() is negligible
     * in comparison to that of dg_hsfb(), so the gain is < 4% */
    nncl = dg_ncsep_direct(ng, nbc);

    nfb = 0;
    /* we need to compute `nfb' if `sys' == 3 */
    hasnfb = needfb = (sys > 1 && !eql);
    if ( needfb ) { /* if we must compute the new fb */
      /* if the connectivity is unchanged, use the old fb */
      if (hasfb && g->c[i] == ng->c[i])
        nfb = fb;
      else
        nfb = dg_hsfb_direct(ng, nncl);
    }

    GETACC();

    if ( acc1 ) { /* accept the move */
      rvn_copy(x[i], xi);
      dg_copy(g, ng);
      ncl = nncl;
      hasfb = hasnfb;
      if (hasfb) fb = nfb;
    } /* if rejected, maintain the old `hasfb' value */

    /* compute fb only when necessary */
    ifb = (!eql && (nstfb == 1 || rnd0() < 1./nstfb));
    if ( !hasfb && (ifb || needfb) ) {
      fb = dg_hsfb_direct(g, ncl);
      hasfb = 1;
    }

    if (eql) {
      if (t >= nequil) {
        printf("equilibrated at t %g, nstfb %d\n", t, nstfb);
        t = 0;
        it = 0;
        eql = 0;
        /* we only compute `fb' after equilibration */
        fb = dg_hsfb_direct(g, ncl);
        hasfb = 1;
      }
      continue;
    }

    /* change `sys', we do this after equilibration, so `fb' is not needed */
    sys0 = sys;
    sys = movesys_metro(sys0, nsystem, fb, ncl, ifb, Z, tacc);
/* // simplest transition probability for debugging
    Z[nuniform - 1] = 1;
    if (sys0 == 0 && ncl == 0) sys = nuniform - 1;
    else if (sys0 == nuniform - 1) sys = 0; */

#ifdef CHECK
    /* verify if fb has been correctly computed */
    //if (hasfb) {
    if (!eql && (ifb || sys >= nsystem - 1)) {
      int fb1 = dg_hsfb(g), nfb1 = dg_hsfb(ng);
      int ncl1 = dg_ncsep(g), nncl1 = dg_ncsep(ng);
      if (ncl1 != ncl || fb != fb1
       || (ncl1 != 0 && fb1 != 0) /* when ncl1 > 0, fb1 must be 0 */
       || (sys < nsystem - 1 && ncl1 + sys + 2 > nsystem) ) {
        printf("t %g, sys %d->%d, ncl %d (%d), nncl %d(%d), "
            "fb %d (%d), hasfb %d, acc %d, nfb %d (%d) hasnfb %d\n",
            t, sys0, sys, ncl, ncl1, nncl, nncl1,
            fb, fb1, hasfb, acc1, nfb, nfb1, hasnfb);
        dg_print(g);
        exit(1);
      }
      if (n <= 5 && sys == nsystem - 2 && (fb1 == 0 || ncl1 != 0)) {
        printf("abnormal last system %d, ncl %d, fb1 %d\n", sys, ncl, fb);
        dg_print(g);
        exit(1);
      }
    }
#endif

    /* accumulate data after equilibration */
    cacc[sys0] += acc1;
    ctot[sys0] += 1;
    if (ifb || sys0 >= nsystem - 1) {
      if (sys == 0) nzsm += (ncl == 0);
      hist[sys] += 1;
      if (sys < nsystem - 1) {
        fbsm[sys] += fb;
      } else {
        WTSYS2(wt, fb);
        fbsm[sys] += fb/wt;
      }
    }

    if (it % nstnmv == 0) { /* compute the acceptance rates */
      rvn_rmcom(x, n); /* remove the origin to the center of mass */
      if (sys == 0) nmove(x, n, xi, g, sg, nrat);
    }

    if (it % nstrep == 0) {
      it = 0;
      REPORT();
    }
  }
  if (it % nstrep != 1) REPORT();
  free(x);
  dg_close(g);
  dg_close(ng);
  dg_close(sg);
}



/* compute the sign of the virial coefficient */
static void mcrat(int n, double nequil, double nsteps,
    double amp[], int usegrand, int nstfb, int lookup,
    int nstnmv, double *nrat, double *tacc, double *fbav)
{
  if (nsteps < nequil) nsteps = nequil;

  printf("D %d, n %d, %s version, %s displacement\n", D, n,
      lookup ? "lookup" : "direct", usegrand ? "Gaussian" : "uniform");
  if (lookup)
    mcrat_lookup(n, nequil, nsteps, amp, usegrand,
        nstnmv, nrat, tacc, fbav);
  else
    mcrat_direct(n, nequil, nsteps, amp, usegrand,
        nstfb, nstnmv, nrat, tacc, fbav);
}



int main(int argc, char **argv)
{
  double fbav1[2], fbav2[2], rvir[2];
  double nrat1[4], nrat2[4], tr1[4], tr2[4], nr1, nr2;
  int i;

  doargs(argc, argv);
  printf("D %d, n %d, nsteps %g/%g, amp %g/%g, nstfb %d, code %d-bit\n",
      D, n, (double) nsteps, nsteps2, mcamp[0], mcamp[1],
      nstfb, (int) sizeof(code_t) * 8);

  mcrat(n, nequil, nsteps, mcamp, usegrand,
      nstfb, lookup, nstnmv, nrat1, tr1, fbav1);
  mcrat(n - 1, nequil, nsteps2, mcamp, usegrand,
      nstfb, lookup, nstnmv, nrat2, tr2, fbav2);
  nr1 = nrat1[3]/nrat1[2];
  nr2 = nrat2[1]/nrat2[0];
  for (i = 0; i < 2; i++)
    rvir[i] = fbav1[i]/fbav2[i] * nr2/nr1 * (1. - n)/(2. - n)/n;
  printf("D %d, n %d, fbav %g/%g, %g/%g; nrat %g/%g = %g\nvir%d/vir%d:\n"
      "    %.7f, %.7f, %.7f (importance sampling)\n"
      "    %.7f, %.7f, %.7f (multiple histogram)\n",
      D, n, fbav1[0], fbav2[0], fbav1[1], fbav2[1],
      nr2, nr1, nr2/nr1, n, n - 1,
      rvir[0], rvir[0] * 2, rvir[0] * ndvol(D),
      rvir[1], rvir[1] * 2, rvir[1] * ndvol(D));

  //mtsave(NULL);
  return 0;
}

