/* compute the ratio of two successive virial coefficients
 * importance sampler coupled to a uniform sampler
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 mcrat.c -lm */

/* somehow float is slower */
/*
#define HAVEREAL 1
typedef float real;
*/

#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"
#include "dgrjw.h"
#include "mcutil.h"

#define NSYS 3

int ncopy = 1; /* number of independent copies */
int n = 7; /* order */
double nequil = 100000; /* number of equilibration steps */
double nsteps = 10000000;
double nsteps2 = 0;
real mcamp[NSYS] = {1.5f, 0.9f, 0.9f};
int gaussdisp = 0; /* normally distributed */
int nstfb = 0; /* interval of evaluting the weight */
int nstnmv = 100; /* frequency of the n move */
int nstrep = 100000000; /* interval of reporting */
int lookup = -1; /* if to use the lookup table */
double Z[NSYS] = {1, -1, -1}; /* inverse weights (etimated partition function) */
int cachesize = 5; /* number of recently visited diagrams */
double ratcr = 0; /* rate of coordinates replacement */
double Zb = 0, Zbvir = 1;
real nrc = (real) 1; /* radius for particle moves */



/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;
  int i;

  ao = argopt_open(0);
  argopt_add(ao, "-m", "%d", &ncopy, "number of independent copies");
  argopt_add(ao, "-n", "%d", &n, "order n");
  argopt_add(ao, "-0", "%lf", &nequil, "number of equilibration steps");
  argopt_add(ao, "-1", "%lf", &nsteps, "number of simulation steps");
  argopt_add(ao, "-2", "%lf", &nsteps2, "number of simulation steps, for n - 1");
  argopt_add(ao, "-a", "%r", &mcamp[0], "MC amplitude for system 0 (biconnected diagrams)");
  argopt_add(ao, "-A", "%r", &mcamp[1], "MC amplitude for system 1 (nonzero biconnected diagrams)");
  argopt_add(ao, "--A2", "%r", &mcamp[2], "MC amplitude for system 2 (importance sampler)");
  argopt_add(ao, "-Z", "%lf", &Z[1], "relative partition function of system 1");
  argopt_add(ao, "--Z2", "%lf", &Z[2], "relative partition function of system 2");
  argopt_add(ao, "-L", "%d", &lookup, "lookup table (-1: automatic)");
  argopt_add(ao, "-w", "%d", &nstfb, "interval of evaluating the weight");
  argopt_add(ao, "-y", "%d", &nstnmv, "interval of the n-move");
  argopt_add(ao, "-q", "%d", &nstrep, "interval of reporting");
  argopt_add(ao, "-G", "%b", &gaussdisp, "normally-distributed displacement");
  argopt_add(ao, "-C", "%d", &cachesize, "cache size of recently visited diagrams");
  argopt_add(ao, "-R", "%lf", &ratcr, "rate of coordinates replacement");
  argopt_add(ao, "-c", "%r", &nrc, "radius of n-move");
  argopt_add(ao, "-B", "%lf", &Zb, "partition function of biconnected configurations");
  argopt_parse(ao, argc, argv);

  /* decide whether to use lookup table or not */
  if (lookup < 0)
    lookup = (n <= DGMAP_NMAX);

  /* interval of computing fb */
  if (nstfb <= 0) {
    if (lookup) {
      nstfb = 1;
    } else { /* TODO: improve this */
      if (D >= 15) nstfb = 10;
      else if (D >= 10) nstfb = 100;
      else for (nstfb = 10, i = 8; i < n; i++) nstfb *= 3;
    }
  }

  /* Monte Carlo move amplitude */
  for (i = 0; i < NSYS; i++) mcamp[i] /= D;

  /* default weights */
  if (Z[1] < 0) {
    static double Z1[] = {1, 1, 1, 1, 0.29*D/3, 0.084*D/3, 0.044*D/3, 0.032*D/3, 0.026*D/3, 0.022*D/3, 0.019*D/3};
    Z[1] = (n > 10) ? 0.2*D/3/n : Z1[n];
  }
  if (Z[2] < 0) {
    static double Z2[] = {1, 1, 1, 1, 0.53, 0.3, 0.2, 0.17, 0.16, 0.158};
    Z[2] = (n > 9) ? 0.15 : Z2[n];
  }

  /* the smaller system converges faster, so we spend less time on it */
  //if (nsteps2 <= 0) nsteps2 = nsteps / 4;

  //if (ratcr < 0) ratcr = 0.1/n;

  /* if Zb is not given, guess it */
  if (Zb <= 0)
    Zb = getZrat(D, n, NULL);
  if (Zb > 0) {
    /* Zbvir = Zb * (1 - n) / n! */
    Zbvir = Zb * (1 - n);
    for (i = 2; i <= n; i++)
      Zbvir *= 1./i;
  }

  argopt_dump(ao);
  argopt_close(ao);
}



/* compute the probability of biconnectivity after adding and
 * removing a particle */
static void nmove(rvn_t *x, int n, rvn_t xi, const dg_t *g,
    double *nrat, real rcp, real rcm)
{
  if (rnd0() < 0.5) {
    /* the rate of biconnectivity after adding a random vertex */
    nrat[0] += 1;
    nrat[1] += rvn_voladd_dock(x, n, xi, rcp);
  } else {
    /* the rate of biconnectivity after removing a random vertex */
    nrat[2] += 1;
    nrat[3] += dgmc_nremove_dock(g, x, n, rcm);
  }
}



#define WTSYS2(wt, fb) if ((wt = abs(fb)) == 0) wt = 1;

/* change `sys' between three systems of biconnected diagrams
 * 1. all
 * 2. those without clique separators
 * 3. those with |fb| as the weight  */
INLINE int movesys(int sys, int fb, int nz, int ifb,
    const double *Z, double * RESTRICT tacc)
{
  double r;
  int wt;

  if (sys == 0) {
    if ( nz ) { /* only a `nz' configuration can switch */
      if (rnd0() < 0.5) { /* switch sys from 0 to 1 */
        if ((r = Z[0]/Z[1]) >= 1 || rnd0() < r)
          sys = 1;
      } else { /* switch sys from 0 to 2 */
        if ( ifb ) {
          WTSYS2(wt, fb);
          if ((r = wt*Z[0]/Z[2]) >= 1 || rnd0() < r)
            sys = 2;
        }
      }
    }
  } else if (sys == 1) {
    if (rnd0() < 0.5) { /* switch sys from 1 to 0 */
      if ((r = Z[1]/Z[0]) >= 1 || rnd0() < r)
        sys = 0;
    } else { /* switch sys from 1 to 2 */
      if ( ifb ) {
        tacc[0] += 1;
        WTSYS2(wt, fb);
        if ((r = wt*Z[1]/Z[2]) >= 1 || rnd0() < r) {
          sys = 2;
          tacc[1] += 1;
        }
      }
    }
  } else { /* sys 2 */
    if (rnd0() < 0.5) { /* switch sys from 2 to 1 */
      tacc[2] += 1;
      WTSYS2(wt, fb);
      if ((r = Z[2]/(Z[1]*wt)) >= 1 || rnd0() < r) {
        sys = 1;
        tacc[3] += 1;
      }
    } else { /* switch sys from 2 to 0 */
      WTSYS2(wt, fb);
      if ((r = Z[2]/(Z[0]*wt)) >= 1 || rnd0() < r)
        sys = 0;
    }
  }
  return sys;
}



/* change `sys' between three systems of biconnected diagrams
 * 1. all
 * 2. those without clique separators
 * 3. those with |fb| as the weight */
INLINE int movesys_heatbath(int sys, int fb,
    const double *Z, double * RESTRICT tacc)
{
  double proba[4], r;
  int i, wt;

  /* compute the relative weights */
  proba[0] = 0;
  proba[1] = 1./Z[0]; /* system 0 accepts all diagrams */
  proba[2] = proba[1] + 1./Z[1];
  WTSYS2(wt, fb);
  proba[3] = proba[2] + wt/Z[2];

  r = rnd0() * proba[3];
  if (r > proba[2]) i = 2;
  else if (r > proba[1]) i = 1;
  else i = 0;

  /* compute the acceptance probabilities */
  if (sys == 1) {
    tacc[0] += 1;
    tacc[1] += (i == 2);
  } else if (sys == 2) {
    tacc[2] += 1;
    tacc[3] += (i == 1);
  }
  return i;
}




/* compute the sign of the virial coefficient,
 * assuming a lookup table */
static void mcrat_lookup(int n, double nequil, double nsteps,
    real amp[], int gaussdisp,
    int nstnmv, double *nrat, double *tacc, double *fbav)
{
  rvn_t *x, *xm, *nx, xi;
  int i, j, pid, m, it, *fb, nfb = 0, wt, nwt, *nz, nnz, nbc;
  dg_t *g, *ng;
  code_t *code, ncode;
  double *nzsm, (*fbsm)[NSYS], (*hist)[NSYS];
  double t, cacc[NSYS], ctot[NSYS], racc, rtot;
  int *sys, *sys0, acc1;
  unqid_t gmapid;
#ifdef UNITVOL
  real rcp = nrc, rcm = nrc;
#else
  real rcp = nrc * (real) pow(trialvol(n, D), 1./D);
  real rcm = nrc * (real) pow(trialvol(n - 1, D), 1./D);
#endif

  die_if (n > DGMAP_NMAX, "no diagram map for n %d\n", n);
  xnew(x, ncopy * n);
  xnew(nx, n);
  xnew(code, ncopy);
  xnew(fb, ncopy);
  xnew(nz, ncopy);
  xnew(sys, ncopy);
  xnew(sys0, ncopy);
  xnew(nzsm, ncopy);
  xnew(fbsm, ncopy);
  xnew(hist, ncopy);
  g = dg_open(n);
  ng = dg_open(n);
  for (m = 0; m < ncopy; m++) {
    xm = x + m * n;
    initx(xm, n);
    mkgraph(g, xm, n);
    die_if (!dg_biconnected(g),
        "m %d, initial diagram not biconnected D %d\n", m, D);
    /* initialize the code */
    dg_encode(g, &code[m]);
    fb[m] = dg_hsfb_lookup(g);
    nz[m] = 1; /* if (lookup), nz means fb != 0,
               otherwise nz means no clique separator */
    sys[m] = 0;
    nzsm[m] = 0;
    for (i = 0; i < NSYS; i++) {
      fbsm[m][i] = 0;
      hist[m][i] = 1e-6;
    }
  }
  nrat[0] = nrat[2] = 1e-6; nrat[1] = nrat[3] = 0;
  tacc[0] = tacc[2] = 1e-6; tacc[1] = tacc[3] = 0;
  for (i = 0; i < NSYS; i++) {
    cacc[i] = 0;
    ctot[i] = 1e-6;
  }
  racc = 0;
  rtot = 1e-6;

  /* equilibration */
  for (t = 1; t <= nequil; t += 1)
    for (m = 0; m < ncopy; m++) {
      BCSTEP(acc1, i, n, g, ng, x + m*n, xi, amp[0], gaussdisp);
    }
  printf("equilibrated at t %g, nstfb %d, rc+ %g, rc- %g\n", nequil, nstfb, rcp, rcm);
  for (m = 0; m < ncopy; m++) {
    mkgraph(g, x + m*n, n);
    dg_encode(g, &code[m]);
    gmapid = dgmap_[n].map[code[m]];
    nbc = dg_biconnected_lookuplow(n, gmapid);
    fb[m] = nbc ? dg_hsfb_lookuplow(n, gmapid) : 0;
    nz[m] = (fb[m] != 0);
  }


  /* main loop */
  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    for (m = 0; m < ncopy; m++) {
      xm = x + m*n;
      sys0[m] = sys[m];

      if (sys[m] == 0 && rnd0() < ratcr) { /* particle replacement */
        rtot += 1;
        dg_decode(g, &code[m]);
        /* regenerate a configuration */
        if ( grepl(xm, nx, g, ng) ) {
          fb[m] = dg_hsfb_lookup(g);
          nz[m] = (fb[m] != 0);
          racc += 1;
        }
      } else {
        /* randomly displace a particle */
        i = (int) (rnd0() * n);
        if (gaussdisp)
          rvn_granddisp(xi, xm[i], amp[sys0[m]]);
        else
          rvn_rnddisp(xi, xm[i], amp[sys0[m]]);

        /* directly change the connectivity bitstring of the graph */
        ncode = code[m];
        /* for j < i pairs */
        for (pid = i - 1, j = 0; j < i; j++, pid += n - j - 1)
          if (rvn_dist2(xi, xm[j]) >= 1)
            ncode &= ~(1u << pid);
          else
            ncode |= (1u << pid);
        /* for j > i pairs */
        for (j = i + 1, pid = n*i - j*i/2; j < n; j++, pid++)
          if (rvn_dist2(xi, xm[j]) >= 1)
            ncode &= ~(1u << pid);
          else
            ncode |= (1u << pid);

        if (ncode == code[m]) { /* topology unchanged */
          acc1 = 1;
          rvn_copy(xm[i], xi);
        } else {
          gmapid = dgmap_[n].map[ncode];
          nbc = dg_biconnected_lookuplow(n, gmapid);
          /* the weight `fb' is computed from the lookup table
           * for a small n, so it can be updated in every step */
          nfb = nbc ? dg_hsfb_lookuplow(n, gmapid) : 0;
          nnz = (nfb != 0);

          if (sys0[m] == 0) {
            acc1 = nbc;
          } else if (sys0[m] == 1) {
            acc1 = nnz;
          } else { /* sys0[m] == 2 */
            if (nnz) {
              WTSYS2(wt, fb[m]);
              WTSYS2(nwt, nfb);
              if (nwt >= wt) acc1 = 1;
              else acc1 = (rnd0() < 1.*nwt/wt);
            } else acc1 = 0;
          }

          if ( acc1 ) { /* accept the move */
            rvn_copy(xm[i], xi);
            code[m] = ncode;
            nz[m] = nnz;
            fb[m] = nfb;
          }
        }

        /* change `sys', this costs about 16% of simulation time */
        if (nz[m]) /* if fb == 0, then it belongs to system 0 */
          sys[m] = movesys_heatbath(sys0[m], fb[m], Z, tacc);

        /* the cost of accumulating acc is negligible */
        cacc[sys0[m]] += acc1;
        ctot[sys0[m]] += 1;

#ifdef CHECK
        /* verify if fb has been correctly computed */
        {
          int fb1, nfb1;
          dg_decode(g, &code[m]);
          dg_decode(ng, &ncode);
          fb1 = dg_hsfb(g), nfb1 = dg_hsfb(ng);
          if (fb[m] != fb1) {
            printf("t %g, sys %d, nz %d (%d), fb %d (%d), "
                "acc %d, nfb %d (%d)\n",
                t, sys[m], nz[m], fb[m] != 0,
                fb[m], fb1, acc1, nfb, nfb1);
            dg_print(g);
            exit(1);
          }
        }
#endif
      }
    }

    /* accumulate data, doing this every step
     * costs about 5% of the simulation time */
    for (m = 0; m < ncopy; m++) {
      if (sys[m] == 0) nzsm[m] += nz[m];
      hist[m][sys[m]] += 1;
      if (sys[m] <= 1) {
        fbsm[m][sys[m]] += fb[m];
      } else {
        WTSYS2(wt, fb[m]);
        fbsm[m][sys[m]] += fb[m]/wt;
      }
    }

    if (it % nstnmv == 0) { /* compute the n-move rates */
      for (m = 0; m < ncopy; m++) {
        xm = x + m*n;
        rvn_rmcom(xm, n); /* move the center of mass to the origin */
        if (sys[m] == 0) {
          dg_decode(g, &code[m]);
          nmove(xm, n, xi, g, nrat, rcp, rcm);
        }
      }
    }

#define REPORT(gmapid) { \
      double ta1 = tacc[1]/tacc[0], ta2 = tacc[3]/tacc[2]; \
      double hist0 = 0, hist1 = 0, hist2 = 0, nzsm0 = 0; \
      double fbsm0 = 0, fbsm1 = 0, fbsm2 = 0; \
      double fb0, dfb0 = 0, fb0sm = 0, fb0sm2 = 0; \
      double fb1, dfb1 = 0, fb1sm = 0, fb1sm2 = 0; \
      for (m = 0; m < ncopy; m++) { \
        nzsm0 += nzsm[m]; \
        hist0 += hist[m][0]; \
        hist1 += hist[m][1]; \
        hist2 += hist[m][2]; \
        fbsm0 += fbsm[m][0]; \
        fbsm1 += fbsm[m][1]; \
        fbsm2 += fbsm[m][2]; \
        fb0 = fbsm[m][2]/hist[m][2]*ta1/ta2*Z[2]/Z[1]*nzsm[m]/hist[m][0]; \
        fb1 = (fbsm[m][0]+fbsm[m][1])/hist[m][0]/(1 + hist[m][1]/nzsm[m]); \
        fb0sm += fb0; fb0sm2 += fb0 * fb0; \
        fb1sm += fb1; fb1sm2 += fb1 * fb1; \
        if (ncopy > 1) \
          printf("m %d, %+.6f/%+.6f|%+.6f %+.6f %d\n", m, \
            fbsm[m][0]/hist[m][0], fbsm[m][1]/hist[m][0]*Z[1]/Z[0], \
            fb1, fb0, fb[m]); \
      } \
      if (ncopy > 1) { /* compute the standard deviation */ \
        fb0sm /= ncopy; \
        dfb0 = sqrt((fb0sm2/ncopy - fb0sm*fb0sm)/(ncopy - 1)); \
        fb1sm /= ncopy; \
        dfb1 = sqrt((fb1sm2/ncopy - fb1sm*fb1sm)/(ncopy - 1)); \
      } \
      fbav[0] = fbsm2/hist2 * (ta1/ta2)*(Z[2]/Z[1]) * (nzsm0/hist0); \
      /* multiple histogram: fb = sum_k fbsm[k] / sum_k Z0/Zk hist[k] */ \
      fbav[1] = (fbsm0 + fbsm1) / (hist0 + hist1*hist0/nzsm0); \
      printf("D %d, n %d, t%11g, cacc %.4f/%.4f/%.4f, " \
          "fb %+.6f/%+.6f|%+.6f(%.6f) %+.6f(%.6f)\n" \
          "    nz %.7f, his1/0 %.7f, 2/1 %.7f; " \
          "tacc %.7f %.7f; n+: %.7f, n-: %.7f, racc %.6f id %d\n", \
          D, n, t, cacc[0]/ctot[0], cacc[1]/ctot[1], cacc[2]/ctot[2], \
          fbsm0/hist0, fbsm1/hist0*Z[1]/Z[0], fbav[1], dfb1, fbav[0], \
          dfb0, nzsm0/hist0, hist1/hist0, hist2/hist1, ta1, ta2, \
          nrat[1]/nrat[0], nrat[3]/nrat[2], \
          racc/rtot, (int) gmapid); \
    }

    if (it % nstrep == 0) {
      it = 0;
      REPORT(gmapid);
    }
  }
  if (it % nstrep != 1) REPORT(gmapid);
  free(x);
  free(nx);
  free(code);
  dg_close(g);
  dg_close(ng);
  free(fb);
  free(nz);
  free(sys);
  free(sys0);
  free(nzsm);
  free(fbsm);
  free(hist);
}



#include "dgque.h"

/* compute fb with a short history list
 * nocsep: no clique separator in the graph */
INLINE int dg_hsfb_que(dgque_t *q, const dg_t *g, int nocsep,
    int *ned, int *degs, int *neval)
{
  int fb = 0;
  static code_t c[DG_NMAX/2];

  dg_encode(g, c);
  if (dgque_find(q, c, &fb) < 0) {
    /* hsfb_mixed() is much faster than hsfb_rjw() in high dimensions D
     * because when the number of edges ~ n, rhsc() is faster than hsfb() */
    fb = dg_hsfb_mixed0(g, nocsep, ned, degs);
    dgque_add(q, c, fb);
    if (neval) (*neval)++;
  }
  return fb;
}



/* compute the sign of the virial coefficient by importance sampling
 * direct version without lookup table */
INLINE void mcrat_direct(int n, double nequil, double nsteps,
    real amp[], int gaussdisp, int nstfb,
    int nstnmv, double *nrat, double *tacc, double *fbav)
{
  rvn_t **x, *nx, xi;
  int i, m, it, *fb, nfb = 0, wt, nwt, *nz, nnz, nbc;
  dg_t **g, *ng;
  double *nzsm, (*fbsm)[NSYS], (*hist)[NSYS];
  double t, cacc[NSYS], ctot[NSYS], racc, rtot;
  int *sys, *sys0, acc1;
  int *ifb; /* if fb is randomly computed */
  int *hasfb; /* if fb has been computed for the step */
  int hasnfb = 0; /* if fb has been computed for the MC trial */
  int neval[3] = {0, 0, 0};
  dgque_t **que;
#ifdef UNITVOL
  real rcp = nrc, rcm = nrc;
#else
  real rcp = nrc * (real) pow(trialvol(n, D), 1./D);
  real rcm = nrc * (real) pow(trialvol(n - 1, D), 1./D);
#endif

  xnew(x, ncopy);
  xnew(nx, n + 1);
  xnew(g, ncopy);
  xnew(fb, ncopy);
  xnew(ifb, ncopy);
  xnew(hasfb, ncopy);
  xnew(nz, ncopy);
  xnew(sys, ncopy);
  xnew(sys0, ncopy);
  xnew(nzsm, ncopy);
  xnew(fbsm, ncopy);
  xnew(hist, ncopy);
  xnew(que, ncopy);
  ng = dg_open(n);
  for (m = 0; m < ncopy; m++) {
    xnew(x[m], n + 1);
    initx(x[m], n);
    g[m] = dg_open(n);
    mkgraph(g[m], x[m], n);
    die_if (!dg_biconnected(g[m]),
        "m %d, initial diagram not biconnected D %d\n", m, D);
    que[m] = dgque_open(cachesize, (n * (n - 1) / 2 + 31) / 32);
    fb[m] = dg_hsfb_que(que[m], g[m], 0, NULL, NULL, NULL);
    hasfb[m] = 1;
    nz[m] = 1; /* if (lookup), nz means fb != 0,
               otherwise nz means no clique separator */
    sys[m] = 0;
    nzsm[m] = 0;
    for (i = 0; i < NSYS; i++) {
      fbsm[m][i] = 0;
      hist[m][i] = 1e-6;
    }
    ifb[m] = 0;
    hasfb[m] = 0;
  }
  nrat[0] = nrat[2] = 1e-6; nrat[1] = nrat[3] = 0;
  tacc[0] = tacc[2] = 1e-6; tacc[1] = tacc[3] = 0;
  for (i = 0; i < NSYS; i++) {
    cacc[i] = 0;
    ctot[i] = 1e-6;
  }
  racc = 0;
  rtot = 1e-6;

  /* equilibration */
  for (t = 1; t <= nequil; t += 1)
    for (m = 0; m < ncopy; m++) {
      BCSTEP(acc1, i, n, g[m], ng, x[m], xi, amp[0], gaussdisp);
    }
  printf("equilibrated at t %g, nstfb %d, rc+ %g, rc- %g\n", nequil, nstfb, rcp, rcm);
  for (m = 0; m < ncopy; m++) {
    nbc = dg_biconnected(g[m]);
    nz[m] = nbc ? (dg_cliquesep(g[m]) == 0) : 0;
    fb[m] = nz[m] ? dg_hsfb_que(que[m], g[m], 1, NULL, NULL, NULL) : 0;
    hasfb[m] = 1;
  }

  /* main loop */
  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    for (m = 0; m < ncopy; m++) {
      /* save the system, 0: biconnected; 1: nocsep; 2: w ~ 1/|fb| */
      sys0[m] = sys[m];

#ifdef CHECK
      if ( hasfb[m] ) {
        int nz1 = (dg_cliquesep(g[m]) == 0), ncl = dg_ncsep(g[m]);
        int fb1 = dg_hsfb(g[m]);
        if (nz[m] != nz1 || fb[m] != fb1) {
          printf("PRE1 t %g, sys %d, hasfb %d, nz %d (%d), ncl %d, fb %d (%d)\n",
              t, sys[m], hasfb[m], nz[m], nz1, ncl, fb[m], fb1);
          dg_print(g[m]);
          exit(1);
        }
      } else {
        int nz1 = (dg_cliquesep(g[m]) == 0), ncl = dg_ncsep(g[m]);
        int fb1 = dg_hsfb(g[m]);
        if (nz[m] != nz1) {
          printf("PRE2 t %g, sys %d, hasfb %d, nz %d (%d), ncl %d, fb %d (%d)\n",
              t, sys[m], hasfb[m], nz[m], nz1, ncl, fb[m], fb1);
          dg_print(g[m]);
          exit(1);
        }
      }
#endif

      if (sys0[m] == 0 && rnd0() < ratcr) { /* replace all coordinates */
        rtot += 1;
        /* regularly regenerate a configuration */
        if ( grepl(x[m], nx, g[m], ng) ) {
          nz[m] = (dg_cliquesep(g[m]) == 0);
          hasfb[m] = 0;
          fb[m] = 0;
          racc += 1;
        }
      } else {  /* randomly displace a particle */
        DISPRNDI(i, n, x[m], xi, amp[sys0[m]], gaussdisp);
        UPDGRAPH(i, n, g[m], ng, x[m], xi);

        /* try to avoid the expensive computation of fb */
        nbc = dg_biconnected(ng);
        if ( !nbc ) {
          nnz = 0;
          nfb = 0;
          hasnfb = 1;
        } else { /* biconnected */
          int ned = -1;
          static int degs[DG_NMAX];
          int sc, err;

          /* detect special cases, including ring diagrams
           * and clique separators */
          sc = dg_rhsc_spec0(ng, 0, &ned, degs, &err);
          if (err == 0) { /* special case worked */
            hasnfb = 1;
            if (sc != 0) {
              nnz = 1;
              nfb = DG_SC2FB(sc, ned);
            } else {
              nnz = 0;
              nfb = 0;
            }
          } else { /* general case */
            nnz = 1; /* no clique separator */
            /* we need to compute `nfb' if `sys' == 3 */
            if ( sys0[m] >= 2 ) { /* if we must compute the new fb */
              hasnfb = 1;
              /* if the connectivity is unchanged, use the old fb */
              if (hasfb[m] && g[m]->c[i] == ng->c[i]) {
                nfb = fb[m];
              } else {
                nfb = dg_hsfb_que(que[m], ng, 1, &ned, degs, &neval[sys0[m]]);
              }
            } else { /* indicate fb has not been computed */
              hasnfb = nfb = 0;
            }
          } /* end of the general case */
        } /* end if biconnected */

        if (sys0[m] == 0) {
          acc1 = nbc;
        } else if (sys0[m] == 1) {
          acc1 = nnz;
        } else { /* sys0[m] == 2 */
          if (nnz) {
            WTSYS2(wt, fb[m]);
            WTSYS2(nwt, nfb);
            if (nwt >= wt) acc1 = 1;
            else acc1 = (rnd0() < 1.*nwt/wt);
          } else acc1 = 0;
        }

        if ( acc1 ) { /* accept the move */
          rvn_copy(x[m][i], xi);
          dg_copy(g[m], ng);
          nz[m] = nnz;
          hasfb[m] = hasnfb;
          if (hasfb[m]) fb[m] = nfb;
        } /* if rejected, maintain the old `hasfb' value */

        cacc[sys0[m]] += acc1;
        ctot[sys0[m]] += 1;
      }

      /* compute fb only when necessary
       * this usually happens after an accepted move */
      ifb[m] = (nstfb == 1 || rnd0() < 1./nstfb);
      if ( !hasfb[m] && (ifb[m] || sys0[m] >= 2) ) {
        fb[m] = nz[m] ? dg_hsfb_que(que[m], g[m], 1, NULL, NULL, &neval[sys0[m]]) : 0;
        hasfb[m] = 1;
      }

      /* change `sys' */
      sys[m] = movesys(sys0[m], fb[m], nz[m], ifb[m], Z, tacc);

#ifdef CHECK
      /* check if fb has been correctly computed */
      if (hasfb[m]) {
      //if ((ifb[m] || sys[m] >= 2)) {
        int nz1 = (dg_cliquesep(g[m]) == 0), ncl = dg_ncsep(g[m]);
        int fb1 = dg_hsfb(g[m]), nfb1 = dg_hsfb(ng);
        if (!hasfb[m] || nz[m] != nz1 || fb[m] != fb1) {
          printf("t %g, m %d, acc %d, ifb %d, sys %d/%d, hasfb %d, nz %d (%d), "
              "ncl %d, fb %d (%d), hasnfb %d, nfb %d (%d)\n",
              t, m, acc1, ifb[m], sys0[m], sys[m], hasfb[m], nz[m], nz1,
              ncl, fb[m], fb1, hasnfb, nfb, nfb1);
          dg_print(g[m]);
          exit(1);
        }
      }
#endif

      /* accumulate data */
      if (ifb[m] || sys0[m] >= 2) {
        if (sys[m] == 0) nzsm[m] += nz[m];
        hist[m][sys[m]] += 1;
        if (sys[m] <= 1) {
          fbsm[m][sys[m]] += fb[m];
        } else {
          WTSYS2(wt, fb[m]);
          fbsm[m][sys[m]] += fb[m]/wt;
        }
      }

      if (it % nstnmv == 0) { /* compute the acceptance rates */
        rvn_rmcom(x[m], n); /* remove the origin to the center of mass */
        if (sys[m] == 0) nmove(x[m], n, xi, g[m], nrat, rcp, rcm);
      }
    } /* loop over copies */

#define REPORT_DIRECT() { REPORT(0); \
    printf("eval %d/%d/%d\n", neval[0], neval[1], neval[2]); }
    if (it % nstrep == 0) {
      it = 0;
      REPORT_DIRECT();
    }
  }
  if (it % nstrep != 1) REPORT_DIRECT();
  for (m = 0; m < ncopy; m++) {
    dg_close(g[m]);
    free(x[m]);
    dgque_close(que[m]);
  }
  free(x);
  free(nx);
  free(g);
  dg_close(ng);
  free(fb);
  free(ifb);
  free(hasfb);
  free(nz);
  free(sys);
  free(sys0);
  free(nzsm);
  free(fbsm);
  free(hist);
  free(que);
}



/* compute the sign of the virial coefficient */
static void mcrat(int n, double nequil, double nsteps,
    real amp[], int gaussdisp, int nstfb, int lookup,
    int nstnmv, double *nrat, double *tacc, double *fbav)
{
  if (nsteps < nequil) nsteps = nequil;

  printf("D %d, n %d, %s version, %s displacement\n", D, n,
      lookup ? "lookup" : "direct", gaussdisp ? "Gaussian" : "uniform");
  if (lookup)
    mcrat_lookup(n, nequil, nsteps, amp, gaussdisp,
        nstnmv, nrat, tacc, fbav);
  else
    mcrat_direct(n, nequil, nsteps, amp, gaussdisp,
        nstfb, nstnmv, nrat, tacc, fbav);
}



int main(int argc, char **argv)
{
  double fbav1[2], fbav2[2], rvir[2];
  double nrat1[4], nrat2[4], tr1[4], tr2[4], nr1, nr2;
  double vol;
  int i;

  doargs(argc, argv);
  printf("D %d, n %d, nsteps %g/%g, amp %g/%g/%g, nstfb %d, code %d-bit, Z %g/%g\n",
      D, n, (double) nsteps, nsteps2, mcamp[0], mcamp[1], mcamp[2],
      nstfb, (int) sizeof(code_t) * 8, Zb, Zbvir);

  mcrat(n, nequil, nsteps, mcamp, gaussdisp,
      nstfb, lookup, nstnmv, nrat1, tr1, fbav1);

  if (nsteps2 <= 0) /* a single simulation */
  {
    for (i = 0; i < 2; i++)
      rvir[i] = fbav1[i] * Zbvir;
    printf("D %d, n %d, fbav %.8f, %.8f; nrat %g(+), %g(-), Z %g/%g\nvir:\n"
        "    %.12f, %.12f, %.12f (importance sampling)\n"
        "    %.12f, %.12f, %.12f (multiple histogram)\n",
        D, n, fbav1[0], fbav1[1], nrat1[3]/nrat1[2], nrat1[1]/nrat1[0], Zb, Zbvir,
        rvir[0], rvir[0] * pow(2, n - 1), rvir[0] * pow(ndvol(D), n - 1),
        rvir[1], rvir[1] * pow(2, n - 1), rvir[1] * pow(ndvol(D), n - 1));
  }
  else /* a second simulation */
  {
#ifdef UNITVOL
    vol = pow(nrc, D);
#else
    vol = pow(nrc, D) * trialvol(n - 1, D);
#endif
    mcrat(n - 1, nequil, nsteps2, mcamp, gaussdisp,
        nstfb, lookup, nstnmv, nrat2, tr2, fbav2);
    nr1 = nrat1[3]/nrat1[2]; /* remove a vertex n --> n - 1 */
    nr2 = nrat2[1]/nrat2[0]; /* add a vertex, n - 1 --> n */
    for (i = 0; i < 2; i++)
      rvir[i] = fbav1[i]/fbav2[i] * vol * nr2/nr1 * (1. - n)/(2. - n)/n;
    printf("D %d, n %d, fbav %.8f/%.8f, %.8f/%.8f; nrat %g/%g = %g\nvir%d/vir%d:\n"
        "    %.12f, %.12f, %.12f (importance sampling)\n"
        "    %.12f, %.12f, %.12f (multiple histogram)\n",
        D, n, fbav1[0], fbav2[0], fbav1[1], fbav2[1],
        nr2, nr1, nr2/nr1, n, n - 1,
        rvir[0], rvir[0] * 2, rvir[0] * ndvol(D),
        rvir[1], rvir[1] * 2, rvir[1] * ndvol(D));
  }

  mtsave(NULL);
  return 0;
}

