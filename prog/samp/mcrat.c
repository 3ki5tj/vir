/* compute the ratio of two successive virial coefficients
 * importance sampler coupled to a uniform sampler
 * define `D' in compiling to change the default dimension:
 *  gcc -DD=4 -O3 mcrat.c -lm */



#define ZCOM_PICK
#define ZCOM_UTIL
#define ZCOM_ARGOPT
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"
#include "dgrjw.h"
#include "mcutil.h"

#define NSYS 3

int n = 7; /* order */
double nequil = 100000; /* number of equilibration steps */
double nsteps = 10000000, nsteps2 = -1;
real mcamp[NSYS] = {1.5f, 0.9f, 0.9f};
int usegrand = 0; /* normally distributed */
int nstfb = 0; /* interval of evaluting the weight */
int nstnmv = 100; /* frequency of the n move */
int nstrep = 100000000; /* interval of reporting */
int lookup = -1; /* if to use the lookup table */
double Z[NSYS] = {1, -1, -1}; /* inverse weights (etimated partition function) */
int cachesize = 10; /* number of recently visited diagrams */


/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;
  int i;

  ao = argopt_open(0);
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
  argopt_add(ao, "-G", "%b", &usegrand, "normally-distributed displacement");
  argopt_add(ao, "-C", "%d", &cachesize, "cache size of recently visited diagrams");
  argopt_parse(ao, argc, argv);

  /* decide whether to use lookup table or not */
  if (lookup < 0)
    lookup = (n <= DGMAP_NMAX);

  /* interval of computing fb */
  if (nstfb <= 0) {
    if (lookup) {
      nstfb = 1;
    } else {
      for (nstfb = 10, i = 8; i < n; i++) nstfb *= 3;
    }
  }

  /* Monte Carlo move amplitude */
  for (i = 0; i < NSYS; i++) mcamp[i] /= D;

  /* default weights */
  if (Z[1] < 0) {
    static double Z1[] = {1, 1, 1, 1, 0.29, 0.084, 0.044, 0.032, 0.026, 0.022, 0.019};
    Z[1] = (n > 10) ? 0.2/n : Z1[n];
  }
  if (Z[2] < 0) {
    static double Z2[] = {1, 1, 1, 1, 0.53, 0.3, 0.2, 0.17, 0.16};
    Z[2] = (n > 8) ? 0.15 : Z2[n];
  }

  /* the smaller system converges faster, so we spend less time on it */
  if (nsteps2 <= 0) nsteps2 = nsteps / 4;

  argopt_dump(ao);
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



#define WTSYS2(wt, fb) if ((wt = abs(fb)) == 0) wt = 1;

/* change `sys' between three systems of biconnected diagrams
 * 1. all
 * 2. those without clique separators
 * 3. those with |fb| as the weight  */
INLINE int movesys(int sys, int fb, int nz, int ifb,
    double Z[], double tacc[])
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



/* compute the sign of the virial coefficient,
 * assuming a lookup table */
INLINE void mcrat_lookup(int n, double nequil, double nsteps,
    double amp[], int usegrand,
    int nstnmv, double *nrat, double *tacc, double *fbav)
{
  rvn_t *x, xi;
  int i, it, fb = 0, nfb = 0, wt, nwt, nz, nnz, nbc, eql;
  dg_t *g, *ng, *sg;
  double t, fbnz = 0, fbsum[NSYS], hist[NSYS], cacc[NSYS], ctot[NSYS];
  int sys, sys0, acc1;
  unqid_t gmapid;
#ifdef MCLOG
  code_t code[4];
  FILE *fplog = fopen("dg.log", "w");
#endif

  die_if (n > DGMAP_NMAX, "no diagram map for n %d\n", n);
  xnew(x, n + 1);
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.05 / sqrt(D)), (real) (0.05 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  sg = dg_open(n - 1);
  mkgraph(g, x);
  die_if (!dg_biconnected_lookup(g), "initial diagram not biconnected D %d\n", D);
  fb = dg_hsfb_lookup(g);
  nz = 1; /* if (lookup), nz means fb != 0,
             otherwise nz means no clique separator */
  sys = 0;
  eql = 1;
  nrat[0] = nrat[2] = 1e-6; nrat[1] = nrat[3] = 0;
  tacc[0] = tacc[2] = 1e-6; tacc[1] = tacc[3] = 0;
  for (i = 0; i < NSYS; i++) {
    fbsum[i] = 0;
    hist[i] = 1e-6;
    cacc[i] = 0;
    ctot[i] = 1e-6;
  }

  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    /* randomly displace a particle */
    i = disp(x, n, xi, amp[sys], g, ng, usegrand);

    gmapid = dg_getmapid(ng);
    nbc = dg_biconnected_lookuplow(n, gmapid);
    /* the weight `fb' is computed from the lookup table
     * for a small n, so it can be updated in every step */
    nfb = nbc ? dg_hsfb_lookuplow(n, gmapid) : 0;
    nnz = (nfb != 0);

    if (sys == 0) {
      acc1 = nbc;
    } else if (sys == 1) {
      acc1 = nnz;
    } else { /* sys == 2 */
      if (nnz) {
        WTSYS2(wt, fb);
        WTSYS2(nwt, nfb);
        if (nwt >= wt) acc1 = 1;
        else acc1 = (rnd0() < 1.*nwt/wt);
      } else acc1 = 0;
    }

    if ( acc1 ) { /* accept the move */
      rvn_copy(x[i], xi);
      dg_copy(g, ng);
      nz = nnz;
      fb = nfb;
    }

    /* change `sys' */
    sys = movesys(sys0 = sys, fb, nz, 1, Z, tacc);

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
    if (sys == 0) fbnz += nz;
    hist[sys] += 1;
    if (sys <= 1) {
      fbsum[sys] += fb;
    } else {
      WTSYS2(wt, fb);
      fbsum[sys] += fb/wt;
    }

#ifdef MCLOG
    { static code_t c0 = 0;
      dg_encode(g, code);
      gmapid = dg_getmapid(g);
      if (code[0] != c0) {
        fprintf(fplog, "%x %d %d\n", code[0], gmapid, fb);
        c0 = code[0];
      } }
#endif

    if (it % nstnmv == 0) { /* compute the acceptance rates */
      rvn_rmcom(x, n); /* remove the origin to the center of mass */
      if (sys == 0) nmove(x, n, xi, g, sg, nrat);
    }

#define REPORT(gmapid) { \
      double ta1 = tacc[1]/tacc[0], ta2 = tacc[3]/tacc[2]; \
      fbav[0] = fbsum[2]/hist[2] * (ta1/ta2)*(Z[2]/Z[1]) * (fbnz/hist[0]); \
      /* multiple histogram: fb = sum_k fbsum[k] / sum_k Z0/Zk hist[k] */ \
      fbav[1] = (fbsum[0] + fbsum[1]) / (hist[0] + hist[1]*hist[0]/fbnz); \
      fprintf(stderr, "D %d, n %d, t%11g, cacc %.4f/%.4f/%.4f, " \
          "fb %+.8f/%+.8f|%+.8f %+.8f, %d, id %d\n" \
          "    nz %.7f, his1/0 %.7f, 2/1 %.7f; " \
          "tacc %.7f %.7f; n+: %.7f, n-: %.7f\n", \
        D, n, t, cacc[0]/ctot[0], cacc[1]/ctot[1], cacc[2]/ctot[2], \
        fbsum[0]/hist[0], fbsum[1]/hist[0] * Z[1]/Z[0], fbav[1], fbav[0], \
        fb, (int) gmapid, fbnz/hist[0], hist[1]/hist[0], hist[2]/hist[1], \
        ta1, ta2, nrat[1]/nrat[0], nrat[3]/nrat[2]); }

    if (it % nstrep == 0) {
      it = 0;
      REPORT(gmapid);
    }
  }
#ifdef MCLOG
  fclose(fplog);
#endif
  if (it % nstrep != 1) REPORT(gmapid);
  free(x);
  dg_close(g);
  dg_close(ng);
  dg_close(sg);
}



#include "dgque.h"

/* compute fb with a lookup table */
INLINE int dg_hsfb_direct(const dg_t *g, int nz)
{
  if (nz) {
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
  int i, it, fb = 0, nfb = 0, wt, nwt, nz, nnz, nbc, eql;
  dg_t *g, *ng, *sg;
  double t, fbnz = 0, fbsum[NSYS], hist[NSYS], cacc[NSYS], ctot[NSYS];
  int sys, sys0, acc1;
  int needfb = 0; /* need to compute fb */
  int hasfb = 1; /* if fb has been computed for the step */
  int hasnfb = 1; /* if fb has been computed for the MC trial */
  int ifb = 1; /* if fb is randomly computed */

  xnew(x, n + 1);
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.05 / sqrt(D)), (real) (0.05 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  sg = dg_open(n - 1);
  mkgraph(g, x);
  die_if (!dg_biconnected(g), "initial diagram not biconnected D %d\n", D);

  fb = dg_hsfb_direct(g, 1);
  hasfb = 1;
  nz = 1; /* if (lookup), nz means fb != 0,
             otherwise nz means no clique separator */
  sys = 0;
  eql = 1;
  nrat[0] = nrat[2] = 1e-6; nrat[1] = nrat[3] = 0;
  tacc[0] = tacc[2] = 1e-6; tacc[1] = tacc[3] = 0;
  for (i = 0; i < NSYS; i++) {
    fbsum[i] = 0;
    hist[i] = 1e-6;
    cacc[i] = 0;
    ctot[i] = 1e-6;
  }

  if (nsteps < nequil) nsteps = nequil;

  for (it = 1, t = 1; t <= nsteps; t += 1, it++) {
    /* randomly displace a particle */
    i = disp(x, n, xi, amp[sys], g, ng, usegrand);

    nbc = dg_biconnected(ng);
    /* since computing fb is expensive,
     * we only try to find if there is a clique separator,
     * if so, fb == 0, otherwise fb is likely but not always nonzero */
    nnz = nbc ? (dg_cliquesep(ng) == 0) : 0;

    nfb = 0;
    /* we need to compute `nfb' if `sys' == 3 */
    hasnfb = needfb = (sys > 1 && !eql);
    if ( needfb ) { /* if we must compute the new fb */
      /* if the connectivity is unchanged, use the old fb */
      if (hasfb && g->c[i] == ng->c[i])
        nfb = fb;
      else
        nfb = dg_hsfb_direct(ng, nnz);
    }

    if (sys == 0) {
      acc1 = nbc;
    } else if (sys == 1) {
      acc1 = nnz;
    } else { /* sys == 2 */
      if (nnz) {
        WTSYS2(wt, fb);
        WTSYS2(nwt, nfb);
        if (nwt >= wt) acc1 = 1;
        else acc1 = (rnd0() < 1.*nwt/wt);
      } else acc1 = 0;
    }

    if ( acc1 ) { /* accept the move */
      rvn_copy(x[i], xi);
      dg_copy(g, ng);
      nz = nnz;
      hasfb = hasnfb;
      if (hasfb) fb = nfb;
    } /* if rejected, maintain the old `hasfb' value */

    /* compute fb only when necessary */
    ifb = (!eql && (nstfb == 1 || rnd0() < 1./nstfb));
    if ( !hasfb && (ifb || needfb) ) {
      fb = dg_hsfb_direct(g, nz);
      hasfb = 1;
    }

#ifdef CHECK
    /* verify if fb has been correctly computed */
    if (hasfb) {
    //if (!eql && (ifb || sys >= 2)) {
      int nz1 = (dg_cliquesep(g) == 0);
      int fb1 = dg_hsfb(g), nfb1 = dg_hsfb(ng);
      if (!hasfb || nz1 != nz || fb != fb1) {
        printf("ifb %d, sys %d, hasfb %d (1), nz %d (%d), fb %d (%d), "
            "acc %d, hasnfb %d, nfb %d (%d)\n",
            ifb, sys, hasfb, nz, nz1, fb, fb1, acc1, hasnfb, nfb, nfb1);
        dg_print(g);
        exit(1);
      }
    }
#endif

    if (eql) {
      if (t >= nequil) {
        printf("equilibrated at t %g, nstfb %d\n", t, nstfb);
        t = 0;
        it = 0;
        eql = 0;
        /* we only compute `fb' after equilibration */
        fb = dg_hsfb_direct(g, nz);
        hasfb = 1;
      }
      continue;
    }

    /* change `sys', we do this after equilibration
     * so `fb' is not needed */
    sys = movesys(sys0 = sys, fb, nz, ifb, Z, tacc);

    /* accumulate data after equilibration */
    cacc[sys0] += acc1;
    ctot[sys0] += 1;
    if (ifb || sys0 >= 2) {
      if (sys == 0) fbnz += nz;
      hist[sys] += 1;
      if (sys <= 1) {
        fbsum[sys] += fb;
      } else {
        WTSYS2(wt, fb);
        fbsum[sys] += fb/wt;
      }
    }

    if (it % nstnmv == 0) { /* compute the acceptance rates */
      rvn_rmcom(x, n); /* remove the origin to the center of mass */
      if (sys == 0) nmove(x, n, xi, g, sg, nrat);
    }

    if (it % nstrep == 0) {
      it = 0;
      REPORT(0);
    }
  }
  if (it % nstrep != 1) REPORT(0);
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

  mtsave(NULL);
  return 0;
}

