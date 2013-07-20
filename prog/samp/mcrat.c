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

int n = 3; /* order */
double nequil = 100000; /* number of equilibration */
double nsteps = 10000000, nsteps2 = -1;
real mcamp[NSYS] = {1.5f, 0.9f, 0.9f};
int nstfb = 0; /* interval of evaluting the weight */
int nstnmv = 100; /* frequency of the n move */
int nstrep = 100000000; /* interval of reporting */
int lookup = -1; /* if to use the lookup table */
double Z[NSYS] = {1, -1, -1}; /* inverse weights (etimated partition function) */



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
  argopt_add(ao, "-z", "%d", &nstnmv, "interval of the n-move");
  argopt_add(ao, "-q", "%d", &nstrep, "interval of reporting");
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



/* compute the sign of the virial coefficient */
static double mcrun(int n, double nequil, double nsteps,
    double amp[2], int nstfb, int lookup,
    int nstnmv, double *nrat, double *tacc)
{
  rvn_t *x, xi;
  int i, j, it, fb = 0, nfb = 0, ifb, wt, nwt, nz, nnz, nbc, eql;
  dg_t *g, *ng, *sg;
  double fbnz, fbav, t, r;
  double fbsum[NSYS], hist[NSYS], cacc[NSYS], ctot[NSYS];
  int sys, sys0, acc1, hasfb, hasnfb;

  xnew(x, n + 1);
  for (i = 0; i < n; i++)
    rvn_rnd(x[i], (real) (-0.05 / sqrt(D)), (real) (0.05 / sqrt(D)) );
  g = dg_open(n);
  ng = dg_open(n);
  sg = dg_open(n - 1);
  mkgraph(g, x);
  die_if (!dg_biconnected(g), "initial diagram not biconnected D %d\n", D);
  fb = lookup ? dg_hsfb(g) : dg_hsfbrjw(g);
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
    sys0 = sys;

    i = (int) (rnd0() * n);
    rvn_rnddisp(xi, x[i], amp[sys]);
    dg_copy(ng, g);
    for (j = 0; j < n; j++) {
      if (j == i) continue;
      if (rvn_dist2(xi, x[j]) >= 1)
        dg_unlink(ng, i, j);
      else
        dg_link(ng, i, j);
    }

    if ( lookup ) { /* cheap lookup verseion */
      unqid_t gmapid = dg_getmapid(ng);
      nbc = dg_biconnected_lookuplow(n, gmapid);
      /* the weight `fb' is computed from the lookup table
       * for a small n, so it can be updated in every step */
      nfb = nbc ? dg_hsfb_lookuplow(n, gmapid) : 0;
      nnz = (nfb != 0);
      hasnfb = hasfb = 1;
    } else { /* direct and expensive computation, large n */
      nbc = dg_biconnected(ng);
      /* since computing fb is expensive,
       * we only try to find if a clique separator,
       * if so, fb == 0, otherwise fb is likely nonzero
       * so this is a good approximation of the ensemble */
      nnz = nbc ? (dg_cliquesep(ng) == 0) : 0;
      /* we have to compute `nfb' if `sys' == 3 */
      hasnfb = hasfb = (sys > 1);
      nfb = hasnfb ? dg_hsfbrjw(ng) : 0;
    }

#define WTSYS2(wt, fb) if ((wt = abs(fb)) == 0) wt = 1;
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
      if ( lookup || hasnfb ) {
        fb = nfb;
        hasfb = 1;
      }
    }

    /* compute fb only when necessary */
    ifb = (nstfb == 1 || rnd0() < 1./nstfb);
    if ( !hasfb && ifb ) {
      fb = dg_hsfbrjw(g);
      hasfb = 1;
    }

    /* change `sys' */
    if (sys == 0) { /* `sys' 0 */
      if ( nz ) { /* only a `nz' configuration can switch */
        if (rnd0() < 0.5) { /* switch `sys' from 0 to 1 */
          sys = 1;
        } else { /* switch `sys' from 0 to 2 */
          if ( ifb ) {
            WTSYS2(wt, fb);
            if ((r = wt*Z[0]/Z[2]) >= 1 || rnd0() < r)
              sys = 2;
          }
        }
      }
    } else if (sys == 1) { /* `sys' 1 */
      if (rnd0() < 0.5) { /* switch `sys' from 1 to 0 */
        if ((r = Z[1]/Z[0]) >= 1 || rnd0() < r)
          sys = 0;
      } else { /* switch `sys' from 1 to 2 */
        if ( ifb ) {
          tacc[0] += 1;
          WTSYS2(wt, fb);
          if ((r = wt*Z[1]/Z[2]) >= 1 || rnd0() < r) {
            sys = 2;
            tacc[1] += 1;
          }
        }
      }
    } else { /* `sys' 2 */
      if (rnd0() < 0.5) { /* switch `sys' from 2 to 1 */
        tacc[2] += 1;
        WTSYS2(wt, fb);
        if ((r = Z[2]/(Z[1]*wt)) >= 1 || rnd0() < r) {
          sys = 1;
          tacc[3] += 1;
        }
      } else { /* switch `sys' from 2 to 0 */
        WTSYS2(wt, fb);
        if ((r = Z[2]/(Z[0]*wt)) >= 1 || rnd0() < r)
          sys = 0;
      }
    }

    if (eql) {
      if (t >= nequil) {
        printf("equilibrated at t %g, nstfb %d, lookup %d\n", t, nstfb, lookup);
        t = 0;
        it = 0;
        eql = 0;
      }
      continue;
    }

    cacc[sys0] += acc1;
    ctot[sys0] += 1;

    /* accumulate data after equilibration */
    if (lookup || ifb || sys0 >= 2) {
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

      if (sys == 0) {
        if (rnd0() < 0.5) {
          /* the rate of biconnectivity after adding a random vertex */
          nrat[0] += 1;
          nrat[1] += rvn_voladd(x, n, xi);
        } else {
          /* the rate of biconnectivity after removing a random vertex */
          nrat[2] += 1;
          nrat[3] += dgmc_nremove(g, sg, n, &i);
        }
      }
    }

#define REPORT() { \
      fbav = (fbsum[0] + fbsum[1]) / (1 + hist[1]/fbnz) / hist[0]; \
      fprintf(stderr, "D %d, n %d, t%11g, cacc %.4f/%.4f/%.4f, " \
          "fb %+.8f/%+.8f|%+.8f %+.8f, %d\n" \
          "    nz %.7f, his1/0 %.7f, 2/1 %.7f; " \
          "tacc %.7f %.7f; n+: %.7f, n-: %.7f\n", \
        D, n, t, cacc[0]/ctot[0], cacc[1]/ctot[1], cacc[2]/ctot[2], \
        fbsum[0]/hist[0], fbsum[1]/hist[1] * (hist[1]/hist[0]) * Z[1]/Z[0], fbav, \
        fbsum[2]/hist[2] * (tacc[1]/tacc[0])/(tacc[3]/tacc[2]) * (Z[2]/Z[1]) * (fbnz/hist[0]), fb, \
        fbnz/hist[0], hist[1]/hist[0], hist[2]/hist[1], \
        tacc[1]/tacc[0], tacc[3]/tacc[2], nrat[1]/nrat[0], nrat[3]/nrat[2]); }

    if (it % nstrep == 0) {
      it = 0;
      REPORT();
    }
  }
  free(x);
  dg_close(g);
  dg_close(ng);
  dg_close(sg);

  REPORT();
  return fbav;
}



int main(int argc, char **argv)
{
  double fbav1, fbav2, nrat1[4], nrat2[4], tr1[4], tr2[4], nr1, nr2, rvir;

  doargs(argc, argv);
  printf("D %d, n %d, nsteps %g/%g, amp %g/%g, nstfb %d, code %d-bit\n",
      D, n, (double) nsteps, nsteps2, mcamp[0], mcamp[1],
      nstfb, (int) sizeof(code_t) * 8);

  fbav1 = mcrun(n, nequil, nsteps, mcamp,
      nstfb, lookup, nstnmv, nrat1, tr1);
  fbav2 = mcrun(n - 1, nequil, nsteps2, mcamp,
      nstfb, lookup, nstnmv, nrat2, tr2);
  nr1 = nrat1[3]/nrat1[2];
  nr2 = nrat2[1]/nrat2[0];
  rvir = fbav1/fbav2 * nr2/nr1 * (1. - n)/(2. - n)/n;
  printf("D %d, n %d, fbav %g/%g, nrat %g/%g = %g, vir%d/vir%d %g, %g, %g\n",
      D, n, fbav1, fbav2, nr2, nr1, nr2/nr1,
      n, n - 1, rvir, rvir * 2, rvir * ndvol(D));

  mtsave(NULL);
  return 0;
}

