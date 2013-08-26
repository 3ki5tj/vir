#ifndef MCUTIL_H__
#define MCUTIL_H__
/* utilities for MC sampling */



#define ZCOM_PICK
#define ZCOM_RVN
#include "zcom.h"
#include "dg.h"



/* compute the Ree-Hoover diagram from the coordinates */
INLINE void mkgraph(dg_t *g, rvn_t *x, int n)
{
  int i, j;

  g->n = n;
  dg_empty(g);
  for (i = 0; i < n - 1; i++)
    for (j = i + 1; j < n; j++)
      if (rvn_dist2(x[i], x[j]) < 1)
        DG_LINK(g, i, j);
}



/* remove the center of mass */
INLINE void rvn_rmcom(rvn_t *x, int n)
{
  int i;
  rvn_t xc;

  rvn_zero(xc);
  for (i = 0; i < n; i++) rvn_inc(xc, x[i]);
  rvn_smul(xc, (real) 1./n);
  for (i = 0; i < n; i++) rvn_dec(x[i], xc);
}



/* volume of d-dimensional sphere: vol(d) = pi^(d/2) / (d/2)!
 * we use the recursion: vol(d) = (2 * pi / d) vol(d - 2)
 * with vol(0) = 1, vol(1) = 2 */
INLINE double ndvol(int d)
{
  double vol;

  for (vol = (d % 2) ? 2 : 1; d > 1; d -= 2)
    vol *= 2 * M_PI / d;
  return vol;
}



#define SQRT2OVERPI   0.4501581580785531
#define ACOS13OVERPI  0.3918265520306073
#define SQRT3OVERPI   0.5513288954217921
#define ONEOVERPISQR  0.10132118364233778


/* return B3/B2^2 for d-dimensional system
 * M. Luban and A. Baram, J. Chem. Phys. 76. 3233 (1981) */
INLINE double b3rat(int d)
{
  int i, n;
  double fac, sum;

  if (d % 2 != 0) { /* odd, d = 2 n - 1*/
    /* B3/B2^2 = 2 (1 - 2F1(1/2, (1-d)/2, 3/2, 1/4) / B(1/2, (1+d)/2) ) */
    n = (d + 1) / 2;
    for (fac = 2, i = 1; i < n; i++)
      fac *= 2.*i / (2.*i + 1);
    for (sum = 1, i = n - 1; i > 0; i--)
      sum = 1 + sum * .25 * (i - n) * (2. * i - 1) / i / (2 * i + 1);
    return 2 - 2 * sum / fac;
  } else { /* even, d = 2 n */
    /* B3/B2^2 = 4/3 - n!/Sqrt(Pi)/Gamma(n+1/2) (3/4)^(n-1/2) 2F1(1, n+1, 3/2, 1/4)_n
     * = 4/3 - n!/(2*n-1)!! sqrt(3)/pi (3/2)^(n-1)
     *         sum_{i 0 to n - 1} (n+1)...(n+i)/(3*5...(2*i+1)) (3/2)^i */
    n = d / 2;
    for (fac = SQRT3OVERPI, i = 2; i <= n; i++)
      fac *= 1.5 * i / (2*i - 1);
    for (sum = 1, i = n - 1; i > 0; i--)
      sum = 1 + sum * 0.5 * (n + i) / (1 + 2 * i);
    return 4./3 - fac * sum;
  }
}



/* the partition function of 3-vertices biconnected configurations
 * divided by the square of the unit spherial volume */
INLINE double Z3rat(int d)
{
  return 3./4 * b3rat(d);
}



INLINE double B4rat_ring(int d)
{
  int n, i;
  double x, y, fac, sum;

  if (d % 2 == 1) {
    /* M. Luban and A. Baram, J. Chem. Phys. 76, 3233 (1982) */
    n = (d + 1) / 2;
    fac = 16./3;
    for (i = 2; i <= n; i++) {
      x = 2*i - 1;
      fac *= 8./3 * x * x * x * (i - 1) / ((6.*i - 5) * (6.*i - 7) * i * i);
    }
    for (sum = 1, i = n - 1; i > 0; i--)
      sum = 1 + sum * (i - 0.5) * (i - n) / (i + n) / (i + n);
    return fac * sum;
  } else {
    /* C. G. Joslin, J. Chem. Phys. 77, 2701 (1982) */
    n = d / 2;
    /* y = (2*n)!! / (2*n - 1)!!; */
    for (y = 1, i = 1; i <= n; i++)
      y *= i/(i - .5);
    for (x = 1, sum = 0, i = 1; i <= n; i++) {
      if (i > 1) x *= (i - 1.) / (i - .5);
      y *= 4. * (i * (n + i - 1)) / ((n + 2 * i - 1) * (n + 2 * i));
      sum += x * (4 + y) / i;
    }
    return 8 - sum * 8 * ONEOVERPISQR;
  }
  return 0;
}



INLINE double B4rat_diamond(int d)
{
  /* Clisby and McCoy J. Stat. Phys. */
  static double tab[] = {-8,
    -14./3,
    -8 + 8*SQRT3OVERPI + 20./3*ONEOVERPISQR,
    -6347./3360,
    -8 + 12*SQRT3OVERPI + 173./135*ONEOVERPISQR,
    -20830913./24600576.,
    -8 + 72./5*SQRT3OVERPI - 193229./37800.*ONEOVERPISQR,
    -87059799799./217947045888.,
    -8 + 558./35*SQRT3OVERPI - 76667881./7276500.*ONEOVERPISQR,
    -332647803264707./1711029608251392.,
    -8 + 594./35*SQRT3OVERPI - 9653909./654885.*ONEOVERPISQR,
    -865035021570458459./8949618140032008192.,
    -8 + 972./55*SQRT3OVERPI - 182221984415./10188962784.*ONEOVERPISQR
  };
  if (d <= 12) return tab[d];
  die_if (d > 12, "d %d is not supported\n", d);
  return 0;
}



/* fully connected Mayer diagram
 * approximate for n > 12 */
INLINE double B4rat_full(int d)
{
  /* Even: Clisby and McCoy, J. Stat. Phys. 114 1343 (2004)
   * Odd: Lyberg J. Stat. Phys. 119 747 (2005) */
  static double tab[] = {0, 0,
    8 - 12*SQRT3OVERPI + 8 * ONEOVERPISQR,
    (-89./280 - 219./1120*SQRT2OVERPI + 4131./2240*ACOS13OVERPI) * 4,
    8 - 18*SQRT3OVERPI + 238./9 * ONEOVERPISQR,
    (-163547./128128 - 3888425./8200192. * SQRT2OVERPI
            + 67183425./16400384.*ACOS13OVERPI) * 4,
    8 - 108./5*SQRT3OVERPI + 37259./900 * ONEOVERPISQR,
    (-283003297./141892608. - 159966456685./217947045888.*SQRT2OVERPI
            + 292926667005./48432676864.*ACOS13OVERPI) * 4,
    8 - 837./35*SQRT3OVERPI + 5765723./110250 * ONEOVERPISQR,
    (-88041062201./34810986496. - 2698457589952103./2851716013752320.*SQRT2OVERPI
            + 8656066770083523./1140686405500928.*ACOS13OVERPI) * 4,
    8 - 891./35*SQRT3OVERPI + 41696314./694575 * ONEOVERPISQR,
    (-66555106087399./22760055898112. - 16554115383300832799./14916030233386680320.*SQRT2OVERPI
            + 52251492946866520923./5966412093354672128.*ACOS13OVERPI) * 4,
    8 - 1458./55*SQRT3OVERPI + 88060381669./1344697200. * ONEOVERPISQR
  };
  if (d <= 12) return tab[d];
  /* use the approximate results */
  return 8 * 0.27433 * pow(0.66658, d - 2) * sqrt(3 / (d + 1));
}



INLINE double B4rat(int d)
{
  return -3*B4rat_ring(d)/8 - 3*B4rat_diamond(d)/4 - B4rat_full(d)/8;
}



/* the partition function of 4-vertices biconnected configurations
 * divided by the cube of the unit spherial volume */
INLINE double Z4rat(int d)
{
  return (3 * B4rat_ring(d) - 2 * B4rat_full(d)) / 8;
}



/* load the partition function of biconnected configurations
 * of the d-dimensional hard-sphere system from file */
INLINE int loadZrat(int d, int n, double *Z, const char *fn)
{
  char fndef[32], s[512];
  FILE *fp;
  int i, dim, nmax, n1;
  double Zr, x;

  if (fn == NULL) {
    sprintf(fndef, "ZrD%d.dat", d);
    fn = fndef;
  }
  xfopen(fp, fn, "r", return -1);
  /* handle the information line */
  if (fgets(s, sizeof s, fp) == NULL) {
    fprintf(stderr, "%s, no tag line\n%s", fn, s);
    fclose(fp);
    return -1;
  }
  if (s[0] == '#') {
    if (sscanf(s + 1, "%d %d", &dim, &nmax) != 2 || dim != d || nmax < n) {
      fprintf(stderr, "%s, d %d vs %d, nmax %d vs %d\n%s", fn, dim, d, nmax, n, s);
      fclose(fp);
      return -1;
    }
  } else { /* no info. line, give back the first line */
    rewind(fp);
  }

  for (i = 1; i < n; i++) {
    if (fgets(s, sizeof s, fp) == NULL)
      break;
    if (3 != sscanf(s, "%d%lf%lf", &n1, &Zr, &x) || i != n1) {
      fprintf(stderr, "%s ends on line %d\n%s", fn, i, s);
      break;
    }
  }
  *Z = x;
  fclose(fp);
  return 0;
}



/* return the sum of biconnected configurations for the d-dimensional
 * hard-sphere system */
INLINE double getZrat(int d, int n, const char *fn)
{
  double Z = 0;

  if (n <= 2) return 1;
  else if (n == 3) return Z3rat(d);
  else if (n == 4) {
    Z = Z4rat(d);
    /* the result is exact for d <= 12 */
    if (d <= 12) return Z;
  }
  /* load it from file */
  loadZrat(d, n, &Z, fn);
  return Z;
}



/* initialize a fully-connected configuration */
static void initx(rvn_t *x, int n)
{
  int i;
  real a = (real) (0.5 / sqrt(D));

  for (i = 0; i < n; i++)
    rvn_rnd(x[i], -a, a);
}



/* randomly displace a vertex */
#define DISPRNDI(i, nv, x, xi, amp, gauss) { \
  i = (int) (rnd0() * (nv)); \
  if (gauss) rvn_granddisp(xi, (x)[i], amp); \
  else rvn_rnddisp(xi, (x)[i], amp); }



/* construct a new graph `ng', with vertex i displaced */
#define UPDGRAPH(i, nv, g, ng, x, xi) { int j; \
  ng->n = (nv); \
  dg_copy(ng, g); \
  for (j = 0; j < (nv); j++) { \
    if (j == i) continue; \
    if (rvn_dist2(xi, (x)[j]) < 1) { DG_LINK(ng, i, j); } \
    else { DG_UNLINK(ng, i, j); } \
  } }



/* A Monte Carlo step of sampling biconnected configurations
 * This is a macro for the function-version is slower */
#define BCSTEP(acc, i, nv, g, ng, x, xi, amp, gauss) { \
  DISPRNDI(i, nv, x, xi, amp, gauss); \
  UPDGRAPH(i, nv, g, ng, x, xi); \
  if ( dg_biconnected(ng) ) { \
    rvn_copy((x)[i], xi); \
    dg_copy(g, ng); \
    acc = 1; \
  } else { acc = 0; } }



/* get sig2 that optimizes the acceptance rates */
INLINE double getgsig2(int n)
{
#if D == 2
  return 0.03 + 0.03*n;
#elif D == 3
  return 0.03 + 0.02*n;
#elif D == 8
  return 0.024 + 0.006*n;
#else
  return 0.024 + 0.048*n/D;
#endif
}



/* compute the logarithmic weight of a configuration `x'
 * without the center of mass */
INLINE double getloggwt(rvn_t *x, int n)
{
  int i;
  double d2sm, sig2;
  rvn_t xc;

  rvn_zero(xc);
  for (i = 0; i < n; i++)
    rvn_inc(xc, x[i]);
  rvn_smul(xc, 1./n);
  sig2 = getgsig2(n);
  for (d2sm = 0, i = 0; i < n; i++)
    d2sm += rvn_dist2(x[i], xc);
  return -0.5 * d2sm / sig2;
}



/* generate a configuration */
INLINE double gengconf(rvn_t *x, int n)
{
  double d2sm = 0, sig, sig2;
  int i;
  rvn_t xc;

  sig2 = getgsig2(n);
  sig = sqrt(sig2);
  rvn_zero(xc);
  for (i = 0; i < n; i++) {
    rvn_grand(x[i], 0, sig);
    d2sm += rvn_sqr(x[i]);
    rvn_inc(xc, x[i]);
  }
  d2sm -= rvn_sqr(xc) / n;
  return -0.5 * d2sm / sig2;
}



/* randomly replace a configuration */
INLINE int grepl(rvn_t *x, rvn_t *nx, dg_t *g, dg_t *ng)
{
  double logwt, lognwt;
  int n = g->n;

  lognwt = gengconf(nx, n);
  mkgraph(ng, nx, n);
  if ( !dg_biconnected(ng) ) return 0;
  logwt = getloggwt(x, n);
  if (logwt > lognwt || rnd0() < exp(logwt - lognwt)) {
    dg_copy(g, ng);
    rvn_ncopy(x, nx, n);
    return 1;
  } else return 0;
}


#if 0
/* compute the relative probability of adding a vertex
 * that preserves the biconnectivity of the graph
 * return the trial volume if successful, or 0 otherwise */
INLINE double rvn_voladd(rvn_t *x, int n, rvn_t xi)
{
  int j, deg = 0;
  real rad, r2, r2m = 0, r2n = 0, vol;

  /* naively, in a trial, the vertex should be added uniformly to
   * a very large volume Vmax, and the probability is computed as
   *  r = (V_bi / Vmax) * (Vmax / Vunit)             (1)
   * where Vmax/Vunit is the normalization to a uniform sphere
   * and B2 = Vunit / 2; (V_bi / Vmax) is the acceptance ratio.
   * We can however the shrink the sampling volume to a sphere
   * with the radius being the second largest distance from any
   * vertex to the center of mass (because there should be at
   * least two vertices connected to the new vertex to make the
   * graph biconnected). Thus
   *  r = (V_bi / V_samp) * (Vsamp / Vunit)         (2)
   * and (V_bi / V_samp) is the new and enlarged acceptance
   * probability. */
  for (j = 0; j < n; j++)
    if ((r2 = rvn_sqr(x[j])) > r2m) {
      r2n = r2m; /* second largest volume */
      r2m = r2;
    } else if (r2 > r2n) {
      r2n = r2;
    }
  if (r2n <= 0) r2n = r2m;

  /* we only need to sample the second largest sphere */
  rad = (real) (sqrt(r2n) + 1);
  /* compute the relative volume to the unit sphere */
  for (vol = 1, j = 0; j < D; j++) vol *= rad;
  /* see if adding a vertex leaves the graph biconnected */
  rvn_rndball(xi, rad);
  /* biconnectivity means xi is connected two vertices */
  for (j = 0; j < n; j++)
    if (rvn_dist2(xi, x[j]) < 1)
      if (++deg >= 2)
        return vol; /* see Eq. (2) above */

  return 0;
}



/* return if removing vertex i leaves the diagram biconnected */
INLINE int dgmc_nremove(const dg_t *g, int n, int *i)
{
  code_t vs = ((code_t) 1u << n) - 1u;
  int j;

  *i = (int) (rnd0() * n);
  vs ^= (code_t) 1u << (*i);
  for (j = 0; j < n; j++)
    /* see if removing i and j leaves the diagram connected */
    if ( j != (*i) && !dg_connectedvs(g, vs ^ ((code_t) 1u << j)))
      return 0;
  return 1;
}



/* return if removing vertex i leaves the diagram biconnected */
INLINE int dgmc_nremove1(const dg_t *g, dg_t *sg, int n, int *i)
{
  *i = (int) (rnd0() * n);
  dg_shrink1(sg, g, *i);
  return dg_biconnected(sg);
}
#endif



/* recommended trial volume for docked-vertex moves between n and n + 1 */
INLINE double trialvol(int n, int d)
{
  double vol = 1;
  static double ab[][2] = {{0.7, 1.8},
    {0.70, 1.8}, {0.70, 1.8}, {0.82, 1.9}, {0.90, 2.0}, {0.94, 2.0},
    {0.95, 2.0}, {0.95, 1.8}, {0.95, 1.8}, {0.95, 2.6}, {0.95, 3.0}};
  if (d <= 10) vol = ab[d][0] * n - ab[d][1];
  /* TODO: the optimal trial vol should be the ratio of partition function
   * it is currently hard to obtain a good formula for d > 10, use unity */
  return (vol > 1) ? vol : 1;
}



/* compute the relative probability of adding a vertex (n + 1)
 * that preserves the biconnectivity of the graph
 * return the trial volume if successful, or 0 otherwise */
INLINE int rvn_voladd_dock(rvn_t *x, int n, rvn_t xi, real rc)
{
  int j, k, deg = 0;

  /* choose a dock */
  k = (int) (rnd0() * n);
  rvn_rndball(xi, rc);
  rvn_inc(xi, x[k]); /* attach to the dock */

  /* biconnectivity means xi is connected two vertices */
  for (j = 0; j < n; j++)
    if (rvn_dist2(xi, x[j]) < 1)
      if (++deg >= 2)
        return 1;
  return 0;
}



/* return if removing vertex i leaves the diagram biconnected */
INLINE int dgmc_nremove_dock(const dg_t *g, rvn_t *x, int n, real rc)
{
  code_t vs = ((code_t) 1u << n) - 1u;
  int i, j, k;

  j = (int) (rnd0() * n * (n - 1));
  i = j % n;
  k = j / n; /* choose a dock */
  if (k >= i) k++;
  if (rvn_dist2(x[i], x[k]) > rc * rc)
    return 0;
  vs ^= (code_t) 1u << i;
  for (j = 0; j < n; j++)
    /* see if removing i and j leaves the diagram connected */
    if ( j != i && !dg_connectedvs(g, vs ^ ((code_t) 1u << j) ) )
      return 0;
  return 1;
}



#endif

