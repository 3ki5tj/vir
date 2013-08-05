#ifndef MCUTIL_H__
#define MCUTIL_H__
/* utilities for MC sampling */



/* compute a Ree-Hoover diagram from the coordinates */
INLINE void mkgraph(dg_t *g, rvn_t *x)
{
  int i, j, n = g->n;

  dg_empty(g);
  for (i = 0; i < n - 1; i++)
    for (j = i + 1; j < n; j++)
      if (rvn_dist2(x[i], x[j]) < 1)
        dg_link(g, i, j);
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



/* volume of D-dimensional sphere: pi^(n/2) / (n/2)! */
INLINE double ndvol(int n)
{
  double vol;

  for (vol = (n % 2) ? 2 : 1; n > 1; n -= 2)
    vol *= 2 * M_PI / n;
  return vol;
}



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
  mkgraph(ng, nx);
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
INLINE int dgmc_nremove_dock(const dg_t *g, rvn_t *x, int n, int *i, real rc)
{
  code_t vs = ((code_t) 1u << n) - 1u;
  int j, k;

  *i = (int) (rnd0() * n);
  k = (int) (rnd0() * (n - 1)); /* choose a dock */
  if (k >= *i) k++;
  if (rvn_dist2(x[*i], x[k]) > rc * rc)
    return 0;
  vs ^= (code_t) 1u << (*i);
  for (j = 0; j < n; j++)
    /* see if removing i and j leaves the diagram connected */
    if ( j != (*i) && !dg_connectedvs(g, vs ^ ((code_t) 1u << j)))
      return 0;
  return 1;
}



#endif

