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
INLINE int dgmc_nremove(const dg_t *g, dg_t *sg, int n, int *i)
{
  *i = (int) (rnd0() * n);
  dg_shrink1(sg, g, *i);
  return dg_biconnected(sg);
}



#endif

