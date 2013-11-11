/* old mcutil.h routines */


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
  dgword_t vs = MKBITSMASK(n);
  int j;

  *i = (int) (rnd0() * n);
  vs ^= MKBIT(*i);
  for (j = 0; j < n; j++)
    /* see if removing i and j leaves the diagram connected */
    if ( j != (*i) && !dg_connectedvs(g, vs ^ MKBIT(j)))
      return 0;
  return 1;
}



/* return if removing vertex i leaves the diagram biconnected */
INLINE int dgmc_nremove1(const dg_t *g, dg_t *sg, int n, int *i)
{
  *i = (int) (rnd0() * n);
  dg_remove1(sg, g, *i);
  return dg_biconnected(sg);
}
#endif



/* recommended trial volume for docked-vertex moves between n and n + 1 */
INLINE double trialvol(int n, int d)
{
  double vol = 1;
  static const double ab[][2] = {{0.7, 1.8},
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
  dgword_t vs = MKBITSMASK(n);
  int i, j, k;

  j = (int) (rnd0() * n * (n - 1));
  i = j % n;
  k = j / n; /* choose a dock */
  if (k >= i) k++;
  if (rvn_dist2(x[i], x[k]) > rc * rc)
    return 0;
  vs ^= MKBIT(i);
  for (j = 0; j < n; j++)
    /* see if removing i and j leaves the diagram connected */
    if ( j != i && !dg_connectedvs(g, vs ^ MKBIT(j) ) )
      return 0;
  return 1;
}



