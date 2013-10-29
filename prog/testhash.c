#ifndef N
#define N 10
#endif

#include "dghash.h"
#include "testutil.h"


#ifndef HASHBITS
#define HASHBITS 0
#endif

#ifndef BLKSZ
#define BLKSZ 0
#endif

#ifndef MEMMAX
#define MEMMAX 0x40000000
#endif


/* this function simulates dghash_fbnr_lookup0()
 * it is used to test the coverage of the hash function */
INLINE double dghash_dummy_lookup(const dg_t *g)
{
  static dghash_t *hash;
  static code_t c[DGHASH_CWORDS];
#pragma omp threadprivate(hash, c);
  DG_DEFN_(g);
  int pos, ipos;
  dgls_t *ls;

  if (hash == NULL) {
    hash = dghash_open(DG_N_, HASHBITS, BLKSZ, MEMMAX);
    hash->dostat = 1; /* turn on statistics */
  }

  ls = hash->ls + dghash_getid(g, c, hash->cwords, hash->bits);
  pos = dgls_find(ls, c, hash->cwords, &ipos);
  hash->tot += 1;
  hash->hits += (pos >= 0);
  if (fmod(hash->tot, 1000000) < 0.1)
    dghash_printstat(hash, stderr);
  if (pos >= 0) { /* entry exists */
    return ls->arr[pos].fb;
  }
  dgls_add(ls, ipos, c, hash->cwords, 1.0, 1.0,
      hash->blksz, &hash->mem, hash->memmax);
  return 1;
}



/* test coverage randomly */
static void test_rndcover(int n, int nsamp, int nedmax)
{
  int t, ned, eql = 1, nequil = 1000, isamp = 0, good = 0, tot = 0;
  dg_t *g, *ng;
  double sum = 0, tsum = 0, rnp = 0.1;
  clock_t t0;

  printf("test coverage for n %d, nsamp %d, nedmax %d\n",
      n, nsamp, nedmax);
  g = dg_open(n);
  ng = dg_open(n);
  dg_full(g);
  ned = dg_nedges(g);

  for (t = 1; ; t++) {
    dg_rndswitchedge(g, &ned, nedmax, rnp);
    if (eql && t >= nequil) {
      t = 0;
      eql = 0;
      continue;
    }
    if (t % 5 != 0) continue;
    adjustrnp(ned, nedmax, t, 1000000, &good, &tot, &rnp);
    if ( ned > nedmax ) continue;
    //if ( dg_cliquesep(g) ) continue;

    t0 = clock();
    sum += dghash_dummy_lookup(g);
    tsum += clock() - t0;
    if (++isamp >= nsamp) break;
  }
  tsum /= CLOCKS_PER_SEC;
  printf("canonical label, n %d, samples %d, steps %d, nedmax %d; "
      "sum %g, time used: %gs/%d = %gmcs\n",
      n, nsamp, t, nedmax, sum / nsamp,
      tsum, nsamp, tsum / nsamp * 1e6);
  dg_close(g);
  dg_close(ng);
}




int main(int argc, char **argv)
{
  int n = 9, nsamp = 40000000, nedmax = 1000000;

#ifdef N
  n = N;
#endif
  if (argc >= 2) n = atoi(argv[1]);
  if (argc >= 3) nsamp = atoi(argv[2]);
  if (argc >= 4) nedmax = atoi(argv[3]);
  test_rndcover(n, nsamp, nedmax);
  return 0;
}
