#ifndef N
#define N 10
#endif

#include "dghash.h"
#include "testutil.h"


#ifndef HASHBITS
#define HASHBITS 0
#endif

#ifndef BLKMEM
#define BLKMEM 0
#endif

#ifndef MEMMAX
#define MEMMAX 0x40000000
#endif


/* this function simulates dghash_fbnr_lookup0()
 * it is used to test the coverage of the hash function */
INLINE double dghash_dummy_lookup(const dg_t *g)
{
  static dghash_t *hash;
  static dgword_t c[DG_CWORDS];
  DG_DEFN_(g);
  static dgvs_t ng_c[DG_NMAX];
  static dg_t ng[1];
  int pos, level = -1, enumiso = 1;
  dgword_t hashid;
  dgls_t *ls, *ls1;

  if (hash == NULL) {
    hash = dghash_open(DG_N_, HASHBITS, BLKMEM, MEMMAX, level, enumiso, 0);
    hash->dostat = 1; /* turn on statistics */
  }

  ng->n = g->n;
  ng->c = ng_c;

  /* find the representative graph of according to the automorphism level */
  dg_repiso(ng, g, hash->level);
  dg_encode(ng, c);
  DGHASH_GETID(hashid, c, DGHASH_CWORDS_(hash->cwords), hash->bits);
  ls = hash->ls + hashid;
  ls1 = dgls_find(ls, c, hash->cwords, &pos);
  hash->tot += 1;
  hash->hits += (ls1 != NULL);
  if (fmod(hash->tot, 1000000) < 0.1)
    dghash_printstat(hash, stderr);
  if (ls1 != NULL) { /* entry exists */
    return ls1->arr[pos].fb;
  }
  dgls_add(ls, c, 1.0, 1.0, hash);
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

