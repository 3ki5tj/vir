/* count the number of connected and biconnected graphs */
#define ZCOM_PICK
#define ZCOM_UTIL
#include "zcom.h"
#include "dg.h"
#include "dgdb.h"

int n = 5; /* order */

typedef unsigned long long ulong_t;
typedef long long long_t;


/* pre-filter nonunique diagrams */
INLINE int filter(dg_t *d)
{
  int i, i0, deg, n = d->n, once = 0;

  /* we demand the vertices are sorted by the degrees */
  for (i0 = 0; i0 < n - 1; i0++) {
    /* 1. vertex i0 must have the highest degree */
    deg = dg_deg(d, i0);
    for (i = i0 + 1; i < n; i++) {
      if (dg_deg(d, i) > deg) return 1;
    }
  }
  return 0;
}


/* count degeneracy */
static void run(int n)
{
  int i, j;
  long_t ipr, npr, ngr, gid, even = 0, odd = 0, unique = 0, conn = 0, bconn = 0;
  dg_t *d;
  dgdb_t *db;

  npr = n * (n - 1) / 2;
  die_if (npr >= sizeof(ulong_t) * 8, "too many diagrams, n %d, npr %d\n", n, npr);
  ngr = 1 << npr;


  d = dg_open(n);
  db = dgdb_open(n);
  printf("# of pairs %llu, enumerating %llu\n", npr, ngr);
  conn = 0;
  for (gid = 0; gid < ngr; gid++) {
    /* construct a graph according to gid */
    dg_empty(d);
    for (ipr = 0, i = 0; i < n - 1; i++)
      for (j = i + 1; j < n; j++) {
        if ( (gid >> ipr) & 0x1u )
          dg_link(d, i, j);
        ipr++;
      }

    if ( dg_connected(d) ) conn++;
    if ( !dg_biconnected(d) ) continue;
    if (dg_nedges(d) % 2 == 0)
      even++;
    else
      odd++;
#if 0
    if (filter(d)) continue;
    //if (n <= 5) dg_print(d);
    unique++;
    if (dgdb_add(db, d) == 0) continue;
#endif
  }
  bconn = odd + even;
  printf("conn %lld (%g), biconn %lld (%g)\n", conn, 1.*conn/ngr, bconn, 1.*bconn/ngr);
  printf("unique %lld/%d, even %lld, odd %lld, tot %lld, diff %lld, ratio %g\n",
      unique, db->tot, even, odd, bconn, even - odd, 1.*even/bconn);
  dgdb_close(db);
  dg_close(d);
}



int main(int argc, char **argv)
{
  if (argc > 1) n = atoi(argv[1]);
  run(n);
  return 0;
}
