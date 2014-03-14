/** look up table */
class DiagramMap {
  public static final int NMAX = 7;
  int n = 0;
  short nUnique = 0;
  short [] map;
  int [] first;
  static final int nfirst[] = {1,
    1, 2, 4, 11, 34, 156, 1044, 12346};

  DiagramMap(int nv) {
    init(nv);
  }

  int init(int nv) {
    if (nv < 1 || nv > NMAX) return -1;
    if (nUnique > 0) return 0; // initialized
    n = nv;
    int npr = n * (n - 1) / 2;
    int ngr = 1 << npr;
    map = new short [ngr];
    for (int c = 0; c < ngr; c++)
      map[c] = -1;
    first = new int [ nfirst[n] ];

    if (n == 1) {
      map[0] = 0;
      first[0] = 0;
      return 0;
    }

    long startTime = System.nanoTime();
    System.out.printf("n %d: initializing map\n", n);
    // compute all permutations
    getPerm(n);
    long [][] masks = new long [npm][npr];
    for (int ipm = 0; ipm < npm; ipm++) {
      for (int ipr = 0, i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++, ipr++) {
          masks[ipm][ipr] = Bits.makeBit(
              getPairIndex(pm[ipm*n + i], pm[ipm*n + j], n) );
        }
      }
    }
    pm = null;

    // loop over all diagrams
    short gid = 0;
    int c;
    for (c = 0; c < ngr; c++) {
      if (map[c] >= 0) continue;
      first[gid] = c;
      for (int ipm = 0; ipm < npm; ipm++) {
        long c1 = 0;
        for (int ipr = 0; ipr < npr; ipr++) {
          if (((c >> ipr) & 1) != 0)
            c1 |= masks[ipm][ipr];
        }
        if (map[(int) c1] < 0) {
          map[(int) c1] = gid;
        } else {
          if (map[(int) c1] != gid) {
            System.out.printf("map conflict at %d, %d vs %d\n",
                (int) c1, map[(int) c1], gid);
          }
        }
      }
      gid++;
    }
    masks = null;
    nUnique = gid;
    long endTime = System.nanoTime();
    System.out.printf("n %d: initialized map, time %fs\n",
        n, (endTime - startTime)*1e-9);
    return 0;
  }

  /** Get the pair index from 0 to n*(n - 1)/2 - 1 */
  int getPairIndex(int i, int j, int n)
  {
    if (i > j) { int i1 = i; i = j; j = i1; }
    return n*i - (i + 1)*(i + 2)/2 + j;
  }

  int npm, pm[];
  int [] st = new int [NMAX + 2];
  boolean [] used = new boolean [NMAX + 2];

  /** Compute all permutations of n */
  int getPerm(int n) {
    int npp, ipp, i, top;

    for (npm = 1, i = 2; i <= n; i++)
      npm *= i;
    npp = npm * n; /* each permutation needs n numbers */
    pm = new int[npp];
    for (i = 0; i < n; i++) {
      st[i] = -1;
      used[i] = false;
    }
    for (ipp = 0, top = 0; ; ) {
      if (top >= n) {
        for (i = 0; i < n; i++)
          pm[ipp++] = st[i];
      } else {
        for (i = st[top]; ++i < n; )
          if ( !used[i] ) break;
        if (i < n) {
          used[i] = true;
          st[top] = i;
          st[++top] = -1;
          continue;
        }
      }
      // exhausted this level
      if (--top < 0) break;
      used[ st[top] ] = false;
    }
    return npm;
  }
}



