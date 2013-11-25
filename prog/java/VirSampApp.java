/** Computing the virial coefficients of the 2D hard sphere fluid
 *  by Monte Carlo simulation */
import java.applet.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.*;
import static java.lang.Math.*;
import java.text.*;
import javax.swing.*;
import java.util.Hashtable;
import java.util.Random;



class Diagram {
  public static final int NMAX = 63; // maximal number of vertices
  int n;
  long maskn;
  long c[];
  long RHSC_c[]; // temporary c used by RHSC iteration
  int degs[];

  /** Construct a graph */
  Diagram(int nv) {
    if (nv > NMAX) nv = NMAX;
    n = nv;
    maskn = Bits.makeBitsMask(n);
    c = new long [n];
    RHSC_c = new long [n];
    degs = new int [n];
  }

  /** Construct from D-dimensional coordinates */
  Diagram(int nv, int dim, double r[][]) {
    this(nv);
    fromCoordinates(dim, r);
  }

  /** Update the graph form D-dimensional coordinates */
  void fromCoordinates(int dim, double r[][]) {
    int i, j, k;
    double r2, dx;
    empty();
    for (i = 1; i < n; i++) {
      for (j = 0; j < i; j++) {
        r2 = 0;
        for (k = 0; k < dim; k++) {
          dx = r[i][k] - r[j][k];
          r2 += dx * dx;
        }
        if (r2 < 1) link(i, j);
      }
    }
  }

  /** Check if two vertices i and j are linked */
  boolean isLinkedc(long cc[], int i, int j) {
    return ((cc[i] >> j) & 1L) != 0;
  }

  /** Check if two vertices i and j are linked */
  boolean isLinked(int i, int j) { return isLinkedc(c, i, j); }

  /** Link two vertices i and j */
  void linkc(long cc[], int i, int j) {
    cc[i] |= 1L << j;
    cc[j] |= 1L << i;
  }

  /** Link two vertices i and j */
  void link(int i, int j) { linkc(c, i, j); }

  /** Unlink two vertices i and j */
  void unlinkc(long cc[], int i, int j) {
    cc[i] &= ~(1L << j);
    cc[j] &= ~(1L << i);
  }

  /** Unlink two vertices i and j */
  void unlink(int i, int j) { unlinkc(c, i, j); }

  /** Copy from src */
  void copy(Diagram src) {
    n = src.n;
    for (int i = 0; i < n; i ++)
      c[i] = src.c[i];
  }

  /** Remove all edges in the graph */
  void empty() {
    for (int i = 0; i < n; i ++)
      c[i] = 0;
  }

  /** Print the diagram */
  void printc(long cc[]) {
    String s = "";

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if ( isLinkedc(cc, i, j) )
          s += "* ";
        else
          s += "  ";
      }
      s += "\n";
    }
    System.out.println(s + "biconnected: " + biconnected());
  }

  void print() { printc(c); }

  /** Add all edges to the graph */
  void full() {
    for (int i = 0; i < n; i++) {
      c[i] = maskn ^ Bits.makeBit(i);
    }
  }

  /** Check if the subgraph of `vs' is connected */
  boolean connectedvsc(long cc[], long vs) {
    long stack = vs & (-vs);
    while (stack != 0) {
      long bk = stack & (-stack);
      int k = Bits.bit2id(bk);
      vs ^= bk; // remove k from the to-do list
      stack = (stack | cc[k]) & vs;
      if (stack == vs)
        return true;
    }
    return false;
  }

  /** Check if the subgraph of `vs' is connected */
  boolean connectedvs(long vs) {
    return connectedvsc(c, vs);
  }

  /** Check if the graph is connected */
  boolean connectedc(long cc[]) {
    return connectedvsc(cc, maskn);
  }

  /** Check if the graph is connected */
  boolean connected() {
    return connectedvsc(c, maskn);
  }

  /** Check if the graph is biconnected */
  boolean biconnectedc(int nn, long cc[]) {
    if (nn > 2) {
      long mask = Bits.makeBitsMask(nn);
      for (long b = 1; (b & mask) != 0; b <<= 1)
        if ( !connectedvsc(cc, mask ^ b) )
          return false;
      return true;
    } else if (nn == 2) {
      return ((cc[1] & 0x1l) != 0) ? true : false;
    } else return true;
  }

  /** Check if the graph is biconnected */
  boolean biconnected() {
    return biconnectedc(n, c);
  }

  /** Degree of vertex i */
  int getDegreec(long cc[], int i) {
    return Bits.bitCount(cc[i]);
  }

  /** Degrees of all vertices */
  int getDegreesc(long cc[], int degs[]) {
    int nedges = 0;
    for (int i = 0; i < n; i++) {
      degs[i] = Bits.bitCount(cc[i]);
      nedges += degs[i];
    }
    return nedges / 2;
  }

  /** Get the number of edges */
  int getNumEdgesc(long cc[]) {
    return getDegreesc(cc, degs);
  }

  /** Get the number of edges */
  int getNumEdges() { return getNumEdgesc(c); }

  /** ====================== Star content ====================== */

  /** Ree-Hoover formula to convert the star content of
   *  a smaller diagram to that of a larger diagram */
  double getRHSCIter(int n, int n0, double sc) {
    if (n < n0) return 0;
    for (int i = n0; i < n; i++) {
      sc *= (i - 1) * ((i % 2 != 0) ? -1 : 1);
    }
    return sc;
  }

  boolean errRHSCSpec0;

  /** Special Ree-Hoover star content */
  double getRHSCSpec0c(long cc[]) {
    errRHSCSpec0 = false;
    int nedges = getDegreesc(cc, degs);

    if (nedges == n) {
      return 1;
    } else if (nedges == n + 1) {
      int i, j;
      for (i = 0; i < n; i++)
        if (degs[i] == 3) break;
      for (j = i + 1; j < n; j++)
        if (degs[j] == 3) break;
      if ( j < n && !isLinkedc(cc, i, j) )
        return 1;
      else
        return 0;
    }

    int nedgesinv = n * (n - 1) / 2 - nedges;
    if (nedgesinv < 3) {
      if (nedgesinv == 0) { // fully connected
        return getRHSCIter(n, 2, 1);
      } else if (nedgesinv == 1) {
        return 0;
      } else {
        // if there are two wiggly lines, the SC is nonzero only if
        // the two wiggly lines are not connected to the same vertex
        for (int i = 0; i < n; i++)
          if (degs[i] <= n - 3) return 0;
        return getRHSCIter(n, 4, 1);
      }
    }

    errRHSCSpec0 = true; // not a special case
    return 0;
  }

  /** Recursively compute Ree-Hoover star content
   *  starting from the edge (i, j) */
  double getRHSCRecur(long cc[], int sgn, int i, int j) {
    long maskn = Bits.makeBitsMask(n), avsi, avs, bi, bj;
    double sc = 0;

    if (++j >= n) j = (++i) + 1;
    for (; i < n - 1; j = (++i) + 1) {
      if (Bits.bitCount( cc[i] ) <= 2) continue;
      bi = Bits.makeBit(i);
      avsi = maskn ^ bi;
      for (; j < n; j++) {
        // try to remove the edge i-j
        if ( (cc[j] & bi) == 0 || Bits.bitCount(cc[j]) <= 2 )
          continue;
        // unlink i and j
        cc[j] &= ~bi;
        bj = Bits.makeBit(j);
        cc[i] &= ~bj;
        avs = avsi ^ bj;
        if ( biconnectedc(n, cc) ) {
          sc += sgn + getRHSCRecur(cc, -sgn, i, j);
        }
        // link back i and j
        cc[i] |= bj;
        cc[j] |= bi;
      }
    }
    return sc;
  }

  /** Compute Ree-Hoover star content directly */
  double getRHSCDirectLowc(long cc[]) {
    for (int i = 0; i < n; i++)
      RHSC_c[i] = cc[i];
    return 1 + getRHSCRecur(RHSC_c, -1, 0, 0);
  }

  /** Compute Ree-Hoover star content directly */
  double getRHSCDirectc(long cc[]) {
    double sc = getRHSCSpec0c(cc);
    if (errRHSCSpec0)
      sc = getRHSCDirectLowc(cc);
    //System.out.println("sc " + sc + ", spec " + !errRHSCSpec0);
    return sc;
  }

  /** Compute Ree-Hoover star content directly */
  double getRHSCDirect(long cc[]) {
    return getRHSCDirectc(cc);
  }

  double getFbDirectc(long cc[]) {
    double sc = getRHSCDirectc(cc);
    return SC2FB(sc, getNumEdgesc(cc));
  }

  /** Convert star content to sum of clusters */
  double SC2FB(double sc, int ned) {
    return sc * (1 - 2 * (ned % 2));
  }

  /** Convert sum of clusters to star content */
  double FB2SC(double fb, int ned) {
    return fb * (1 - 2 * (ned % 2));
  }



  /** ==================  Clique Separator ================= */
  long [] csFillc = new long [NMAX];
  int [] csPerm = new int [NMAX];

  int [] csL2 = new int [NMAX]; // the label l times 2
  long [] csReach = new long [NMAX]; // csReach[label] gives a set of vertices
  int [] csCnt = new int [NMAX * 2];

  /** Compute a minimal order and the corresponding fill-in of a graph.
   *  The minimal order is an order of removing vertices such than the
   *  resulting edges by contraction is locally a minimum.
   *  the algorithm used here constructs a hierarchical spanning tree
   *  ``Algorithmic aspects of vertex elimination on graphs''
   *  Donald J. Rose, R. Endre Tarjan and George S. Lueker
   *  SIAM J. Comput. Vol 5, No. 2, June 1976 */
  void getMinimalOrderc(long c[]) {
    int i, j, k, v = 0, w, z, l;
    long numbered, reached, bv, bw, bz, r;

    for (i = 0; i < n; i++) csFillc[i] = c[i];
    for (i = 0; i < n; i++) csL2[i] = 0; // l(i) * 2
    reached = numbered = 0;
    k = 1; // number of different labels
    for (i = n - 1; i >= 0; i--) {
      // select the vertex with the largest label
      bw = numbered;
      for (r = ~numbered & maskn; r != 0; r ^= bv) {
        bv = r & (-r);
        v = Bits.bit2id(bv);
        if ( csL2[v]/2 == k - 1 ) { // found it
          numbered |= bv;
          break;
        }
      }
      csPerm[i] = v;
      // csReach[l] gives the set of vertcies with the same label `l'
      // the label `l' is the hierarchy level of the spanning tree
      // vertices with the same label are treated as the same
      for (l = 0; l < k; l++) csReach[l] = 0;
      // set all unnumbered vertices as unreached
      reached = numbered;

      // handle immediate neighbors of v
      for (r = c[v] & ~numbered & maskn; r != 0; r ^= bw) {
        bw = r & (-r);
        w = Bits.bit2id(bw);
        // the immediate neighbors of w are the starting points
        // of the search, we group them by the labels
        csReach[ csL2[w]/2 ] |= bw;
        reached |= bw;
        // only update the label after we have done with it
        // the csL2[w]++ operation only applies to a vertex once
        csL2[w]++;

        linkc(csFillc, v, w);
      }

      // search paths from v
      // the loop is over all labels, from smaller labels to larger ones.
      // during the search, only vertices of equal or larger labels are
      // produced, so the one-pass search is sufficient
      for (j = 0; j < k; j++) {
        while ( csReach[j] != 0 ) { // if the vertex set is not empty
          // delete the first w from csReach[j]
          bw = csReach[j] & (-csReach[j]);
          w = Bits.bit2id(bw);
          csReach[j] ^= bw; // remove w from csReach[j]
          while ( (r = (c[w] & ~reached)) != 0 ) {
            bz = r & (-r);
            z = Bits.bit2id(bz);
            reached |= bz;
            if ((l = csL2[z]/2) > j) {
              csReach[l] |= bz;
              csL2[z]++;
              linkc(csFillc, v, z);
            } else { // lower label encountered, count it as label j
              csReach[j] |= bz;
            }
          }
        }
      }

      // re-assign labels of the vertices by counting sort
      for (l = 0; l < 2*n; l++) csCnt[l] = 0;
      // accumulate the number of visits of each label
      for (r = ~numbered & maskn; r != 0; r ^= bw) {
        bw = r & (-r);
        w = Bits.bit2id(bw);
        csCnt[ csL2[w] ] = 1;
      }
      // count the number of different labels
      for (k = 0, l = 0; l < 2*n; l++)
        if ( csCnt[l] != 0 ) // compute the new label
          csCnt[l] = k++; // csCnt[l] is now the new label
      for (r = ~numbered & maskn; r != 0; r ^= bw) {
        bw = r & (-r);
        w = Bits.bit2id(bw);
        csL2[w] = 2 * csCnt[ csL2[w] ]; // set the new label
      }
    }
  }

  /** Decompose a diagram by clique separators.
   *  A clique separator is a fully-connected subgraph
   *  `g' is the input diagram, `f' is the fill-in diagram
   *  `a' is the elimination order, a[0] is the first vertex to eliminate
   *  return the number of cliques, `cl' is the array of cliques
   *  `stop1' means stop the search after the first clique separator
   *  The algorithm first find a minimal order of elimination.
   *  Using this order on a graph with a clique separator, at least
   *  one part of the graph is eliminated before the clique
   *  ``Decomposition by clique separators'' Robert E. Tarjan,
   *  Discrete Mathematics 55 (1985) 221-232 */
  long decompCliqueSepc(long c[]) {
    long unvisited = maskn;
    long cc, r, cb, bw;

    for (int i = 0; i < n; i++) {
      int v = csPerm[i];
      unvisited ^= Bits.makeBit(v); // remove the `v' bit
      // compute C(v), the set of succeeding vertices that
      // are adjacent to v
      cc = unvisited & csFillc[v];
      // test if C(v) is a clique, a fully-connected subgraph
      for (r = cc; r != 0; r ^= bw) {
        bw = r & (-r);
        int w = Bits.bit2id(bw);
        // c ^ bw is the set of vertices connected to `w'
        // in `c', if `c' is a clique
        cb = cc ^ bw;
        if ((c[w] & cb) != cb) // not a clique
          break; // break the loop prematurally, r != 0
      }
      if (r == 0) { // if the loop is completed, `cc' is a clique
        if (unvisited == cc) { // clique `cc' == the rest vertices
          return 0;            // so it is not a separator
        } else { // found a clique `cc'
          return cc;
        }
      }
    }
    return 0;
  }

  /** General algorithm for clique separators */
  long getCliqueSepc(long c[]) {
    getMinimalOrderc(c);
    return decompCliqueSepc(c);
  }

  /** Clique separator of two vertices */
  long getCSep2c(long cc[]) {
    // loop over the first vertex
    for (int i = 1; i < n; i++) {
      long bi = Bits.makeBit(i), b;
      long maski = maskn ^ bi;
      long cci = cc[i] & Bits.makeBitsMask(i);
      // loop over vertices connected to i
      // with indices lower than i
      for (; cci != 0; cci ^= b) {
        b = cci & (-cci);
        if ( !connectedvsc(cc, maski ^ b) ) {
          //printc(cc); System.out.printf("csep2: %d and %d\n", i, Bits.bit2id(b));
          return bi ^ b;
        }
      }
    }
    return 0;
  }

  /** Clique separator, best strategy */
  long getCSepc(long cc[]) {
    long cs;
    return ((cs = getCSep2c(cc)) != 0) ? cs : getCliqueSepc(cc);
  }

  /** Clique separator, best strategy */
  long getCSep() { return getCSepc(c); }



  /** ====================== Wheatley's method ====================== */
  public static final int RJW_NMAX = 22;
  public static final double RJW_FBDIRTY = -1e301;

  public static boolean isFbRJWInvalid(double x) {
    return (x < -1e300);
  }

  /** Increment ms1 in the limits of bits of ms1 | ms2 */
  long incBits(long ms1, long ms2) {
    long nms2, lbit, hbits;
    nms2 = -ms2;
    lbit = ms2 & nms2; // the lowest bits
    hbits = ms2 ^ nms2; // collect bits higher than lbit
    ms1 &= hbits; // wipe out bits lower or equal to the lowest bit
    ms1 |= lbit; // add the lowest bit
    return ms1;
  }

  /** Compute Boltzmann weight, the sum of all diagrams
   *  `cc' is the connectivity matrix, `vs' is the vertex set */
  int getFqRJWLow(long cc[], long vs) {
    long w, bi;
    // if there is a bond, i.e., r(i, j) < 1, then fq = 0
    for (w = vs; w != 0; w ^= bi) {
      // if cc[i] share vertices with vs, there is bond
      // the Boltzmann weight = \prod_(ij) e_ij, and it allows no clash
      // therefore return zero immediately
      bi = w & (-w);
      int i = Bits.bit2id(bi);
      if ((cc[i] & vs) != 0) return 0;
    }
    return 1;
  }

  /** Compute the sum of connected diagrams by Wheatley's recursion formula
   *  `cc' is the connectivity matrix,  `vs' is the vertex set */
  double getFcRJWLow(long cc[], long vs, double [] fcArr, double [] fqArr)
  {
    double fc, fc1, fq2;
    long ms1, ms2, b1;

    b1 = vs & (-vs);
    if ( (vs ^ b1) == 0 ) return 1; // only one vertex
    if ( isFbRJWInvalid(fc = fqArr[(int) vs]) )
      fqArr[(int) vs] = fc = getFqRJWLow(cc, vs); // start with fq
    vs ^= b1; // remove vertex b1 from vs, vs = vs - {b1}
    // before the above statement, `vs' is the set of all vertices
    // since `b1' is fixed, after excluding `b1' in the statement,
    // `vs' is the set of *variable* vertices
    // `ms1' is the subset of variable vertices
    // `ms2' is the complement set of variable vertices
    // loop over subsets ms1 of vs, stops when the complement set ms2 == 0
    for (ms1 = 0; (ms2 = ms1 ^ vs) != 0; ) {
      if ( isFbRJWInvalid(fq2 = fqArr[(int) ms2]) )
        fqArr[(int) ms2] = fq2 = getFqRJWLow(cc, ms2); // fq of the complement set
      if (fq2 != 0) {
        long vs1 = ms1 | b1; // vs1 = ms1 + {b1}
        if ( isFbRJWInvalid(fc1 = fcArr[(int) vs1]) )
          fcArr[(int) vs1] = fc1 = getFcRJWLow(cc, vs1, fcArr, fqArr); // recursion
        fc -= fc1 * fq2;
      }
      // update the subset `ms1'
      ms1 = incBits(ms1, ms2);
    }
    return fc;
  }

  /** Compute the sum of all connected diagrams without the articulation point
   *  at vertices lower than `v', `vs' is the vertex set
   *  if v = 0, it returns fc; if v = n, it returns fb */
  double getFbRJWLow(long [] cc, int v, long vs, double faArr[], double fbArr[])
  {
    int i;
    double fb, fa;
    long r, b, bv = Bits.makeBit(v);
    int id;

    if ((i = Bits.bitCount(vs)) <= 1) {
      return 1;
    } else if (i == 2) {
      return ((cc[ Bits.bitFirst(vs) ] & vs) != 0) ? -1 : 0;
    }

    // start with the sum of connected diagrams, the first 2^n numbers of
    // fbArr and faArr are used for saving fcArr and fqArr, respectively
    if ( isFbRJWInvalid(fb = fbArr[(int) vs]) ) {
      fb = getFcRJWLow(cc, vs, fbArr, faArr);
      fbArr[(int) vs] = fb;
    }
    // remove diagrams with the lowest articulation points at i < v
    for (r = vs & (bv - 1); r != 0; r ^= b) {
      b = r & (-r);
      i = Bits.bit2id(b);
      id = ((i + 1) << n)  + (int) vs; // (i + 1) * 2^n + vs
      if ( isFbRJWInvalid(fa = faArr[id]) ) {
        fa = getFaRJWLow(cc, i, vs, faArr, fbArr);
        faArr[id] = fa;
      }
      fbArr[id] = (fb -= fa);
    }
    return fb;
  }

  /** Compute the sum of all connected diagrams with the articulation point at v
   *  and no articulation point at any vertex lower than v */
  double getFaRJWLow(long [] cc, int v, long vs, double faArr[], double fbArr[])
  {
    double fa = 0, fb, fa2, fb2;
    long ms1, ms2, vs1, bv = Bits.makeBit(v), b1, b1v;
    int id0, id;

    b1 = vs & (-vs); // lowest vertex
    if ( b1 == bv ) { // if vertex 1 coincide with `v', find the next lowest
      b1 = vs ^ bv; // remove `bv' from the vertex set
      b1 = b1 & (-b1);
    }
    b1v = b1 ^ bv;
    vs ^= b1v; // remove the fixed vertices `b1' and `bv' from `vs'
    // `vs' is the set of *variable* vertices from now on
    if ( vs == 0 ) return 0; // no articulated diagram with only two vertices
    id0 = (v + 1) << n; // (v + 1) * 2^n
    // `id0' is the offset for vertex v
    // loop over subsets of vs, stops when vs == vs1
    for (ms1 = 0; (ms2 = (ms1 ^ vs)) != 0; ) {
      vs1 = ms1 | b1v; // add the two fixed vertices */
      id = id0 + (int) vs1;
      if ( isFbRJWInvalid(fb = fbArr[id]) ) { // compute fb if necessary
        fb = getFbRJWLow(cc, v + 1, vs1, faArr, fbArr);
        fbArr[id] = fb;
      }
      if ( fb != 0 ) {
        long vs2 = ms2 | bv; // unused variable vertices + the articulation point `bv'
                             // `|' is equivalent to `+' here
        id = id0 + (int) vs2;
        fb2 = fbArr[id];
        if ( isFbRJWInvalid(fb2 = fbArr[id]) ) {
          fb2 = getFbRJWLow(cc, v + 1, vs2, faArr, fbArr);
          fbArr[id] = fb2; /* save the fb value */
        }
        if ( isFbRJWInvalid(fa2 = faArr[id]) ) {
          fa2 = getFaRJWLow(cc, v, vs2, faArr, fbArr);
          faArr[id] = fa2; /* save the fa value */
        }
        fa += fb * (fb2 + fa2);
      }
      ms1 = incBits(ms1, ms2); // update the subset `ms1'
    }
    return fa;
  }

  int nmaxRJW = -1;
  double faRJWArr[], fbRJWArr[];

  /** Compute the sum of biconnected diagrams by Wheatley's method.
   *  This is a low level function and the test of clique separator
   *  is not done here. */
  double getFbRJWc(long c[]) {
    int size;

    /* the memory requirement is 2^(n + 1) * (n + 1) * sizeof(double) */
    if (n > nmaxRJW) {
      nmaxRJW = n;
      size = (nmaxRJW + 1) << nmaxRJW; /* (nmax + 1) * 2^nmax */
      faRJWArr = new double [size];
      fbRJWArr = new double [size];
    }
    size = ((n + 1) << n); /* (n + 1) * 2^n */
    for (int i = 0; i < size; i++) faRJWArr[i] = RJW_FBDIRTY;
    for (int i = 0; i < size; i++) fbRJWArr[i] = RJW_FBDIRTY;
    return getFbRJWLow(c, n, maskn, faRJWArr, fbRJWArr);
  }

  /** Compute the sum of biconnected diagrams by Wheatley's method. */
  double getFbRJW() { return getFbRJWc(c); }



  /** Compute the hard-sphere weight, using the best strategy */
  double getFbc(long c[]) {
    // detect clique separators first
    if ( getCSepc(c) != 0 ) return 0;
    int nedges = getNumEdgesc(c);
    if (nedges < 2*n - 2 || n > RJW_NMAX) {
      double sc = getRHSCDirectc(c);
      return SC2FB(sc, nedges);
    } else {
      return getFbRJWc(c);
    }
  }

  /** Compute the hard-sphere weight, using the best strategy */
  double getFb() { return getFbc(c); }



  /** ====================== Ring content ====================== */
  boolean errRingSpec0;

  /** Special ring content */
  double getRingSpec0(long cc[]) {
    errRingSpec0 = false;
    int nedges = getDegreesc(cc, degs);

    if (nedges == n) {
      return 1;
    } else if (nedges == n + 1) {
      int i, j;
      for (i = 0; i < n; i++)
        if (degs[i] == 3) break;
      for (j = i + 1; j < n; j++)
        if (degs[j] == 3) break;
      if ( j < n && isLinkedc(cc, i, j) )
        return 1;
      else
        return 0;
    }

    // fully-connected diagram
    if (n * (n - 1) / 2 == nedges) {
      double x = 1;
      for (int i = 3; i < n; i++) x *= i;
      return x;
    }

    errRingSpec0 = true; // cannot be handled by special cases
    return 0;
  }

  /** Compute the ring content directly */
  double getRingDirectLow(long c[]) {
    if (n <= 2) {
      if (n <= 1) return 1;
      else return (double) (c[1] & 0x1L);
    }

    int [] st = new int [n];
    int top, sttop, root = 0;
    long maskn = Bits.makeBitsMask(n);
    long unused, cc, b, croot, ccp;
    double cnt = 0;

    st[0] = root;
    st[1] = 0;
    top = 1;
    unused = maskn ^ Bits.makeBit(root);
    croot = c[root];
    ccp = unused & croot;
    sttop = st[top];

    while (true) {
      long ccp0 = ccp;
      cc = ccp & ~Bits.makeBitsMask(sttop + 1);
      if (cc != 0) {
        b = cc & (-cc);
        sttop = Bits.bit2id(b);
        unused ^= b;
        if ((ccp = unused & c[(int) sttop]) != 0) { // there are still unused vertices
          if (top < n - 2) {
            st[top++] = sttop;
            sttop = 0;
            continue;
          }
          if ((croot & ccp) != 0) cnt++;
        }
        ccp = ccp0;
        unused ^= b;
      } else {
        if (--top == 0) break;
        sttop = st[top];
        unused ^= Bits.makeBit(sttop);
        ccp = unused & c[(int) st[top - 1]];
      }
    }
    return cnt / 2.;
  }

  /** Compute the ring content */
  double getNrc(long cc[]) {
    double nr = getRingSpec0(cc);
    if ( errRingSpec0 )
      nr = getRingDirectLow(cc);
    return nr;
  }

  /** Compute the ring content */
  double getNr() { return getNrc(c); }

  public static final double [][] BringArr = new double [][] {
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
    {1, 1, 1, 0.78200443791154128382, -1.37886106172259, 2.3892313272351075694, -4.1982838514716088621, 7.4641229611560658631, -13.417394242700023459, 24.348397179977765504, -44.545169259816196107, 82.065204605387569236, -152.10101217798660649, 283.38605242980672054, -530.41494476983759846, 996.80180301564444021, -1880.0113187851478152, 3557.1654496456511541, -6749.9065867469811732, 12841.623492027556613, -24488.69260296270962, 46799.755884672915104, -89613.665306633004656, 171904.63565088797016, -330310.10004020167227, 635656.12620385592595, -1.2250107850950705632e6, 2.3639114083850431951e6, -4.5672717993355103991e6, 8.8344909846425886751e6, -1.7106975673159133898e7, 3.3159117037815562106e7, -6.4334552630282360775e7, 1.2493183627597122774e8, -2.4281047727332852294e8, 4.7228805020026995818e8, -9.1933359327104603302e8, 1.7908070779996959054e9, -3.4907451743635388126e9, 6.8087363102059654382e9, -1.3288633624365976061e10, 2.5950574065009191677e10, -5.0705457033020012524e10, 9.9127141482284093514e10, -1.9388776225287457383e11, 3.794185831455085721e11, -7.4282795429486973759e11, 1.4549596043419985311e12, -2.8510104890733708662e12, 5.5888760133346984222e12, -1.0960262307716214717e13, 2.1502139129810996513e13, -4.219883276465209213e13, 8.2845974054818577379e13, -1.6270055891629477687e14, 3.1963042357048416046e14, -6.2812147670603438279e14, 1.2347269785213119943e15, -2.4278705498979379314e15, 4.775327576882795128e15, -9.3950658164588087833e15, 1.8488925059895925343e16, -3.6394437989399663759e16, 7.1658284424048065756e16, -1.4112452872736966119e17},
    {1, 1, 1, 0.625, -0.97142857142856133707, 1.5234002976190477212, -2.4609257409257409257, 4.0715923664551594239, -6.8723786226059556444, 11.792633697111110532, -20.515380917802643582, 36.105301340712559083, -64.170978978549113475, 115.02318363499727409, -207.69571956858503784, 377.45838894939517474, -689.89349434611755891, 1267.337145083867775, -2338.6642738139117711, 4333.2333031038861879, -8058.5257699675396924, 15036.692010354224822, -28143.302723473549791, 52821.731637650413991, -99395.596775980278516, 187478.74455223572446, -354397.29120775361195, 671295.24828310594615, -1.2739725213605080975e6, 2.4220021366160856354e6, -4.6121880185054193584e6, 8.7965476295311151423e6, -1.680154337162424471e7, 3.2135147456502856256e7, -6.1541968648285512514e7, 1.1800246425093422074e8, -2.265220944047539378e8, 4.3531470950574894135e8, -8.3742309514171835524e8, 1.6125482738709026246e9, -3.1080340007628270426e9, 5.9957649114628493491e9, -1.1576336646658294223e10, 2.2369103275300677989e10, -4.3257474939004569946e10, 8.3713256757896221007e10, -1.6211905221101504277e11, 3.1417219436930284155e11, -6.0923272729202173639e11, 1.182140628124325393e12, -2.2951698007626344311e12, 4.4587199703728595457e12, -8.6665407957403010497e12, 1.6854374088804954812e13, -3.2794594237786493545e13, 6.3842028067308531087e13, -1.2434212425236683907e14, 2.4228687830475620211e14, -4.7231820712145690548e14, 9.2114114429666171084e14, -1.7972074771438847361e15, 3.5078816625988685937e15, -6.8495297151644711518e15, 1.3379514680230657677e16, -2.6144400241936380718e16},
    {1, 1, 1, 0.50633999020064525907, -0.69438017667213604531, 0.98599967699780875388, -1.4651463463075830714, 2.2565364446882362145, -3.5772355975449756908, 5.8054310190772562567, -9.6052016864994833701, 16.150393850052390499, -27.528885536279188025, 47.475762466285764141, -82.708444656938887561, 145.36721004075767479, -257.49209066658808571, 459.26236626577368487, -824.20758018816925993, 1487.369504193234552, -2697.5820617790040304, 4914.7687114916236431, -8991.4499853945542026, 16512.098314053713175, -30428.941341658340459, 56255.534932446045974, -104311.59981021280526, 193952.49058152443986, -361551.19068144482565, 675586.0358775971396, -1.2652026465466142966e6, 2.3743604408772985643e6, -4.4646205780373089465e6, 8.4105234106775603681e6, -1.5871415568074590961e7, 2.9999913287724961698e7, -5.6793292055234553228e7, 1.0767401821532755371e8, -2.0442243032799735772e8, 3.8861481259962333208e8, -7.3969897865749838719e8, 1.4096434815615429643e9, -2.6894145014950211641e9, 5.1366225397048116683e9, -9.8208215820920956454e9, 1.8795264625775225714e10, -3.6004849447460804972e10, 6.9034809402779858254e10, -1.3248098727037514173e11, 2.5444972260501746215e11, -4.8910153024821689853e11, 9.4087326316452329633e11, -1.811279722484592991e12, 3.489395186370453167e12, -6.7268792364385091316e12, 1.2976726997260131464e13, -2.5049240551329774343e13, 4.8382927692934685337e13, -9.3507992796343782766e13, 1.808238121687899583e14, -3.4986824106748775731e14, 6.7731011238376454079e14, -1.3118894404880784097e15, 2.5422977149629423057e15, -4.9291071597310628163e15},
    {1, 1, 1, 0.41406249996117346446, -0.50149850149850128539, 0.64502178947238587305, -0.88199011707820265568, 1.2647920527492500631, -1.883480333493030903, 2.8912619336065928848, -4.5499534013640348976, 7.3097834506431895482, -11.950257063760239127, 19.829958432687626496, -33.331618212844747394, 56.658632810891311519, -97.266586375246105995, 168.44621149921964682, -294.00125950536335825, 516.74987412782196061, -914.02517660280107172, 1626.0243992982772637, -2907.806749465840095, 5224.9274733185020458, -9429.7791734081962123, 17087.539532654639328, -31080.032217735504816, 56726.88760888220943, -103871.43389761298004, 190769.02141007242324, -351348.7169882566004, 648798.89604931471622, -1.201026602094174751e6, 2.2284446094009061801e6, -4.143804253881530401e6, 7.7213141718230419047e6, -1.4415465176216116827e7, 2.696288291382688678e7, -5.0519905363944395945e7, 9.4815371745993158893e7, -1.7822939746484517011e8, 3.3552988890250712934e8, -6.3256278588047220784e8, 1.1941750771644904037e9, -2.2573433799324589191e9, 4.272357592574466988e9, -8.0957036392454813479e9, 1.5358076083724743542e10, -2.9167134078877197531e10, 5.5450593369903008696e10, -1.0552516637377689952e11, 2.0101398281431241998e11, -3.8326645507420016291e11, 7.3141609753377144256e11, -1.3970214147750147726e12, 2.6705635294813111512e12, -5.1091804642568485755e12, 9.7821856047046982473e12, -1.8743269672066907801e13, 3.5939221276622395576e13, -6.8959709671989818317e13, 1.3240846354470178549e14, -2.5440181815542563407e14, 4.8910246029824407206e14, -9.4090570756506062967e14},
    {1, 1, 1, 0.34094132157410764421, -0.36500591619672690899, 0.42535294677155897631, -0.53535302664649732815, 0.71492813937869453628, -1.0002179032158071086, 1.4524488548389463454, -2.1742003565963222403, 3.337671415903535519, -5.233648125617236454, 8.3565663807413897303, -13.552956752744749811, 22.281767746554780849, -37.07306555780786914, 62.339948139861757086, -105.82181399906036162, 181.16028240675762808, -312.51411302417519924, 542.85671496390197943, -948.94522892242585053, 1668.4115085927507444, -2948.9389578937204951, 5237.7850558545183865, -9345.189806659423621, 16743.413719432829507, -30115.306815748769426, 54362.949552663081076, -98466.283936069925408, 178915.17704859820172, -326059.51506919579451, 595880.59277005593248, -1.091851640451983693e6, 2.0056049794898347687e6, -3.6927102389844689668e6, 6.8141010965895239144e6, -1.2600423428886073601e7, 2.3346859994324247928e7, -4.3340667157311658682e7, 8.0602145242186194148e7, -1.5015660872632633796e8, 2.8019130573667104177e8, -5.2365476245029091126e8, 9.8013439959778911187e8, -1.837165557936104521e9, 3.4483135152562796751e9, -6.4809147699092510883e9, 1.2195899106394356454e10, -2.2978278800838558072e10, 4.3343782710949026115e10, -8.1850706547826730203e10, 1.5473412768448396311e11, -2.9281991029322730538e11, 5.5468906289338115784e11, -1.0517624734309431208e12, 1.9961362373356756773e12, -3.7918683273193482881e12, 7.2093100523253436427e12, -1.3718259880213348706e13, 2.6125060655782125489e13, -4.9791729821281482989e13, 9.4970330668593388197e13, -1.8127554587702130078e14},
    {1, 1, 1, 0.2822265624991557025, -0.26726421153665735682, 0.2822500168346356815, -0.32704852825993981792, 0.40677668202210587098, -0.53471025960618799387, 0.73457266448259976935, -1.0460109758977890636, 1.5344248790910509045, -2.3078644753565804336, 3.5458912252980553691, -5.5489930429943498722, 8.8235857515933530074, -14.228969016688380802, 23.232658813085514674, -38.356135461407691244, 63.956526941996893992, -107.60326211426380568, 182.51234644505422847, -311.86723651551290322, 536.51619406159516796, -928.73158510982566059, 1616.8849757275557224, -2829.8342424368629236, 4977.0033776174799979, -8793.2577835977419462, 15601.696218054126314, -27791.533383239763362, 49689.204691857347623, -89149.883915867223817, 160471.54871877694769, -289741.75261504772741, 524667.51838557755186, -952681.79476689752804, 1.7343560811069406422e6, -3.1651661020953527206e6, 5.7898618284039829888e6, -1.0614584897767047042e7, 1.9500882912140648031e7, -3.5898690985379611351e7, 6.6211841692922287988e7, -1.2234544748358683287e8, 2.2646438673260016295e8, -4.1989346596846015624e8, 7.7978458299356484378e8, -1.4503660976814265949e9, 2.7015996948984096586e9, -5.0394046907090675259e9, 9.4130064646253657273e9, -1.7605393818883706045e10, 3.2969324902061130827e10, -6.1816031891191247003e10, 1.1603790568659894831e11, -2.1806616074479802148e11, 4.1025057984042947023e11, -7.7262086821676232348e11, 1.4565472067175496362e12, -2.7485836432505338709e12, 5.1916534187927472323e12, -9.815253458917189451e12, 1.8573075130153892285e13, -3.5175521924467524104e13},
    {1, 1, 1, 0.23461360602847632038, -0.19664006045241387512, 0.18823110439019702039, -0.20082731974121805394, 0.23266497471493014904, -0.2873786083076671224, 0.37351135842978069413, -0.50597151799603459672, 0.70927600660445385296, -1.023282379197919688, 1.5129069855433603952, -2.2845040973089262669, 3.5135438559666039689, -5.4916061684283162462, 8.7065851145942891973, -13.980276044460495556, 22.705563681357050277, -37.257232941322963756, 61.706665935549447556, -103.07047360771125339, 173.50070371980709113, -294.14128341001691422, 501.94286552485962931, -861.74695606246630553, 1487.7836012651836789, -2582.0347901339404311, 4502.8981896957038743, -7888.4302449086914879, 13878.148630626839986, -24513.24021011482998, 43460.448513386861226, -77324.502310619596141, 138032.935802955583, -247178.88651541747999, 443945.7196001179607, -799596.28556486168812, 1.4440151482334532023e6, -2.6144187635368807773e6, 4.7449048055204892822e6, -8.6313579893004501151e6, 1.5735617795062819899e7, -2.8747389547351032462e7, 5.2623864417762541069e7, -9.651596094188613057e7, 1.7734231760240383909e8, -3.2642950853735615121e8, 6.0186515188528691385e8, -1.1115079146241571317e9, 2.0559001656662676496e9, -3.8083923934791250186e9, 7.064918008047118893e9, -1.3124272767982845365e10, 2.441317275638827361e10, -4.5470883154821945077e10, 8.4797521817611694098e10, -1.5832699301411939172e11, 2.9595873771798965622e11, -5.538536396118648744e11, 1.037601457009005625e12, -1.9459096765987735056e12, 3.6530606498711938651e12, -6.8646660813149397233e12},
    {1, 1, 1, 0.19570922851562456384, -0.14525009248390591775, 0.12604671890423504457, -0.1238419819033978679, 0.13365162390930462289, -0.1551253606712447962, 0.19075829876215911274, -0.24583322361909826787, 0.32932208863422782594, -0.45574915305194262353, 0.64841330169082779274, -0.94477648182598386627, 1.4054344414728338323, -2.1290953425603984596, 3.2777119648965146492, -5.1188732865849353587, 8.0976825453533219442, -12.959246904761724936, 20.958422479615559046, -34.220620555169566491, 56.365184093750703918, -93.586894718271557363, 156.53972587540875469, -263.63030435399017533, 446.79671481075161957, -761.68287158899326009, 1305.6103788741236707, -2249.4227005137481715, 3894.0775101426509441, -6771.5067894653257623, 11824.860980111046403, -20731.456908504426688, 36482.854749293262767, -64429.278896354557257, 114164.3942170139745, -202934.08077515409735, 361814.50855095965158, -646932.02629182962754, 1.1598791595433609126e6, -2.0849338305641075844e6, 3.7570341880480378535e6, -6.7861396638912146507e6, 1.2285133309794885242e7, -2.2288148246058331251e7, 4.0519626232025111e7, -7.381029939807391363e7, 1.3470807538494267381e8, -2.4629894241193830414e8, 4.5112070832886944722e8, -8.2766598738624122993e8, 1.5209749247117581163e9, -2.7994170071227297961e9, 5.1602072537129173213e9, -9.52571141471685368e9, 1.7609067775010692697e10, -3.2595906123833890151e10, 6.0416788966999626555e10, -1.1212468550107392737e11, 2.0834155781851765285e11, -3.8758286355007612198e11, 7.2185662727653964775e11, -1.3459209114698272619e12},
    {1, 1, 1, 0.16372846233138877116, -0.10764450681598256957, 0.084695752523439593125, -0.076638072292621459531, 0.077050641693442401354, -0.084040640467297161215, 0.097781287825475905038, -0.11988339155831191659, 0.15347561389772533458, -0.20374009865721775509, 0.27894461404286234869, -0.39219111321653375047, 0.56430223489484801875, -0.82857271989672378715, 1.2386165154328646968, -1.8813911985418622335, 2.8989313181906542448, -4.5248063542769926709, 7.1455930057678486325, -11.405033769459266834, 18.381351469444050066, -29.890452880393542228, 49.006615114006341672, -80.960297070913635178, 134.69224618536133423, -225.55364505665858391, 380.01407196372383502, -643.89812693809340298, 1096.8422539611920233, -1877.7479597607341111, 3229.7289214711031446, -5579.7078413867073812, 9679.759382014488365, -16858.755026887463374, 29471.513835169849881, -51702.40834983735643, 91006.484617994469941, -160699.60982191321103, 284623.77199544950292, -505567.4365934253877, 900493.59826459229749, -1.6081322041519014747e6, 2.8790685661929251227e6, -5.1668307757733796778e6, 9.2938202452998225345e6, -1.6754083628050778967e7, 3.0266665297614350073e7, -5.4788496630252789946e7, 9.9371233414692247362e7, -1.8057017361938436539e8, 3.287114038479471559e8, -5.9942937781870932972e8, 1.0949352120943943943e9, -2.0032727192883248757e9, 3.6708653482707862712e9, -6.7367477764394902823e9, 1.238122137486749784e10, -2.2786995719208434721e10, 4.1995337479270305788e10, -7.7497288681815320394e10, 1.4319415020735759874e11, -2.6491108543316580097e11},
    {1, 1, 1, 0.13731002807617245675, -0.079998435432622863026, 0.057076340189332427268, -0.047568532537063377116, 0.044555267825510723957, -0.045670155503234880305, 0.050277820719810351605, -0.058645538520512384944, 0.071750129961821253052, -0.091368793474534019993, 0.12038149183642120015, -0.16332259471666709578, 0.22729839335971466551, -0.32348389066104371064, 0.4695611455954245793, -0.69370595650354666662, 1.0411392160320223331, -1.5849509267590625988, 2.4440790590807180038, -3.8133274032620948024, 6.0137333615473439839, -9.5774842433993838736, 15.391732409441030462, -24.943224044603219818, 40.736146337695899242, -67.008693268020317337, 110.96659319226022568, -184.91431383206088891, 309.95009384968718733, -522.39376970471446353, 885.00445744089911004, -1506.6176820332069202, 2576.6228028278871404, -4425.667730081487987, 7632.8336398948188003, -13215.348620071256969, 22965.252346837946601, -40048.287296963345788, 70071.749366733135316, -122992.79378687015885, 216536.06511760431197, -382327.19239493937631, 676922.52178235370051, -1.2016822121072317051e6, 2.1386443085114469753e6, -3.8154005771830439173e6, 6.8226257811590348153e6, -1.2227355148173660343e7, 2.1960648698003849172e7, -3.9523349347276146945e7, 7.1272919299589054705e7, -1.2877327544717277245e8, 2.3309191788682095495e8, -4.2266906577300378392e8, 7.6774762547749162178e8, -1.3968689234045773919e9, 2.5455853911077492477e9, -4.6461306965612120701e9, 8.4926835823568237617e9, -1.5546310024926072869e10, 2.8498279925629139886e10, -5.2311832073466076052e10},
    {1, 1, 1, 0.11539768253791998759, -0.059595693920953327063, 0.038560042584181697651, -0.029601258388251370253, 0.025831853886284408108, -0.024884127249910599589, 0.025921113164862781882, -0.028765701214898231778, 0.033633806344924294882, -0.041086028279538424801, 0.052093124474031492974, -0.068198866668757017986, 0.09180498206957758477, -0.12663755095815564666, 0.1784996905316204056, -0.25648602830106124736, 0.37495002358876109057, -0.55670726361161278787, 0.83827952186494968542, -1.2785258144871271431, 1.9729220426057845076, -3.077304446860862425, 4.847539109577550207, -7.7061059796909984086, 12.354352585584475054, -19.962544543122583396, 32.492937179652532799, -53.251112878258231488, 87.830380373973466154, -145.73525779065927962, 243.18194417546908577, -407.94502283064539153, 687.77213167619046942, -1165.0377385671769681, 1982.3368215533603736, -3387.3117301184422058, 5811.3785688371191478, -10008.35008474911005, 17299.139454674193989, -30004.770399257530706, 52214.391654241577158, -91150.548437078622173, 159601.4582808671829, -280263.29690067667718, 493508.12058334853291, -871307.18526635358948, 1.5422346316133758348e6, -2.7364518727869638422e6, 4.8667832703388291969e6, -8.6750907880691128307e6, 1.5496980354529048589e7, -2.7741274424056621584e7, 4.9759858579674163839e7, -8.9428131288432462105e7, 1.6102080118199405624e8, -2.9045219327045425074e8, 5.2483845135896588063e8, -9.4996962810772997583e8, 1.7222766154984299238e9, -3.1273923922685348428e9, 5.6875652587645906771e9, -1.0358914170338160083e10},
    {1, 1, 1, 0.097160577774047894284, -0.044489304711848987499, 0.026107393718048188287, -0.01846158833762984492, 0.015010531576627374659, -0.013589618028839288944, 0.013394727486675770096, -0.014142436863164101454, 0.01580314865306305773, -0.018518626233765637399, 0.022595543569549581401, -0.02854516673891972808, 0.037167478622618648066, -0.049693708911697481333, 0.068016405376967165413, -0.095057088600432767664, 0.13535425210703964899, -0.19600771778079049512, 0.28820322057419633039, -0.42968709643120620566, 0.64880476153540112791, -0.99112559009594692702, 1.5303653825118023032, -2.386481574255137904, 3.755791090397898787, -5.9613229702574183719, 9.5373619828804883177, -15.371985185446484702, 24.948323300610968138, -40.754479600631123886, 66.982554156111916022, -110.72485914189644273, 184.02787307919895464, -307.43042751067082374, 516.07771806894215957, -870.31815149828059726, 1474.1231052704816418, -2507.1941164362304171, 4281.0790491836661406, -7337.5051558775308508, 12621.12084170196654, -21783.684764973791168, 37720.9770548000782, -65522.615438301828013, 114155.97779261117105, -199457.90858167834241, 349460.70851802707426, -613892.44233764573431, 1.0811556157022560862e6, -1.9087271222040878974e6, 3.377686267390843107e6, -5.9906910290121438789e6, 1.0648319079565205146e7, -1.896697138828715919e7, 3.3852891204684011907e7, -6.0540231950053594048e7, 1.0847113044104125995e8, -1.9470571068243768909e8, 3.5011575842936937585e8, -6.3064990700351426991e8, 1.1378503840881820407e9, -2.056267786886630105e9},
    {1, 1, 1, 0.081937911911672368206, -0.033273059118077850774, 0.017709990269181201228, -0.011536578244680296405, 0.0087397455282398924425, -0.0074364138128112906458, 0.0069357340082769337077, -0.0069671814817587201521, 0.0074404467863708110556, -0.0083640104240459279442, 0.0098210789341494660101, -0.011972488130544960578, 0.015078508683869035779, -0.019540718081146762645, 0.025971176178057423659, -0.035302738019545512606, 0.048963753811181020512, -0.069155195062790343547, 0.09929219663057919341, -0.14471125444303326262, 0.21380906117507064791, -0.31988627643601708709, 0.48414817363671405833, -0.74061269543888884, 1.1441774392589525302, -1.7839408013174318491, 2.8052984689035918876, -4.4467548270845100956, 7.1015040349393430163, -11.420860693706405486, 18.488650369719615877, -30.116333958825208774, 49.344236230348363409, -81.295742088624178208, 134.63784467064935057, -224.08643567678422787, 374.71687165004016718, -629.40267335620651982, 1061.6891443316088633, -1798.1343306544021375, 3057.1830158281975483, -5216.9828749811888948, 8933.9848475415622161, -15350.856229676140109, 26461.811319779118508, -45755.989083886125215, 79352.960771222090623, -138011.02850017368207, 240686.42163409960253, -420853.5368501419548, 737749.61175817647136, -1.2964179039708298283e6, 2.2834982689806796835e6, -4.0312492804592120196e6, 7.1322654141043154781e6, -1.2645370593806656422e7, 2.2465746614157305356e7, -3.9991337788833783526e7, 7.1324508418201143798e7, -1.2744205613041519165e8, 2.281195445535548642e8, -4.0903810590717471361e8},
    {1, 1, 1, 0.069199353456497190463, -0.02492497222920343137, 0.012033897722087095248, -0.0072216354303137387169, 0.0050975717881861449332, -0.0040765200730976272598, 0.0035977132733782204158, -0.0034385081704644259029, 0.0035094416794035393322, -0.0037844858376044318677, 0.0042764611257589471236, -0.0050306857040742801538, 0.0061283892852068308353, -0.0076979268568932630203, 0.0099349346339129547578, -0.013134963436498472749, 0.017744963862091113804, -0.024444151485421462327, 0.034271312608154042688, -0.048826137541625754641, 0.070589306839302498865, -0.10343413068234208118, 0.15344875186465247616, -0.23026430956848982998, 0.34921124026840765624, -0.53483712210518932111, 0.82667271888434252745, -1.2887257679962880862, 2.0251806229572744453, -3.2064693678530796555, 5.1127411084828228204, -8.2066223281331305388, 13.255463020190428757, -21.537468404795222829, 35.190491875105382392, -57.804209857929214453, 95.428771103774377393, -158.29820882781462631, 263.7840279169563239, -441.47185412464814656, 741.91184019011620975, -1251.742498061966622, 2119.8989938079567414, -3603.1456342156134753, 6145.3759612586842324, -10516.052686829097026, 18052.450844008398575, -31084.485243866137991, 53681.355968691354147, -92966.555926021172822, 161438.60995052462083, -281074.83676746508146, 490601.95094676039224, -858401.31755929853093, 1.5054570885198500262e6, -2.6462342807702818613e6, 4.6616201929713854238e6, -8.2292944810288860509e6, 1.4557118905687578391e7, -2.5801609591797089626e7, 4.5819411547008938357e7, -8.1518770108300024726e7},
    {1, 1, 1, 0.058516072473299034634, -0.018698425967719565003, 0.0081893165286372846847, -0.0045275422379006386874, 0.0029778720366884498915, -0.0022382084407525521554, 0.0018691776881328260524, -0.0016997175988047057012, 0.0016579593198567016466, -0.001715135785901713988, 0.0018651390445895907828, -0.0021172543882578491246, 0.0024948179469188858267, -0.0030374761413135211642, 0.0038066710908476069834, -0.0048950515244157275195, 0.0064414600809598808343, -0.0086543545278779987133, 0.011848315562207886572, -0.016501115848236691363, 0.023343349817377463874, -0.03349994529582072017, 0.048714793565191377437, -0.071709217735014579042, 0.10675699249403352474, -0.16061137584274222684, 0.24400691616912146919, -0.37410383479713767273, 0.57848460878080863122, -0.90171699738426774226, 1.4161770636776556301, -2.2399699827778843112, 3.5667199674621061776, -5.715278056787373908, 9.2129685147397704386, -14.935495218454299303, 24.342874310642522629, -39.878612093353222833, 65.647228678781532377, -108.56775276844910477, 180.34340133990630928, -300.83462250806748921, 503.85143057708928535, -847.12700855704553917, 1429.5355872486947138, -2420.8913768471471398, 4113.6474194418246569, -7012.8013344787809461, 11992.602482542295206, -20570.302482188818255, 35385.399350604044875, -61040.387160355886698, 105578.70008163004658, -183087.8980526507373, 318293.58899897960112, -554681.55853921664096, 968883.81356902882441, -1.6962041994871847435e6, 2.9759868213359586066e6, -5.232386682395759413e6, 9.2184008494275097441e6, -1.6273111672509771522e7},
    {1, 1, 1, 0.04953911760821938498, -0.01404556909416572729, 0.0055805237967120084941, -0.0028424219744794633507, 0.0017420324628161726442, -0.0012306226142688047244, 0.00097250621551600731156, -0.00084140399505948772736, 0.00078439250703418904756, -0.00077842360579616410693, 0.00081463963178622386623, -0.00089237629143969409306, 0.0010170955885908530134, -0.0012002824695819829732, 0.0014606902911163550176, -0.001826917586179471947, 0.0023416800456415366259, -0.0030685209877449386819, 0.0041022090071962202202, -0.0055848317393514181879, 0.0077307880267656900489, -0.010865784789838080237, 0.015488007600769743783, -0.022364587721029365057, 0.032684537817115247709, -0.048302463661848537501, 0.072128879276180051499, -0.1087583477138019252, 0.16548510287170978931, -0.25395262598158697126, 0.39284463511357720139, -0.61229353936366417064, 0.96113202875476546084, -1.5188685471605816299, 2.4155394072283489136, -3.864739648367965684, 6.218777202695745628, -10.061085389056302258, 16.361578725443696558, -26.738664990141328073, 43.902513292362942245, -72.407229732981259441, 119.93093947272482261, -199.46043649854646904, 333.02981550318198015, -558.13542428117999829, 938.77127460687673989, -1584.4601401834693542, 2683.1528643447713354, -4558.2350576371443101, 7767.5307252976481757, -13275.621411185670958, 22754.416031524452859, -39108.51454515167343, 67395.343599612834386, -116439.87284517219269, 201673.64184933542532, -350135.66889574124562, 609297.4085035656936, -1.0626632874615805017e6, 1.8573961548010871104e6, -3.2533180208778021605e6},
    {1, 1, 1, 0.041983009340329622701, -0.010562886691557685918, 0.003807420436343224454, -0.0017867124952553092665, 0.0010203611650301472629, -0.00067748839062685632327, 0.00050662961017895873723, -0.00041705338810021370645, 0.00037158192598575038328, -0.00035375038834854632947, 0.00035627457634967941481, -0.00037660776991685152338, 0.0004151949429340168803, -0.00047492209601952565172, 0.00056122897599258073494, -0.00068273242377592523749, 0.00085239583677811975952, -0.0010894181832761570783, 0.0014221669397907167855, -0.0018926877167422613439, 0.0025636389834284667841, -0.0035289946597165626259, 0.0049306422653646288647, -0.0069842598278246151786, 0.010019873630484475052, -0.014545760013155750972, 0.021349644336832275272, -0.031659765449074844872, 0.047402459253262178601, -0.071616029832540397882, 0.10911874775695411266, -0.16759174189672601176, 0.25934191788257318022, -0.40418371076562417784, 0.63416842291420869216, -1.0013754711691928899, 1.5907955478026689585, -2.5417099698927672856, 4.0832925532557877067, -6.5940948458771665931, 10.701760299481947407, -17.450700408104210673, 28.584914179341085575, -47.026430571352791315, 77.687003235228823374, -128.84900435282445214, 214.52099315606123083, -358.46640637093338719, 601.11121837109384133, -1011.4170504171297846, 1707.3377588085211781, -2891.1464704799738831, 4910.5801157216059715, -8364.9006934204877213, 14289.260102296972332, -24475.841205763775362, 42034.370521712990739, -72372.337380276063777, 124912.45234100440081, -216107.40625642226133, 374741.39951329832337, -651269.38127384044118},
    {1, 1, 1, 0.035613117215689271689, -0.0079522073628554073707, 0.0026005545533497107682, -0.0011243725134781395398, 0.00059833928641783864385, -0.00037340425102593119079, 0.00026423631420012838841, -0.00020695933917451297446, 0.00017623169797578832545, -0.00016094881336787694525, 0.00015599669634430699271, -0.00015912665350215798308, 0.00016968978693560448818, -0.00018813754996187974203, 0.00021589219784171925353, -0.00025544514711841863401, 0.0003106498900941363691, -0.00038723692612026319261, 0.0004936287518124826496, -0.00064219262795257234217, 0.0008511533639912601787, -0.0011475169031577129018, 0.0015715562556151273212, -0.0021837285805587667084, 0.0030753969620716129923, -0.0043855375645113280831, 0.0063269093251954079309, -0.009227254490704014151, 0.013594484708573405289, -0.020220312526241054932, 0.030345771192660423276, -0.045926777195200734698, 0.070062068370088200063, -0.1076857243861352471, 0.16669245476530804352, -0.25977336721651231894, 0.40742232546799380177, -0.64287768394934242448, 1.0202746466722668187, -1.6281412097016483941, 2.6118159063395129616, -4.2108075489887823573, 6.8212542493110435487, -11.100666078767196125, 18.144101793580467481, -29.78135655619247441, 49.079689823931527427, -81.196569405492816661, 134.8299560479297823, -224.69118375088303932, 375.73199727654316862, -630.38759862922080156, 1061.0168885345725992, -1791.3180565082669429, 3033.2771294263443853, -5151.0537608422476597, 8771.6754962660012207, -14977.229784985523234, 25639.221567892409016, -44001.387079234297616, 75697.509794699879687, -130532.19208984902355},
    {1, 1, 1, 0.03023583290374609317, -0.005992578232490081172, 0.0017780242863392029146, -0.00070829042611361851784, 0.00035123070289825903983, -0.00020602124136061335647, 0.00013795984690562793553, -0.00010281084045297382702, 0.00008367115088246181497, -0.000073306475093931642609, 0.000068377128221096277728, -0.00006730734563372683564, 0.000069426629350724029461, -0.000074609829845269921992, 0.000083138446055544610309, -0.000095678296919886187795, 0.00011333661445713181941, -0.00013779342546011990909, 0.00017152218519573430814, -0.00021813332026358584287, 0.00028289771474509274178, -0.00037354088525014159589, 0.00050144996049045164467, -0.00068351556579046744454, 0.00094495652261234986504, -0.0013236744120737354205, 0.0018770018182472296759, -0.0026922141590445576701, 0.0039029869313897143666, -0.0057152879313096524168, 0.0084483084518689231543, -0.012599469044915481495, 0.01894812490604438946, -0.028721732491531247271, 0.043863225017704829067, -0.067463007426642710017, 0.10445970220173450999, -0.16278122453795060652, 0.25520982493573331859, -0.40244147701483547579, 0.63812208688128885704, -1.0171660946252235345, 1.62954219939382345, -2.6231925189918934311, 4.2422548288829775122, -6.8909980519760105922, 11.241082271154394829, -18.412012770509939693, 30.275573658763439054, -49.970797344275866869, 82.77733826261895035, -137.60043935776397667, 229.50198001323311501, -384.02483476239534446, 644.59833595591200044, -1085.2488250292990662, 1832.4637645802969285, -3102.8819951013234511, 5268.4013196333484385, -8968.8763558729807888, 15307.578779284382378, -26190.852280314978466},
    {1, 1, 1, 0.02569084193601156585, -0.0045198608203655681576, 0.0012167709011302411461, -0.00044660128577326338255, 0.00020637162604542702332, -0.00011377866646321456205, 0.000072099390970602928143, -0.000051122662475072092399, 0.000039763966418665911444, -0.000033421073852968400871, 0.000030000663740789216382, -0.000028497543616721402932, 0.000028432996670250736799, -0.000029617161705345196613, 0.000032047506260311357263, -0.000035872124012300015938, 0.00004139020522524964128, -0.000049080494746067276191, 0.000059658071014180290683, -0.000074166549068353633669, 0.000094119738647238935336, -0.00012171584277824183386, 0.00016016050309067306718, -0.00021415515479928641595, 0.00029063843554356723096, -0.0003999173743201504531, 0.00055740224012728482021, -0.00078628104018089376417, 0.0011216637586843525873, -0.0016170362277698857337, 0.0023543599433989578506, -0.0034599526502296217439, 0.0051295737454469389717, -0.0076682276146561056662, 0.011553595104695942982, -0.017537546481497881116, 0.026809264478663156768, -0.041258412613360736247, 0.063901346604186568537, -0.099573965754743042739, 0.15606210866342957061, -0.24595231165010808385, 0.38967236292768825045, -0.62050312870984141751, 0.99286631825211252728, -1.5960723867334584474, 2.5771951553249744363, -4.1792444279261388065, 6.8050502885568907391, -11.124475571203500446, 18.254830411092024635, -30.065269139591771593, 49.691686559451443432, -82.40984967676464273, 137.11959277124778024, -228.87367943774135714, 383.19652198493814668, -643.47605039919536271, 1083.6430804423825524, -1829.9665001512677797, 3098.5954601235017, -5260.355925166913406},
    {1, 1, 1, 0.021844992591900714933, -0.0034118546001680012658, 0.00083338705143138630752, -0.00028183996382159111595, 0.00012136293626116074392, -0.000062891421461249660095, 0.000037713299556419219225, -0.000025443322622745770069, 0.000018914325219774947714, -0.000015250594958132075228, 0.000013174681429330127094, -0.000012076537192375693143, 0.000011654935774236620538, -0.000011767442719538625626, 0.00001236455090697523721, -0.000013461484447844568407, 0.000015129254102840305497, -0.000017497753972577050393, 0.00002076879794708460943, -0.000025239901893312496394, 0.000031341932287562918602, -0.000039696305629002360679, 0.000051200885482648011174, -0.000067158822724241539405, 0.000089472354215964011349, -0.00012093559177142801414, 0.00016567906783839943175, -0.00022984817607200691229, 0.0003226439086797461167, -0.00045792747439541304527, 0.00065670661429355162966, -0.00095100676893283661907, 0.001389926900607264602, -0.0020491563593922423664, 0.0030459976142357234819, -0.004563183043231462134, 0.0068867927545621242597, -0.010466869145689159009, 0.016014699004842603728, -0.024659546965644108287, 0.038202118764694117811, -0.059525943982012779156, 0.09326742958597379711, -0.14691097916399594438, 0.23258479816455688059, -0.37001524338801993485, 0.59140243961014432138, -0.94949121459445714066, 1.5309711872922742656, -2.4787891080820914722, 4.0294053234542489805, -6.5751735503998858226, 10.769063218709262684, -17.700920847328950339, 29.194888000695273297, -48.312485206803467712, 80.205611733981685856, -133.56620379653064813, 223.09548111259859228, -373.7191669964311886, 627.79870764816679785, -1057.4937725727854359},
    {1, 1, 1, 0.018587394860787753714, -0.0025774041653940592598, 0.00057124543937543412346, -0.00017800398276453319888, 0.000071428380361550358146, -0.000034791523524245270644, 0.000019742902344178404945, -0.000012673287739889189078, 9.0042668831286810996e-6, -6.9648301746636786627e-6, 5.7903867124839831728e-6, -5.1219614176861477186e-6, 4.781416885141840731e-6, -4.6792983398679122364e-6, 4.7744448738915845153e-6, -5.0557970935854056819e-6, 5.5347577942136777574e-6, -6.2433424066841760266e-6, 7.2362783833223306737e-6, -8.5966501409620556805e-6, 0.000010445591772232615874, -0.000012957329479188215867, 0.000016381819599039232142, -0.000021078534307266038682, 0.000027566881170416411004, -0.000036601686325765427222, 0.000049286707917186654491, -0.000067246189877263939838, 0.00009288545940162868689, -0.00012978880773042165238, 0.0001833300325278905414, -0.00026161392254047882081, 0.00037693505248406267529, -0.00054804872272707597201, 0.00080372231075653981699, -0.0011883136114403416352, 0.0017705709240402507539, -0.0026575738570952971802, 0.0040169079932290795277, -0.0061120764633358945197, 0.0093592680316812887286, -0.014418702623174791436, 0.022342150632579974442, -0.034812013798281646371, 0.054530130657406182633, -0.085852186379487158978, 0.13582623859024495935, -0.21589815410868192044, 0.34472083826834235201, -0.55279549191275163874, 0.89016195373112603577, -1.4391773016184673643, 2.3358078139707636676, -3.8052011858102446118, 6.2212717245925937149, -10.206762156740067604, 16.801694071284720248, -27.747632096008315181, 45.968513832317164227, -76.385831513475595168, 127.30374706846084066, -212.76773240486320319},
    {1, 1, 1, 0.015825476716011639242, -0.0019483934987647367808, 0.00039184231453185433133, -0.00011250590058533757095, 0.000042070498404165123252, -0.000019261058599263793716, 0.000010343177802198429892, -6.3173175570003462747e-6, 4.2897790446263136536e-6, -3.1832022647272411187e-6, 2.5468637868793263158e-6, -2.1740098551446620922e-6, 1.9630673238299539002e-6, -1.8621372689502384348e-6, 1.8450159777404422123e-6, -1.9002877956812349828e-6, 2.0263436244626665404e-6, -2.2293879416709299136e-6, 2.5232080646484395213e-6, -2.9302519179617160961e-6, 3.4839719395427569438e-6, -4.2326798066200701928e-6, 5.2454346078138071613e-6, -6.6208320168122683112e-6, 8.5000478730627203545e-6, -0.000011086211388477929149, 0.00001467327985433657132, -0.000019689271089917057991, 0.000026761314232269667374, -0.000036814026864365379967, 0.000051219044557526214284, -0.00007202343264572358174, 0.00010230029643337214909, -0.0001466895307478403797, 0.00021223570081445678771, -0.00030969220353024301904, 0.00045556014842473790182, -0.00067528956090026776357, 0.0010083265245325942185, -0.0015161030544178072818, 0.0022947354824086321934, -0.0034952837707902837912, 0.0053561946724883846483, -0.0082554451557163462371, 0.012794645299450378705, -0.019935156305666035198, 0.031219131579440219513, -0.049129633917512642379, 0.077679194629858623632, -0.12337469695681904889, 0.1968039575444803288, -0.3152522107496333481, 0.5070293254154085753, -0.81864612108453692982, 1.3267478056937905692, -2.1580101866638728848, 3.5223969497845454702, -5.7688893004946757702, 9.4791001358980277371, -15.624899476305670043, 25.834434718141654159, -42.842099093531679387},
    {1, 1, 1, 0.013481792275470638742, -0.0014738428995269865082, 0.00026896131849513613382, -0.000071156746415240331642, 0.000024796081166655516698, -0.000010670571755702283034, 5.4224991264378350923e-6, -3.1512285662939317787e-6, 2.0451552234379526466e-6, -1.4558731546843727702e-6, 1.1210116749120739044e-6, -9.2340780995826881692e-7, 8.0653103534772619242e-7, -7.4156695038986219225e-7, 7.1348599811180967596e-7, -7.1475550909584939853e-7, 7.4239704568923675783e-7, -7.966417378143503303e-7, 8.8044047372276282963e-7, -9.9951657916718016461e-7, 1.1628555408855541627e-6, -1.3836460716111437693e-6, 1.6807791452513583003e-6, -2.0811084062066138357e-6, 2.622800219459016677e-6, -3.3602797607396164654e-6, 4.3715446707880240118e-6, -5.7690197584346757234e-6, 7.7157414892315435034e-6, -0.000010449608563640677756, 0.000014319901435776793283, -0.000019842553321810816618, 0.000027784213538702727993, -0.000039290720002115401989, 0.000056084369123445878004, -0.000080768221124889836219, 0.00011729762722876811906, -0.00017171407237406921172, 0.00025329213813955423691, -0.00037633960745483120609, 0.00056303408297176087613, -0.00084791057068471167349, 0.0012849881638770416908, -0.0019591301298435186882, 0.0030042179726980650083, -0.0046323295196946560939, 0.0071807451394315368051, -0.01118792927031049828, 0.017516749106191354113, -0.027554932745504722495, 0.043542193897878422765, -0.069105671444603893807, 0.11013891374570076378, -0.17624895844414030398, 0.28314531355941354374, -0.45659463428989877359, 0.73898455566249532154, -1.2002462249946346325, 1.956075653001597235, -3.1984057967993465577, 5.2464878499735904756, -8.6327203064439312264},
    {1, 1, 1, 0.011491425285371504744, -0.0011155437842634089305, 0.00018473021327719583824, -0.000045033042752615825016, 0.000014623990506916915616, -5.9152725847100631327e-6, 2.8446316663351492675e-6, -1.572929668899616194e-6, 9.756650299571247108e-7, -6.6629520056020594978e-7, 4.9374078267685311706e-7, -3.9247367602271030365e-7, 3.3158307821817275251e-7, -2.9551177928357767333e-7, 2.760939289592148383e-7, -2.690184128518505521e-7, 2.7217359437254726025e-7, -2.8485726572347073427e-7, 3.0742131099401591757e-7, -3.4116326576104506882e-7, 3.8838662527976403925e-7, -4.5260791333601697427e-7, 5.3892395834578470257e-7, -6.5458286048664831073e-7, 8.0983577499360268079e-7, -1.0191913537340317413e-6, 1.303259060832950092e-6, -1.6914634527897859117e-6, 2.2260565315694258627e-6, -2.9680764942377548278e-6, 4.0062406122903150264e-6, -5.4702832027415200539e-6, 7.5510597181307112342e-6, -0.000010530997225542450552, 0.000014830436712350419087, -0.000021078491192733657364, 0.000030221884242208155121, -0.000043692872680480258559, 0.000063669455366807361233, -0.000093480297049011414932, 0.00013823745633333330635, -0.00020582906123158828642, 0.0003084828147942664015, -0.0004652380015439213084, 0.00070586845720447106181, -0.001077130752066368231, 0.0016527509672705199556, -0.0025494420199828753881, 0.0039526792898241547286, -0.0061583152027326429337, 0.0096399796459538043392, -0.01515858269322220457, 0.023940757116520410088, -0.037970500004307877811, 0.060467192883613288075, -0.096671302325243491129, 0.1551393267019455722, -0.24988373462082581857, 0.40391845867394072341, -0.65514801914726835294, 1.0661736265214627311, -1.7406608983828923118},
    {1, 1, 1, 0.0097998673341361808298, -0.00084482247528329658966, 0.00012695138126397472402, -0.000028516898039743159584, 8.629930113593054816e-6, -3.2811210890086848591e-6, 1.4931876524747987375e-6, -7.85600137716593648e-7, 4.6573499648835086702e-7, -3.0512245280987806176e-7, 2.1759692266355043155e-7, -1.6691403813690725604e-7, 1.3640470973736623065e-7, -1.1783255501726064409e-7, 1.0690416333108368091e-7, -1.0131485256465918073e-7, 9.9844141082384406515e-8, -1.0191980668641231948e-7, 1.0740763847187854502e-7, -1.1652038569700515608e-7, 1.2979868811582156003e-7, -1.4814498149854294555e-7, 1.7290682148622488566e-7, -2.0601671842173609319e-7, 2.5020542304381476719e-7, -3.0931724634789798833e-7, 3.8877182638955358951e-7, -4.9623968562863300505e-7, 6.4263306870988174332e-7, -8.4356507093722588897e-7, 1.121508396585479168e-6, -1.5090048266154776416e-6, 2.0534603648301463076e-6, -2.824344628450444133e-6, 3.9240521729713784923e-6, -5.5043660272387520959e-6, 7.7915287634199759085e-6, -0.000011124594385617980651, 0.000016014355940893798142, -0.000023234278878795590281, 0.000033961415158331279106, -0.000049995667405566767646, 0.000074102334657523196771, -0.0001105493468732560058, 0.00016595305220418434403, -0.00025061471632002127923, 0.00038064008432915813817, -0.00058131267039301093292, 0.00089248082870901668145, -0.0013771896398869330521, 0.0021355572879283461917, -0.0033271534193591487077, 0.0052072007276982396755, -0.0081853158422738317805, 0.012921104161654039051, -0.020480178906236425088, 0.032589513218902637632, -0.052056513482429925828, 0.083458624400166394967, -0.13428106821429342216, 0.2167987148695506656, -0.35119651517120539477},
    {1, 1, 1, 0.008361277029909185384, -0.00064013552644662050983, 0.00008729148755869514054, -0.000018068090269082399126, 5.095539264201114706e-6, -1.8210126692631254333e-6, 7.8423677940668339911e-7, -3.9258998538438024164e-7, 2.2244526999407754721e-7, -1.3980684058587834118e-7, 9.5951941650180667489e-8, -7.102690366575224275e-8, 5.6145430281775940217e-8, -4.7011503391748972529e-8, 4.1417224608910826839e-8, -3.8177991126929349849e-8, 3.6647827378953505823e-8, -3.6487080601972405268e-8, 3.7547909738459973899e-8, -3.9819069931899459417e-8, 4.340362699004193295e-8, -4.8517850287696198323e-8, 5.55068625778798162e-8, -6.4876934392679212321e-8, 7.7347555549021205835e-8, -9.3929658077067076788e-8, 1.1604036849557220905e-7, -1.4567020302795498061e-7, 1.8562667199892283082e-7, -2.3989024162183565841e-7, 3.1413663319462258272e-7, -4.1650676496127600092e-7, 5.5874726622806938237e-7, -7.5790818380236295546e-7, 1.0388823598513020721e-6, -1.4382219155387127696e-6, 2.0099010319683081055e-6, -2.8340575574457236864e-6, 4.0303127936283293975e-6, -5.778156263215584194e-6, 8.3482765911637740619e-6, -0.000012150917895660213332, 0.000017810818494384163562, -0.000026283807750543959765, 0.000039038924425287174577, -0.000058343955866443984865, 0.000087714785665547856265, -0.00013262507687428663933, 0.00020163106689736518909, -0.00030816040656542974618, 0.00047336660510297543373, -0.00073069874201249918581, 0.0011332404846985907365, -0.0017655332922747564656, 0.0027626814425912129674, -0.0043413150581734258381, 0.0068499159807670628178, -0.010850845613398113376, 0.017254412561633772497, -0.027538583197514314365, 0.044109987806396912719, -0.070898607195714275931},
    {1, 1, 1, 0.007137046617635367518, -0.00048527928578194977342, 0.000060051859892217181349, -0.000011453705270310719386, 3.0102260477231540592e-6, -1.0111866646880867968e-6, 4.1210589227781960279e-7, -1.9629378740908958773e-7, 1.0630111817788991669e-7, -6.4093440088811666728e-8, 4.2333691981912245338e-8, -3.0240198737880249022e-8, 2.3122326342177191063e-8, -1.8766156165299348897e-8, 1.605461703676209675e-8, -1.4394143356682741662e-8, 1.3458816245510065351e-8, -1.3069313952194796001e-8, 1.313317290781738014e-8, -1.3614877716278894639e-8, 1.4521627733503920017e-8, -1.5898268407396520405e-8, 1.7828508157433282231e-8, -2.0441466023021566432e-8, 2.3923814039645025573e-8, -2.8538778290764881629e-8, 3.4654330861671558197e-8, -4.2784270402132600458e-8, 5.3647800798858401006e-8, -6.8256005828863342367e-8, 8.803776162255245596e-8, -1.1502384861428130831e-7, 1.5211742016420269185e-7, -2.034932483573281584e-7, 2.7518991542493586944e-7, -3.7599241487183608043e-7, 5.1875380144730800575e-7, -7.2238348474671820439e-7, 1.0148520322317304302e-6, -1.437752387220348396e-6, 2.0532538314667491567e-6, -2.9547490012322332648e-6, 4.283222999089349704e-6, -6.2525220429520948645e-6, 9.1885147813606199291e-6, -0.000013590019344310992683, 0.000020223951574530970419, -0.000030274463529631873618, 0.00004557754940291137172, -0.000068991394186406547619, 0.00010498300948328204556, -0.00016056057836127844344, 0.00024676007519635216243, -0.00038102322615056840578, 0.00059101306317546329748, -0.0009207547477670324288, 0.0014405480066578248235, -0.002263013929732641347, 0.0035691457670529139988, -0.0056507170951231480453, 0.0089795040898144688113, -0.014320573272462727562},
    {1, 1, 1, 0.0060946179483675058477, -0.00036805435483422890211, 0.00004133211665412757725, -7.2642327777421339687e-6, 1.7791811714872033618e-6, -5.6177608858556839411e-7, 2.1666314375385218467e-7, -9.8194980185887884204e-8, 5.0823961110231267443e-8, -2.9397829481468326134e-8, 1.8686822585765901768e-8, -1.2881418367813651325e-8, 9.5272205669149906685e-9, -7.4948742015140493868e-9, 6.2263992347779625999e-9, -5.4297123779642829922e-9, 4.9452006289654454403e-9, -4.6836565354950078319e-9, 4.5959179856758274959e-9, -4.6575251923949007849e-9, 4.8609765090506716482e-9, -5.2121529029763881467e-9, 5.7293119436119457963e-9, -6.4439622128457041016e-9, 7.4034397778412588647e-9, -8.6753577196702794131e-9, 1.0354410554966259797e-8, -1.2572366796267253119e-8, 1.5512544292478112156e-8, -1.9430714935307078924e-8, 2.4685337578999284647e-8, -3.1781431077192694451e-8, 4.1434514627338679417e-8, -5.4664227989065930973e-8, 7.2932063941193590004e-8, -9.8344969405314968675e-8, 1.33957745995759184e-7, -1.8422429564684580418e-7, 2.5567406966578523214e-7, -3.5793068261222317391e-7, 5.0525253361976274155e-7, -7.188730061085315924e-7, 1.0305702314721281379e-6, -1.4881349137805674424e-6, 2.1637792010346804011e-6, -3.1671194975065294049e-6, 4.6652983350081673081e-6, -6.9142882180069572667e-6, 0.000010307770215034394812, -0.000015453727835911629064, 0.000023294891620193944267, -0.000035298786869835577054, 0.000053758608147141702507, -0.000082271101971565760844, 0.0001264980095926018495, -0.00019538309591872075982, 0.00030310325465705805954, -0.00047220573486568771034, 0.00073866716389277262986, -0.0011600746888406061482, 0.0018288925684651387922, -0.0028940334355140382511},
    {1, 1, 1, 0.0052065015981722778671, -0.00027926771943163738037, 0.000028460541918641192869, -4.6092572944765050333e-6, 1.052060293894372591e-6, -3.1224528091923398289e-7, 1.1396275414006039929e-7, -4.9144444363433920578e-8, 2.4310963474250603355e-8, -1.3490262383553546038e-8, 8.2525583816731470096e-9, -5.4896780446610957915e-9, 3.9274017421756382029e-9, -2.9947323744087913754e-9, 2.4158993001775235258e-9, -2.0491457694270995222e-9, 1.8178839391891720393e-9, -1.6792780570958125847e-9, 1.6090899709415533623e-9, -1.5940511142130262851e-9, 1.6279363126823637647e-9, -1.7095831508000921806e-9, 1.8420260909377667773e-9, -2.0323563210610833573e-9, 2.2921480846726229342e-9, -2.6384293277934054942e-9, 3.0952758992980536079e-9, -3.6962057201789480389e-9, 4.4876644666817825617e-9, -5.5340486447464482832e-9, 6.9249314472787841919e-9, -8.7854776098914013968e-9, 1.1291508719196532031e-8, -1.4691389953639084231e-8, 1.9337974475831622782e-8, -2.5735448381566859473e-8, 3.460835261729441881e-8, -4.7003758742330472644e-8, 6.4443224065686800163e-8, -8.9149806490343885716e-8, 1.2438872878941162897e-7, -1.7498082278351230775e-7, 2.4807969660109143438e-7, -3.5435300860949690717e-7, 5.0978532744705655024e-7, -7.3844067740026906316e-7, 1.0767121971604647737e-6, -1.5798844639662979786e-6, 2.33230493839348373e-6, -3.4632070794538709455e-6, 5.171413329216919591e-6, -7.7640358667884337668e-6, 0.000011717313524775964708, -0.000017772564823754492427, 0.000027088019076993255552, -0.000041479829609651094435, 0.000063805840553701468384, -0.000098578514879349635808, 0.00015294676610222400406, -0.00023827328190344061891, 0.00037267568535186139398, -0.00058513182191887220499},
    {1, 1, 1, 0.0044494621633775771521, -0.00021198605656751134623, 0.000019605675306496680175, -2.9258910430589782072e-6, 6.2237032141562690414e-7, -1.7362710578353734471e-7, 5.9969541313139645924e-8, -2.4606516657923455337e-8, 1.1633942846958417716e-8, -6.1932276785720667765e-9, 3.6461423245364528748e-9, -2.3405731869724618205e-9, 1.6197084464637872158e-9, -1.1971382601607358118e-9, 9.3780706761317514562e-10, -7.7368086378926872878e-10, 6.6856171954864624345e-10, -6.0235626478302949753e-10, 5.6361396566787516916e-10, -5.4581156984151247104e-10, 5.4543732823192152718e-10, -5.6099236841003829289e-10, 5.9249243161935893178e-10, -6.4126939457083432481e-10, 7.0997924788614156718e-10, -8.0278169317369189663e-10, 9.2569374291908160972e-10, -1.0871494635182129876e-9, 1.2988284710309212241e-9, -1.5768532887274357894e-9, 1.9435067489300013073e-9, -2.4296936040194244194e-9, 3.0784771709095707929e-9, -3.9501795281469883372e-9, 5.1297686022354030501e-9, -6.7376069857859117208e-9, 8.9451658709813079809e-9, -1.1998105490889361603e-8, 1.6250333046715748556e-8, -2.2214489501373252671e-8, 3.0637127209863842591e-8, -4.26111481790296069e-8, 5.974469892465035079e-8, -8.4415945727946387734e-8, 1.2015899541006158187e-7, -1.7225084500718247521e-7, 2.4860762807906091741e-7, -3.6115846721882561457e-7, 5.2795945937195941676e-7, -7.7645862488598713824e-7, 1.1485568045864354915e-6, -1.7084802921838701472e-6, 2.5550699096891493355e-6, -3.8410292600946108149e-6, 5.8031731755739161099e-6, -8.8101171471949302166e-6, 0.000013437702564380510779, -0.000020588658602783310254, 0.000031683020829599763842, -0.000048962043187185043044, 0.000075974652762397100027, -0.00011835829565204094798},
    {1, 1, 1, 0.0038038399824686267925, -0.00016097622872393066863, 0.000013511180695938878774, -1.8580648378067519825e-6, 3.6832762293180630424e-7, -9.6586730795348118479e-8, 3.1570221797614359165e-8, -1.2325529422408282426e-8, 5.5697012628687505683e-9, -2.8444238250943889344e-9, 1.6116085505058589807e-9, -9.9834070746326614719e-10, 6.6826666538402374052e-10, -4.7875379962924414651e-10, 3.6419164241993796367e-10, -2.9223538915507759956e-10, 2.4597951606243046209e-10, -2.1615558000339391698e-10, 1.9749925913440495591e-10, -1.8696729039065424815e-10, 1.8282466714975181094e-10, -1.8416463868387783362e-10, 1.9065687767829521364e-10, -2.024248747575890589e-10, 2.2000441898870659947e-10, -2.4436121567056613204e-10, 2.7696071053227602503e-10, -3.1989341467773317525e-10, 3.7606783489108515431e-10, -4.4949267010905939026e-10, 5.4568206119546750697e-10, -6.7223439716929491574e-10, 8.396591431281561483e-10, -1.0625611554294989137e-9, 1.3613435935312208115e-9, -1.7646672390447085187e-9, 2.3130185215704435789e-9, -3.0639101896336512509e-9, 4.099496921375037016e-9, -5.5377785616552364521e-9, 7.5491557966400120615e-9, -1.038100447954802994e-8, 1.4394311989174020909e-8, -2.0118530487062521296e-8, 2.8334050517334721956e-8, -4.0196713446102424408e-8, 5.7426549198748103555e-8, -8.2594997241716884016e-8, 1.1956368833925838775e-7, -1.7415729718315257532e-7, 2.5519914909728736284e-7, -3.7611091926323752286e-7, 5.573924293743063389e-7, -8.3047904483901120057e-7, 1.2437622802356347794e-6, -1.8720178329332349156e-6, 2.8312173049835298317e-6, -4.3018715148537725305e-6, 6.5659334922812429427e-6, -0.000010065314314193896838, 0.000015494944383304523871, -0.00002395120294184998688},
    {1, 1, 1, 0.0032529852288394471916, -0.00012228551560911536707, 9.314686242436127398e-6, -1.180400760102167089e-6, 2.180655882844501883e-7, -5.3750917404967591764e-8, 1.6626229395398714021e-8, -6.1763347853352798621e-9, 2.6675165869482722458e-9, -1.3068994686231432489e-9, 7.1261739310592796039e-10, -4.25996952541688759e-10, 2.7582531227953311769e-10, -1.9153657071131163282e-10, 1.4148752765355731513e-10, -1.1042708448152745329e-10, 9.0537468845268699576e-11, -7.7598182483032892146e-11, 6.923431055766938143e-11, -6.4070895616030242832e-11, 6.1305165617589202672e-11, -6.0482252038101535161e-11, 6.1375441382044458408e-11, -6.3923388568450561733e-11, 6.8200844577299946263e-11, -7.4411446184135495612e-11, 8.2897550960111642958e-11, -9.4165992344049711768e-11, 1.0893147010856873935e-10, -1.2818191749434551521e-10, 1.5327317763638003458e-10, -1.8606421100054932041e-10, 2.2910946536569544941e-10, -2.8593281037272863228e-10, 3.6141877844376309988e-10, -4.6237355134031376719e-10, 5.9833286723514529932e-10, -7.8273084279871103364e-10, 1.0345987447890372214e-9, -1.3810453496082591547e-9, 1.8608949753159100736e-9, -2.5300474511149630495e-9, 3.4694098823041317552e-9, -4.7966850488215833362e-9, 6.6839661077792988083e-9, -9.3841069134545260599e-9, 1.3270406368712010791e-8, -1.8896567886511890182e-8, 2.7087647789153072338e-8, -3.9078538148663333209e-8, 5.6725620309383438968e-8, -8.2831440153048352817e-8, 1.2164455052139364427e-7, -1.7963173615314363771e-7, 2.666751635207180202e-7, -3.9793453215307009892e-7, 5.9675318178406015566e-7, -8.992080519447061633e-7, 1.3612558372657412033e-6, -2.069991389411613382e-6, 3.1614383542257588586e-6, -4.8487452329061328393e-6},
    {1, 1, 1, 0.0027827848357431749074, -0.000092926322557167647788, 6.4238874103144278199e-6, -7.5016215509615583353e-7, 1.2915120406997019397e-7, -2.9923592821909117976e-8, 8.7593129751011160178e-9, -3.096112590804942811e-9, 1.2780370427440642387e-9, -6.0069148522839576439e-10, 3.1522087418141360761e-10, -1.8184280860005269771e-10, 1.1388868019996814875e-10, -7.6657291478469529279e-11, 5.498810779298878639e-11, -4.174274144437899114e-11, 3.3336532967410944105e-11, -2.7867590745990404754e-11, 2.4279522411846838422e-11, -2.1964377728866913674e-11, 2.0564699630715830945e-11, -1.9870680980205014134e-11, 1.9765141623176251349e-11, -2.0193839715431758638e-11, 2.1150048121038109596e-11, -2.266785987196250982e-11, 2.4821532202261691046e-11, -2.7729769952976987257e-11, 3.1564867432681410703e-11, -3.6567425023393643178e-11, 4.30681428346562595e-11, -5.1519126520335749284e-11, 6.2538376479522121025e-11, -7.6972862694170024174e-11, 9.5988076692163403361e-11, -1.2119558052904639154e-10, 1.5483540412164171423e-10, -2.0003802044410231199e-10, 2.6120232341330234283e-10, -3.4454346733636825657e-10, 4.58890517242232113e-10, -6.1685305190075710492e-10, 8.3653495231680995126e-10, -1.1440630222714368421e-9, 1.5773338636920395484e-9, -2.1915893253951746895e-9, 3.0677471893177231556e-9, -4.324899078299530668e-9, 6.1391345770222255269e-9, -8.772004082097129335e-9, 1.261371983022142177e-8, -1.8248971351886405399e-8, 2.6557559828893981234e-8, -3.8868822611430325509e-8, 5.7199429622444279289e-8, -8.4620819100075210428e-8, 1.2582885428356157845e-7, -1.8802989715938177478e-7, 2.8232348695226091496e-7, -4.2586674065170019581e-7, 6.4527287580524742852e-7, -9.8196363508993850737e-7},
    {1, 1, 1, 0.0023812663193902382203, -0.000070639112500173676499, 4.4317410455555401885e-6, -4.7690263295378726318e-7, 7.6517364119493214131e-8, -1.6664516158374823408e-8, 4.61634317342102462e-9, -1.5525833349577491196e-9, 6.1253709977586757495e-10, -2.7619367178741766935e-10, 1.3948474996398171647e-10, -7.7649594156428981048e-11, 4.7041439630271669743e-11, -3.0690855265418515641e-11, 2.1378305849191882624e-11, -1.5784845715801848183e-11, 1.2279102628632718825e-11, -1.001155432826462289e-11, 8.5175207628275934536e-12, -7.5323654274141624842e-12, 6.9008423787914653792e-12, -6.5305844069487207866e-12, 6.3673656217361000691e-12, -6.3816448088560797477e-12, 6.5612656606459007441e-12, -6.9077401656492883824e-12, 7.4348158027319788147e-12, -8.1687052328229275855e-12, 9.1497533953692661081e-12, -1.0435586033784613444e-11, 1.2106009402414829897e-11, -1.4270166635205326653e-11, 1.7076745595222356414e-11, -2.072842217873118537e-11, 2.5502270688452544079e-11, -3.1778660777491917812e-11, 4.0082307354842947253e-11, -5.114082145681695762e-11, 6.59685888743594897e-11, -8.5987473512657617423e-11, 1.1320129877572630368e-10, -1.5044920324939157742e-10, 2.0177516488986925475e-10, -2.7296932670607030613e-10, 3.7236442667457873052e-10, -5.1201252970133530759e-10, 7.0943091222466585987e-10, -9.9020298828738781694e-10, 1.3918685933631428797e-9, -1.9697657174573524508e-9, 2.8058364459876162412e-9, -4.0219514372178171515e-9, 5.8001462650625856316e-9, -8.4134667000135719546e-9, 1.2273151992775370239e-8, -1.8001061568508376448e-8, 2.6541230037126457569e-8, -3.9332263830990234603e-8, 5.8574634292167101189e-8, -8.7646431629028971513e-8, 1.3175206786222620027e-7, -1.9893756153605297952e-7},
    {1, 1, 1, 0.0020382654579225329078, -0.000053713980904374603977, 3.0583712776548566902e-6, -3.0328131491740886804e-7, 4.5348626351110489953e-8, -9.2835719918612587436e-9, 2.4337174766510915129e-9, -7.7882087671343661324e-10, 2.9367447547524648662e-10, -1.2703433166116357656e-10, 6.174245013796336345e-11, -3.3168651867218743933e-11, 1.9436872529294955507e-11, -1.2291655922384147937e-11, 8.3142638939232796486e-12, -5.9709829769223807161e-12, 4.5243796777374608963e-12, -3.5979066801324791214e-12, 2.9890462948048985075e-12, -2.5839869035958938933e-12, 2.3164789482142256879e-12, -2.1470290140595472339e-12, 2.0519476129361826156e-12, -2.0174046765434796377e-12, 2.03615390001024163e-12, -2.1057568584178037683e-12, 2.2277099005336683986e-12, -2.4071713903131213169e-12, 2.6531491892205173318e-12, -2.9791076716557810944e-12, 3.404025396522724974e-12, -3.9539989286121862529e-12, 4.6645589847785994265e-12, -5.5839542437457890065e-12, 6.7777794671878246824e-12, -8.3354958876345853226e-12, 1.0379638270460784603e-11, -1.3078861265861283879e-11, 1.6666501769349402159e-11, -2.1467104707014963638e-11, 2.7934497888201102235e-11, -3.6706689292298309217e-11, 4.8685372178823904148e-11, -6.5151576132559204723e-11, 8.7934628614722969842e-11, -1.1966005724436171727e-10, 1.6411484494605985806e-10, -2.2678781502765784493e-10, 3.1567236040976350473e-10, -4.4246362761432104946e-10, 6.2435096317629917706e-10, -8.8671188147743884034e-10, 1.2671761918136294329e-9, -1.8217797333485356427e-9, 2.6343157457357770718e-9, -3.8305957206544009563e-9, 5.6002726062585972144e-9, -8.2303500098600665527e-9, 1.2156806276188831151e-8, -1.8044385471471934629e-8, 2.691031735610299963e-8, -4.0316751484094301729e-8},
    {1, 1, 1, 0.0017451471151975722143, -0.000040856246465855957933, 2.111244680220605035e-6, -1.929283178527154569e-7, 2.6884614127204237484e-8, -5.1733733782136115941e-9, 1.2834504636552827762e-9, -3.9080269518781305543e-10, 1.4084376961227339732e-10, -5.8447561027219037684e-11, 2.7338778981311678225e-11, -1.4172770577073818297e-11, 8.0336079191335287589e-12, -4.9243668131573788704e-12, 3.2345438114862387487e-12, -2.2593844801681008736e-12, 1.6675941529783297241e-12, -1.2934131392683040049e-12, 1.0492793877914627999e-12, -8.867237882166577985e-13, 7.7784630122708372436e-13, -7.060948053627060096e-13, 6.6147286793741367697e-13, -6.379589705723843042e-13, 6.3208131236892071184e-13, -6.4212542258999840924e-13, 6.6770788061098201796e-13, -7.0957822480561870353e-13, 7.6957945086149245308e-13, -8.5073655519795747897e-13, 9.5746771862566088448e-13, -1.0959320792084003847e-12, 1.2745464946307442511e-12, -1.5047248439372399887e-12, 1.8019206921979727783e-12, -2.1870915759734439468e-12, 2.6887561569554687406e-12, -3.3458916925400120063e-12, 4.2120298504393785075e-12, -5.361070341847036464e-12, 6.8955686999893548349e-12, -8.9586034695615764144e-12, 1.1750844153353081718e-11, -1.5555207412400336643e-11, 2.0772630561721771221e-11, -2.7974198363691042534e-11, 3.7977420513792647189e-11, -5.195831400594153487e-11, 7.1616771620366077471e-11, -9.9421531035157166408e-11, 1.3897449344729609094e-10, -1.9555463881258155414e-10, 2.7693316047486899019e-10, -3.9459961726797291755e-10, 5.6561309848722433333e-10, -8.1540701727785554188e-10, 1.1820537387295142324e-9, -1.7227712075312579851e-9, 2.5238837807366945542e-9, -3.7161215118502017604e-9, 5.4981952584958938923e-9, -8.1732397754935209868e-9},
    {1, 1, 1, 0.0014945703859614061844, -0.000031085114817212963493, 1.457850906048001977e-6, -1.2276490144204750852e-7, 1.5943088905846887467e-8, -2.8837808516240803904e-9, 6.7704619531501347333e-10, -1.9615896125047567666e-10, 6.7567838584854496385e-11, -2.6899417576578701669e-11, 1.210893007208024261e-11, -6.0577768115541498944e-12, 3.3214414252923186563e-12, -1.9734323454556346685e-12, 1.2587347588294229406e-12, -8.5519764961258550473e-13, 6.1482825509271543642e-13, -4.6511125661842730853e-13, 3.6845284737521961905e-13, -3.0438184601518294892e-13, 2.6127124497236306203e-13, -2.3228467602850942046e-13, 2.1329972103521201937e-13, -2.0180178423483991173e-13, 1.9627630193296252762e-13, -1.9586827883352609388e-13, 2.0019213155940045144e-13, -2.0923108705527531602e-13, 2.2329445978188047221e-13, -2.4301703018467138575e-13, 2.6939419594905137624e-13, -3.0385298013437310636e-13, 3.4836423945314712656e-13, -4.05606721543557761e-13, 4.7919992553528315084e-13, -5.7403099537233006235e-13, 6.9671230945344027052e-13, -8.5622264028974431103e-13, 1.0648080522332603167e-12, -1.3392524520271411462e-12, 1.7026768748431560075e-12, -2.1870985551323724428e-12, 2.8370866009907680753e-12, -3.715007138410239194e-12, 4.9085818713215626784e-12, -6.5418274285256124847e-12, 8.7909550673092636712e-12, -1.1907576877360773869e-11, 1.6252715837504586403e-11, -2.2346851741490000693e-11, 3.0943857501300450716e-11, -4.3140657998500466731e-11, 6.0540496135391454172e-11, -8.5496929381164645249e-11, 1.2147982597501057027e-10, -1.7362635571509051753e-10, 2.4957343771328640668e-10, -3.6071981458611983287e-10, 5.2414604158814657787e-10, -7.6554517591908704089e-10, 1.1237109486301102046e-9, -1.6574333304015960006e-9},
    {1, 1, 1, 0.0012802907736721624407, -0.000023657217983088667576, 1.006950058672728494e-6, -7.8140172644615655313e-8, 9.4572318355051320361e-9, -1.6079564038651651713e-9, 3.5725767612941521752e-10, -9.848797385172279822e-11, 3.2424041764199413357e-11, -1.2383526303122297127e-11, 5.3648500829553087022e-12, -2.5899841246328160082e-12, 1.3736244792875523748e-12, -7.9107864943444661401e-13, 4.899830855454531471e-13, -3.2379386887940382027e-13, 2.2674781739706320864e-13, -1.6730246732912865897e-13, 1.2941918739023719611e-13, -1.0451416766182553393e-13, 8.7784031508812266164e-14, -7.6437105631747652556e-14, 6.8800993928948292649e-14, -6.3853317562564037345e-14, 6.0966186520433523181e-14, -5.9763308221522351283e-14, 6.0039050684287168907e-14, -6.1713259195886634384e-14, 6.4808026856339233492e-14, -6.9439201024855624972e-14, 7.5819120649808664599e-14, -8.4269379675959978798e-14, 9.524406173271674577e-14, -1.0936533080587334369e-13, 1.2747479001574544186e-13, -1.5070589125404539293e-13, 1.8058516770591011855e-13, -2.1917351803414053774e-13, 2.6926367590808939518e-13, -3.346570404456890796e-13, 4.2055322726296397648e-13, -5.3410050711291766462e-13, 6.8517692389774988416e-13, -8.8750359170740057349e-13, 1.1602383391088199248e-12, -1.5302668085799764866e-12, 2.0355103220592674865e-12, -2.7297218712734194892e-12, 3.6894685735434188143e-12, -5.0243389309355289654e-12, 6.8919242159023722702e-12, -9.5198937987334090047e-12, 1.3238649786398458658e-11, -1.8529815143378725303e-11, 2.6098502150269688722e-11, -3.698141984043252589e-11, 5.2709176892385150201e-11, -7.5550811372176955349e-11, 1.0888348523859001352e-10, -1.5775331487788684775e-10, 2.2972900078453373068e-10, -3.3620538093427612971e-10},
    {1, 1, 1, 0.0010969933645898322679, -0.00001800888147175021106, 6.9569283777551428153e-7, -4.9749740553215121841e-8, 5.6114209442504182457e-9, -8.9681769130964231374e-10, 1.8856590778460636791e-10, -4.9462599855418362349e-11, 1.5563713661674106292e-11, -5.7024956365035679913e-12, 2.3775451076518369556e-12, -1.1076446165641720933e-12, 5.6823643628606408205e-13, -3.1720266404836513338e-13, 1.9078656447687624977e-13, -1.2262828729796069336e-13, 8.3647392144201914965e-14, -6.0196036826707437957e-14, 4.5471107970568423199e-14, -3.5896471000413098415e-14, 2.9502556195989044959e-14, -2.5159856953208088685e-14, 2.2198285721814544671e-14, -2.0209811711949138083e-14, 1.8942206374010925485e-14, -1.8240030230025498363e-14, 1.8011134074022773757e-14, -1.8207537393038078672e-14, 1.881482275239035081e-14, -1.9846924845156424711e-14, 2.1344686347390414446e-14, -2.3377423251854142929e-14, 2.6047306099283797744e-14, -2.9496791726588979864e-14, 3.391973641115359953e-14, -3.957726076023854114e-14, 4.6819988349801495741e-14, -5.6119022251139850377e-14, 6.8109057946601858803e-14, -8.3648499352358414877e-14, 1.0390355024010631169e-13, -1.3046629317279738319e-13, 1.6552117884585171822e-13, -2.120807762325183915e-13, 2.7432103705013552274e-13, -3.5806013570243383717e-13, 4.714452944957089817e-13, -6.2594209832606441925e-13, 8.3776546550456187609e-13, -1.129957950689068616e-12, 1.5354204319705812012e-12, -2.1013492110816715548e-12, 2.8957576196227212455e-12, -4.0170983038513552486e-12, 5.6085140879658628889e-12, -7.8790222949280789443e-12, 1.1135120744315037539e-11, -1.5828111313962603937e-11, 2.2625206536567151096e-11, -3.2516740462159055815e-11, 4.6978363226267392371e-11, -6.8217233641091297809e-11},
    {1, 1, 1, 0.00094015198719015528924, -0.000013712502598520912957, 4.8076936403785695302e-7, -3.1682420978574728721e-8, 3.3303787693028875552e-9, -5.0031857552370966152e-10, 9.9553785100131991845e-11, -2.4847574153926867456e-11, 7.4726198049832979692e-12, -2.6266328351827393879e-12, 1.0539351392199556589e-12, -4.7382488362040832555e-13, 2.3512799344530790335e-13, -1.2722377983603114762e-13, 7.4306855161581431963e-14, -4.6454434350584585523e-14, 3.0865707512794525442e-14, -2.1664469269456440734e-14, 1.5980378479658959057e-14, -1.2332270708023968542e-14, 9.9178724992384952565e-15, -8.2837482621409275422e-15, 7.1640562959966285288e-15, -6.3981722885378227631e-15, 5.886904504563320219e-15, -5.5684125033866652457e-15, 5.4045961901622638861e-15, -5.3732730739395130833e-15, 5.4636951662830354729e-15, -5.6740969776192248484e-15, 6.0105726508994896509e-15, -6.486919986192758338e-15, 7.125294118036239724e-15, -7.9576527808564308506e-15, 9.0280871151489301124e-15, -1.0396241486523087011e-14, 1.214215196095338214e-14, -1.4372994666705654336e-14, 1.7232454506522224081e-14, -2.0913731048536741082e-14, 2.5677633272800084723e-14, -3.1877837704359249281e-14, 3.9996282104143541507e-14, -5.0692966756640952969e-14, 6.4876325564124124947e-14, -8.3803088284661804214e-14, 1.0922059763879971571e-13, -1.4357048868085105919e-13, 1.9028140764373034739e-13, -2.5419143177208684683e-13, 3.4216014097133012145e-13, -4.6395904372518975111e-13, 6.335720563433265129e-13, -8.7110224221475507346e-13, 1.2055781071052994525e-12, -1.6790993334579127651e-12, 2.3529838406484925548e-12, -3.3169153047107575189e-12, 4.7026041191301089627e-12, -6.7042600172991527141e-12, 9.609378028343686419e-12, -1.3845189108191753241e-11},
    {1, 1, 1, 0.0008059101882285013648, -0.000010443570135864372262, 3.3232338291086307232e-7, -2.0181438120031066077e-8, 1.9770672754380586084e-9, -2.7918796001839495453e-10, 5.25727117119948571e-11, -1.2485310767885578481e-11, 3.5887324815185773378e-12, -1.2101592154966555874e-12, 4.6731303696045050382e-13, -2.027422737312877794e-13, 9.7317004691508380162e-14, -5.1039803581902458782e-14, 2.8948048386645630603e-14, -1.7602448305597211639e-14, 1.1392250325614469657e-14, -7.7989784191343644447e-15, 5.6175655121873614649e-15, -4.2378340944086162118e-15, 3.3349324686037391048e-15, -2.7280687076222068951e-15, 2.3126412310301378158e-15, -2.0260928638991860115e-15, 1.8300088943530243232e-15, -1.7003841808513881287e-15, 1.6221659357294611228e-15, -1.5861215740037448323e-15, 1.5870209889351480993e-15, -1.6225951871902103215e-15, 1.6929800627844293372e-15, -1.8004886143315332084e-15, 1.9496323920390910828e-15, -2.1473614825569379621e-15, 2.4035269598707551168e-15, -2.7315992439739405278e-15, 3.1497061478236266869e-15, -3.6820904899260517436e-15, 4.3611341326214495271e-15, -5.2301594689894280058e-15, 6.3473092073043751579e-15, -7.7909327545269730139e-15, 9.6670898167016987748e-15, -1.2120044122972655149e-14, 1.5346999252926215224e-14, -1.9618878622135772475e-14, 2.5309752916159761348e-14, -3.2938689574539949849e-14, 4.3229517301446961891e-14, -5.7196528269916598771e-14, 7.6267877452356373209e-14, -1.024639768705506054e-13, 1.3865641828952172285e-13, -1.8894524465403884608e-13, 2.5921077997743213944e-13, -3.5792383771408850351e-13, 4.9733971185352004114e-13, -6.9526417792268502322e-13, 9.7767492756249423678e-13, -1.3826265428258155993e-12, 1.9660876543503455481e-12, -2.8106960622579917688e-12},
    {1, 1, 1, 0.00069098055058124307365, -7.9557146330700786163e-6, 2.2976580425153366723e-7, -1.2858417112932806177e-8, 1.1739556584900990377e-9, -1.5582950710521919167e-10, 2.7769388299392346304e-11, -6.2750665689353277095e-12, 1.7239041995220674279e-12, -5.5768591856144517709e-13, 2.0725550756073852585e-13, -8.6771090398577676973e-14, 4.0288168777414921415e-14, -2.0481142677313570787e-14, 1.1280134160429370294e-14, -6.6715018953682589802e-15, 4.2057891608373896657e-15, -2.8082264911338428414e-15, 1.9752134213683605911e-15, -1.4566314828416466808e-15, 1.1216579707096242191e-15, -8.9864596983290066926e-16, 7.4672808630663857879e-16, -6.4175275436553069121e-16, 5.6901591784345950848e-16, -5.1935909678370031355e-16, 4.8700380942102752688e-16, -4.6831612188212695476e-16, 4.6108816685349486322e-16, -4.6411833823000437032e-16, 4.7697209337593026319e-16, -4.9985882964811820459e-16, 5.3359018659972702679e-16, -5.796028509845849192e-16, 6.4004049060826042744e-16, -7.1789812123097262341e-16, 8.1724003327154449706e-16, -9.4351081478930664268e-16, 1.1039692867180581945e-15, -1.3082887202222530649e-15, 1.5693853091094376507e-15, -1.9045629442313476304e-15, 2.3370993109390790565e-15, -2.8984511479076849743e-15, 3.6313323584174068899e-15, -4.5940280803224337444e-15, 5.8664662724322114772e-15, -7.5587986924077641261e-15, 9.8235791194727719863e-15, -1.287311850223473267e-14, 1.7004319153523373782e-14, -2.2634354891596341526e-14, 3.0352138586964336587e-14, -4.0992854630013094182e-14, 5.5746310883685452703e-14, -7.6315056533231669664e-14, 1.0514598896917033792e-13, -1.4577085747556744078e-13, 2.033086814248049085e-13, -2.8520972044383070466e-13, 4.0236099254227868082e-13, -5.7073470268754467539e-13},
    {1, 1, 1, 0.00059255945100911678239, -6.0618234646075238927e-6, 1.5889347580977875306e-7, -8.1944589555380676981e-9, 6.972363941406894439e-10, -8.6996406285628550206e-11, 1.4671388220709902388e-11, -3.154543704470359996e-12, 8.2829418226015150666e-13, -2.5706115561939992822e-13, 9.1939896934302420594e-14, -3.7145453280175704817e-14, 1.6682697909954033181e-14, -8.2205224449698000268e-15, 4.3965231188879460396e-15, -2.5291485316022554904e-15, 1.5530510002066307153e-15, -1.0114090557560674847e-15, 6.9467272548306023655e-16, -5.0079018959227980995e-16, 3.773412337906679781e-16, -2.9608909110499961969e-16, 2.4116661541392190156e-16, -2.0331838127476378968e-16, 1.7696857928902003176e-16, -1.5866786837378166833e-16, 1.462412932508294588e-16, -1.3830642169490008769e-16, 1.3399416987026715122e-16, -1.3278466152893183243e-16, 1.3441098089829497946e-16, -1.3880497308847694769e-16, 1.4607087139158532054e-16, -1.5647918138522375526e-16, 1.7047732114608601199e-16, -1.8871628273529977105e-16, 2.1209477166417996115e-16, -2.4182436395774804008e-16, 2.7952154950112177982e-16, -3.2733544158693957491e-16, 3.8812381088636738836e-16, -4.6569544526224089241e-16, 5.6514433074029551714e-16, -6.9331176743171470632e-16, 8.5942768101667550743e-16, -1.0760041056473555661e-15, 1.3600850785802804527e-15, -1.7350023736109053538e-15, 2.2328520481423738096e-15, -2.89800219872444538e-15, 3.7920817067895769321e-15, -5.0011040614570599998e-15, 6.6456808080289950828e-15, -8.8957225047293742498e-15, 1.1991681343519283639e-14, -1.6275364200889047274e-14, 2.2234796829112377021e-14, -3.0569789713844335735e-14, 4.2288108717399264548e-14, -5.8847044109716971486e-14, 8.2362546713490160678e-14, -1.1591925473363881747e-13},
    {1, 1, 1, 0.00050825483040137411556, -4.619738933798477312e-6, 1.0990538540420418003e-7, -5.223317502788960173e-9, 4.1419271993235776328e-10, -4.8578875076911378385e-11, 7.7530211091867917676e-12, -1.5861705243991339941e-12, 3.9806262977319088875e-13, -1.1851646930762265112e-13, 4.0794118043156919075e-14, -1.5904938061352235954e-14, 6.9095672568804829054e-15, -3.3002021310164369358e-15, 1.7139589291617788005e-15, -9.5900543997284881498e-16, 5.736143468927465514e-16, -3.6434906423376646758e-16, 2.4436703808183384974e-16, -1.7220990544570250361e-16, 1.269709338123981925e-16, -9.7578119565258918532e-17, 7.7905511076183630377e-17, -6.4429069995706172751e-17, 5.50508785834704214e-17, -4.848491019047387752e-17, 4.392422063004005473e-17, -4.0854696244479981446e-17, 3.8947918286291358442e-17, -3.7998244426699336635e-17, 3.7885496696122045961e-17, -3.8553086599850812442e-17, 3.999593828840577273e-17, -4.2255098084658858292e-17, 4.541740158582042091e-17, -4.96195075995152632e-17, 5.5056274954849386892e-17, -6.1994015487053384236e-17, 7.0789712842637777383e-17, -8.1917940347906857481e-17, 9.6008030543126373301e-17, -1.1389514944084107967e-16, 1.3669044891590229055e-16, -1.658776034795480203e-16, 2.0344605706333586429e-16, -2.5207560619513796291e-16, 3.1539310118034682506e-16, -3.983308925075319181e-16, 5.0762940969861748924e-16, -6.5254473278779116816e-16, 8.4584885394394526371e-16, -1.1052494493359128944e-15, 1.455413205143132765e-15, -1.9308607196103225481e-15, 2.5801246473126493244e-15, -3.4717453863289912138e-15, 4.702949237329863503e-15, -6.4122562413968042703e-15, 8.7978647256688097108e-15, -1.2144556692254797899e-14, 1.6863213957168042938e-14, -2.3549059763325130732e-14},
    {1, 1, 1, 0.00043602494272315526812, -3.5214218704952672025e-6, 7.6036201766862771032e-8, -3.3301375135655016267e-9, 2.4610199144481841804e-10, -2.7132154972588524171e-11, 4.097903575289893164e-12, -7.977272954774189184e-13, 1.9134170620992516235e-13, -5.4652815778817261911e-14, 1.8104344072688486504e-14, -6.8116153532297310508e-15, 2.8623801555758110493e-15, -1.3251762295408473871e-15, 6.6831852793490494399e-16, -3.6371390423112967342e-16, 2.1190755495186994974e-16, -1.3128063591878305136e-16, 8.5979959202148993283e-17, -5.9231497523403358276e-17, 4.2733322820274866728e-17, -3.216435252709941353e-17, 2.5171641154891592462e-17, -2.0421115653901461932e-17, 1.7128713714978845346e-17, -1.481892264662183617e-17, 1.3195642973239429885e-17, -1.207074545090839307e-17, 1.1323353873744323929e-17, -1.0876061411032361508e-17, 1.0680798775450010159e-17, -1.0710402895132553769e-17, 1.0953696382030665628e-17, -1.1412852772212210647e-17, 1.2102372579634158334e-17, -1.3049326622685343302e-17, 1.4294740321532829229e-17, -1.5896154510155997534e-17, 1.7931540859886263519e-17, -2.0504898855090669698e-17, 2.3754038506095978291e-17, -2.7861282031123823007e-17, 3.3068127333225876797e-17, -3.9695344680128866173e-17, 4.8170579199165090697e-17, -5.9066382029431031227e-17, 7.3152802251515398481e-17, -9.1470399507719467712e-17, 1.1543201572777932098e-16, -1.4696521295654842827e-16, 1.8871244063596761009e-16, -2.4431347204503079082e-16, 3.188055265501452181e-16, -4.1919237144231339129e-16, 5.5525694930927017929e-16, -7.4072624043443715153e-16, 9.9494742572715336777e-16, -1.3453088735094233921e-15, 1.8307499122363153913e-15, -2.5068676621649158449e-15, 3.4533744026295248053e-15, -4.7850259671428849865e-15},
    {1, 1, 1, 0.00037412637622678945482, -2.6847371842171321616e-6, 5.2614676390106546014e-8, -2.1235572212973952686e-9, 1.4625626229989584647e-10, -1.5156824291702747968e-11, 2.1664062261534724241e-12, -4.0127919794102037196e-13, 9.1993188196103853951e-14, -2.5207763446347230836e-14, 8.0362991102400907307e-15, -2.9178058972671942005e-15, 1.1860199841999820845e-15, -5.3222471139071512623e-16, 2.6064833461206031917e-16, -1.3797079069278836893e-16, 7.8299923841231883755e-17, -4.7312096281796670302e-17, 3.0258011823213889489e-17, -2.0376794947314632769e-17, 1.4385256898783693362e-17, -1.0604393109430179501e-17, 8.1347375304469658853e-18, -6.4738965150303369698e-18, 5.3305737243559196042e-18, -4.5301793262410945193e-18, 3.9650237634139888916e-18, -3.567097335609229783e-18, 3.2927191430759595984e-18, -3.1136411908511482566e-18, 3.0117803017516946913e-18, -2.9760572616647168832e-18, 3.0005049060777314347e-18, -3.0831748891264788216e-18, 3.2255784545351277691e-18, -3.4325163843379103347e-18, 3.7122275738128732247e-18, -4.0768362544826075707e-18, 4.5431174393469020587e-18, -5.1336368119287344436e-18, 5.8783611099737478685e-18, -6.8168836550039362874e-18, 8.0014732476279712423e-18, -9.5012409991136308287e-18, 1.1407839372748891166e-17, -1.3843275572654049013e-17, 1.6970658551998375336e-17, -2.1009035690073203589e-17, 2.625395556804598305e-17, -3.3106081156010283639e-17, 4.2111166331370584036e-17, -5.4016134432733027492e-17, 6.9848060873940344646e-17, -9.1025857800465763018e-17, 1.1951882331306682071e-16, -1.5807259683816903769e-16, 2.1053241563081502418e-16, -2.8230732557904033002e-16, 3.8103930023446522701e-16, -5.175712692826110771e-16, 7.073526666610387309e-16, -9.7248753570392767602e-16},
    {1, 1, 1, 0.00032106991320065228108, -2.0472242354735894791e-6, 3.6414567709634014252e-8, -1.3544048212067429874e-9, 8.6935492637776051936e-11, -8.4686823448813684754e-12, 1.1455183178743811196e-12, -2.0189382580648678859e-13, 4.423704022841673308e-14, -1.1628952171111381825e-14, 3.567910877645078161e-15, -1.2501075505029443769e-15, 4.9152032096210472222e-16, -2.1379683714955362067e-16, 1.0167432742141729225e-16, -5.2347910839344692148e-17, 2.893751510124591829e-17, -1.7054097685019868729e-17, 1.0650466078556476946e-17, -7.0113902647270389485e-18, 4.8434376518402467826e-18, -3.4968901006966701771e-18, 2.629424719189188615e-18, -2.0527556050103728364e-18, 1.6592368476871214334e-18, -1.3851582769686564421e-18, 1.1916433528508093819e-18, -1.054341029232129766e-18, 9.5767817639588516141e-19, -8.915604290796369712e-19, 8.4943111404134881912e-19, -8.2710775242592278084e-19, 8.220785725463796872e-19, -8.3308144462542725218e-19, 8.5986459338810348085e-19, -9.0307228089455644932e-19, 9.6422480424188564814e-19, -1.0457788522076419754e-18, 1.1512661058304370206e-18, -1.2855176968039319406e-18, 1.454991593894244815e-18, -1.6682306337010047464e-18, 1.9364921698172655899e-18, -2.2746078249220233397e-18, 2.7021557070605963008e-18, -3.2450605794716778783e-18, 3.9377839066240387307e-18, -4.8263312174251454563e-18, 5.9723971382950219884e-18, -7.4591007614191537194e-18, 9.3989531515013976973e-18, -1.1944970132304310532e-17, 1.5306234104435729693e-17, -1.9769772859555185464e-17, 2.5731441040587748233e-17, -3.3739678807420482229e-17, 4.4557756528460414178e-17, -5.9252652158433181142e-17, 7.9322433430163211838e-17, -1.0687950287123532333e-16, 1.4491516662560848609e-16, -1.9768299289625259386e-16},
    {1, 1, 1, 0.0002755830221393395, -1.5613704799211646253e-6, 2.5207053558749466269e-8, -8.6399770988664781737e-10, 5.1684449964079829033e-11, -4.732645479729434308e-12, 6.0582187706251715282e-13, -1.0159688711346013956e-13, 2.1276372427957743236e-14, -5.3657212938184615138e-15, 1.5843577887787750786e-15, -5.3569773008545807951e-16, 2.0373820197067909664e-16, -8.5899196374824706728e-17, 3.9668821691062878501e-17, -1.9865213697743478776e-17, 1.069652865847244258e-17, -6.1484701355600431379e-18, 3.749545592957534957e-18, -2.4129828645206883756e-18, 1.6310664809501174965e-18, -1.1533471089303676236e-18, 8.50080088322830623e-19, -6.5101450796310457857e-19, 5.165647409525970485e-19, -4.2360918954004001523e-19, 3.5820263418212378193e-19, -3.1169457956830336377e-19, 2.7859062833660241149e-19, -2.5533775065803173287e-19, 2.3961559319182530429e-19, -2.2991375105601433041e-19, 2.252757135474165077e-19, -2.2514319831139518245e-19, 2.2926333284061383661e-19, -2.376372362711075713e-19, 2.5049790759117496193e-19, -2.6831104241075579076e-19, 2.9179616958587502471e-19, -3.2196828813261623682e-19, 3.6020259213571889399e-19, -4.0832731730712599655e-19, 4.6875258201290903699e-19, -5.4464669130764495748e-19, 6.4017615872906841829e-19, -7.6083224578280242749e-19, 9.1387590463002892241e-19, -1.1089457357308149581e-18, 1.3588915076128939202e-18, -1.6809211828823934864e-18, 2.0981855033902638206e-18, -2.6419757169130297398e-18, 3.3547838229131467228e-18, -4.2945807583956057308e-18, 5.5408208392036546495e-18, -7.2029019425902811322e-18, 9.4321318920815461175e-18, -1.2438718816645600515e-17, 1.6515985773138847235e-17, -2.2075009867346030006e-17, 2.9694355639996690855e-17, -4.0191738704813969993e-17},
    {1, 1, 1, 0.0002365779665016125856, -1.1910242182695903805e-6, 1.7451982043814255742e-8, -5.5125596235127083591e-10, 3.0732635888335940767e-11, -2.6452674513004781865e-12, 3.2045393246218046459e-13, -5.1134698050568629828e-14, 1.0234984407989198267e-14, -2.4762462844641443522e-15, 7.0367290139955937603e-16, -2.2959932864415472489e-16, 8.446599647235594778e-17, -3.4518779726141973895e-17, 1.5479818118637831808e-17, -7.5399025810623056061e-18, 3.9546051653222902232e-18, -2.2170938194988905809e-18, 1.3202840531570977051e-18, -8.3058307422414603991e-19, 5.4937435251131810107e-19, -3.8046689225534330423e-19, 2.7487659318498370777e-19, -2.0650136215125603468e-19, 1.6084960044356692007e-19, -1.2957171668413160543e-19, 1.0769366074097499273e-19, -9.2162933790744095972e-20, 8.1057336120669887841e-20, -7.3140545743395149455e-20, 6.7605328126501121419e-20, -6.3921469437625053166e-20, 6.1743945412702092027e-20, -6.0856807884921273367e-20, 6.1138969370757369423e-20, -6.2543972065149609928e-20, 6.5089194875532740766e-20, -6.8851949156713938487e-20, 7.397116579279647673e-20, -8.0654230740418834321e-20, 8.9189187195053089397e-20, -9.9963121235739089526e-20, 1.1348818116439085961e-19, -1.3043743543543646354e-19, 1.5169374317925652808e-19, -1.7841610928060459399e-19, 2.1212977473797438952e-19, -2.5484876260326254916e-19, 3.0924305626719511828e-19, -3.7886745215512981614e-19, 4.6847600823309954121e-19, -5.8445577515447775529e-19, 7.3542740969094779017e-19, -9.330801606305216762e-19, 1.1933372500909931725e-18, -1.5379887388042238536e-18, 1.9969882602273364306e-18, -2.6116959043433573215e-18, 3.4394743556601503367e-18, -4.5602273526351326663e-18, 6.0857356162373036477e-18, -8.1730356716279388929e-18},
    {1, 1, 1, 0.00020312467354562629792, -9.0867089705938457882e-7, 1.2084826403490348404e-8, -3.5177729958590190819e-10, 1.8277383415836560461e-11, -1.4788013496649522588e-12, 1.6953568837387112339e-13, -2.5741037468788813761e-14, 4.9243848358948581097e-15, -1.1429701139219602756e-15, 3.1258187702443852808e-16, -9.8423069245027185465e-17, 3.5024093137129169829e-17, -1.3873862089072775771e-17, 6.041684751890534442e-18, -2.8622919212324494781e-18, 1.4623089369830083671e-18, -7.9960747085722551982e-19, 4.6497747562909161157e-19, -2.8594844652390609697e-19, 1.8507209821622953659e-19, -1.2553059390517309298e-19, 8.8897913675600679612e-20, -6.5513542953766220162e-20, 5.0094626501153919692e-20, -3.9639758311560670945e-20, 3.2383774330943437419e-20, -2.7255823359436158803e-20, 2.3588163940992327178e-20, -2.0954502932726418133e-20, 1.9077558013826204861e-20, -1.7774793874945166488e-20, 1.6925846987983677181e-20, -1.6452638110768085814e-20, 1.6307135253287355199e-20, -1.6463890973760364713e-20, 1.6915691669300030788e-20, -1.7671361184991833284e-20, 1.8755187846146409361e-20, -2.0207718955919508147e-20, 2.2087861462148180521e-20, -2.4476385122761115542e-20, 2.7481074008727417543e-20, -3.1243936164128765698e-20, 3.5951081016738337978e-20, -4.1846133586391057475e-20, 4.9248403729307353874e-20, -5.8577507935736135949e-20, 7.0386806554966737626e-20, -8.5408949761940311197e-20, 1.0461813362757904179e-19, -1.2931551451408642102e-19, 1.6124684762617510466e-19, -2.0276513882686766254e-19, 2.5705641297023496052e-19, -3.2845431210421729198e-19, 4.2288017067646358187e-19, -5.4846097473133496076e-19, 7.1640040070356617278e-19, -9.4221118544374605885e-19, 1.2474651778852852193e-18, -1.6622876251351686809e-18},
    {1, 1, 1, 0.00017442764032395357884, -6.9336434540799645292e-7, 8.3696336103791550426e-9, -2.2451919166326404838e-10, 1.08717599634893139e-11, -8.2684104476464935585e-13, 8.9707515677242026539e-14, -1.2960111467608182149e-14, 2.3696775739552888878e-15, -5.2765309947705376226e-16, 1.3887671406472767383e-16, -4.2198395403253270475e-17, 1.4525286045829142111e-17, -5.5771480161285826765e-18, 2.3584309903703094469e-18, -1.0867635433360644428e-18, 5.4081424364676611688e-19, -2.8843146052292552274e-19, 1.6378325415187531918e-19, -9.8461279419705225524e-20, 6.2357200387924914556e-20, -4.1424322018317219667e-20, 2.8755340222281287957e-20, -2.0787985377107965269e-20, 1.5603982909716602612e-20, -1.212899912471952687e-20, 9.7395296645452956191e-21, -8.061865152165526224e-21, 6.8654525342525812515e-21, -6.0044020161596308028e-21, 5.3844070403865816163e-21, -4.943512638507607427e-21, 4.6406595165524086936e-21, -4.4487212290242591467e-21, 4.3502129014856508246e-21, -4.3346375668633174799e-21, 4.3968729382790644758e-21, -4.5362509147007354458e-21, 4.7561296088238230359e-21, -5.0638490025365912348e-21, 5.4710219017998914297e-21, -5.9941564457867442508e-21, 6.6556442257472791123e-21, -7.4851855175777822972e-21, 8.5217655329237937552e-21, -9.8163482211731229569e-21, 1.1435523134884457963e-20, -1.3466434063055767619e-20, 1.602344618121865997e-20, -1.9257186173568130951e-20, 2.3366838055620865066e-20, -2.8616926019009166432e-20, 3.536030710422506863e-20, -4.4069792043352681835e-20, 5.5381800377753257742e-20, -7.0156863472621386737e-20, 8.956380132863298601e-20, -1.1519728532430699699e-19, 1.4924265204407054841e-19, -1.9470782806579973924e-19, 2.5575090576841846667e-19, -3.3814451772264486487e-19},
    {1, 1, 1, 0.00014980626608987507376, -5.2915493233224840099e-7, 5.7974948755571244996e-9, -1.4332031835887061805e-10, 6.467773044567583041e-12, -4.6238485647186138708e-13, 4.7475126271679588022e-14, -6.526212240410651944e-15, 1.1405029836096221656e-15, -2.4363075362843613293e-16, 6.1711364565193247017e-17, -1.8095272196097654169e-17, 6.0249396001576330766e-18, -2.242318129826649964e-18, 9.2078574994169756584e-19, -4.126924656777161443e-19, 2.0004490065728540838e-19, -1.0405880215070243586e-19, 5.7700225504982696798e-20, -3.3908897295309137206e-20, 2.1013710738231920211e-20, -1.3671990122145995691e-20, 9.3028489966895726143e-21, -6.5972699507871505743e-21, 4.8612766946268841153e-21, -3.7118421390892992264e-21, 2.9296725678075848602e-21, -2.3849668400242272072e-21, 1.9985490137676196007e-21, -1.720809420960401801e-21, 1.5199300735844869513e-21, -1.3751097392477074097e-21, 1.2725640288359261857e-21, -1.2031104255449495307e-21, 1.1606840309276459135e-21, -1.1414154993647183575e-21, 1.1430592525722412628e-21, -1.1646483443070123034e-21, 1.2063036159914216927e-21, -1.2691556657163991515e-21, 1.3553577224218833434e-21, -1.4681809086183464487e-21, 1.6121936425634168398e-21, -1.7935360424434146424e-21, 2.0203095922520869891e-21, -2.3031132674717774351e-21, 2.6557711078458393336e-21, -3.0963144329601715597e-21, 3.6483065680927871968e-21, -4.3426318689975092392e-21, 5.2199178907401232253e-21, -6.3338252233075248271e-21, 7.7555316164228503782e-21, -9.5798666855349053973e-21, 1.1933736722411804633e-20, -1.4987738946544568049e-20, 1.8972234166502318191e-20, -2.4199674413617518585e-20, 3.1095737591090049633e-20, -4.0242906418320617498e-20, 5.2441692766924511076e-20, -6.8796968917855915743e-20},
    {1, 1, 1, 0.00012867809466539903219, -4.0389489938325588264e-7, 4.0164267373984808064e-9, -9.1501520557222332941e-11, 3.848366000640537183e-12, -2.5861406722897617507e-13, 2.512873352787035635e-14, -3.2868582803894375327e-15, 5.4899835483215873815e-16, -1.1250796537054165442e-16, 2.7426380310613726472e-17, -7.7607196765561104731e-18, 2.4994730059485851974e-18, -9.0167511474647358001e-19, 3.5955211341920982219e-19, -1.5674218471868423765e-19, 7.4007333718328814386e-20, -3.7547666617178101229e-20, 2.0330754149456075579e-20, -1.1679650459893686692e-20, 7.0825051410075557633e-21, -4.5131117730898640004e-21, 3.010103834892890786e-21, -2.0940361170005267999e-21, 1.5147232480193261351e-21, -1.1361145441939606659e-21, 8.8139033849931980515e-22, -7.0566288115230708694e-22, 5.8187347787324639371e-22, -4.932463974637660619e-22, 4.2911873539165580016e-22, -3.8256674833682817446e-22, 3.4901793232639339127e-22, -3.2541978723961400014e-22, 3.0973171894928289361e-22, -3.0060965406506969413e-22, 2.972088664221158722e-22, -2.9906167035499487467e-22, 3.0600447699753756414e-22, -3.1813924531353325948e-22, 3.358208287867148263e-22, -3.5966592103358975088e-22, 3.905822770079821855e-22, -4.2981926435598897394e-22, 4.7904300535102712318e-22, -5.404417295117783087e-22, 6.1686976985968924933e-22, -7.1204222762437893518e-22, 8.3079710154020008251e-22, -9.7944815726784217916e-22, 1.1662607241149062414e-21, -1.4020949616456161824e-21, 1.7012783657668105875e-21, -2.0827934108887293532e-21, 2.5719001467260223158e-21, -3.2023614321697189774e-21, 4.019506254427617758e-21, -5.084462841531079941e-21, 6.4800304542685360183e-21, -8.3188548563793440882e-21, 1.0754853711670077679e-20, -1.3999242758174671804e-20},
    {1, 1, 1, 0.00011054452969064008141, -3.0833008420859430302e-7, 2.782932341416131162e-9, -5.8426900861489402473e-11, 2.2901417582386374818e-12, -1.4466562903326558015e-13, 1.3302704456956101493e-14, -1.6556395083570814904e-15, 2.6430829451499547741e-16, -5.1963651026685709399e-17, 1.2190938441932535325e-17, -3.3289267295174038214e-18, 1.0370737004233212007e-18, -3.6263390112031612892e-19, 1.4042052899427216594e-19, -5.9540268510015828194e-20, 2.7383414425409514574e-20, -1.355041760252563798e-20, 7.1646517740169934498e-21, -4.0235709771973638804e-21, 2.3874633279521393759e-21, -1.4899995030770930733e-21, 9.7412049229822913335e-22, -6.6476752148963034857e-22, 4.7204344082566799447e-22, -3.4779271476285059208e-22, 2.6520593442119758026e-22, -2.0882282386052218864e-22, 1.6943693765583730369e-22, -1.4140372533930947603e-22, 1.2117056583572635103e-22, -1.0644931747264123164e-22, 9.5737404100000767006e-23, -8.8033554749062542287e-23, 8.2665287158181723871e-23, -7.9182263300498153597e-23, 7.7289520791180635076e-23, -7.6805538273685501225e-23, 7.7636290861724573473e-23, -7.976005645380683659e-23, 8.3219898410526453429e-23, -8.8122106344730445314e-23, 9.4639781481921693447e-23, -1.0302141631839554592e-22, 1.1360487199798120713e-22, -1.2683769268480425626e-22, 1.4330528848859115652e-22, -1.6376923898676496235e-22, 1.8921889920699696329e-22, -2.2094073096937641564e-22, 2.6061147057891471744e-22, -3.1042356585655497482e-22, 3.7325453264001491769e-22, -4.5289636178122751661e-22, 5.5436737545513261393e-22, -6.8433773192444617503e-22, 8.5171217736973804672e-22, -1.0684311725856824527e-21, 1.3505763843949213324e-21, -1.7199019070757414149e-21, 2.2059630727557166811e-21, -2.8490870074690194191e-21},
    {1, 1, 1, 0.000094978652550838763452, -2.3540917811506464842e-7, 1.9285311290364162571e-9, -3.7312919938461818052e-11, 1.3630462627023715252e-12, -8.0935875600256191861e-14, 7.043231536910739667e-15, -8.3409094360934686284e-16, 1.2726631619761049868e-16, -2.4003748277835095583e-17, 5.4196204688531434253e-18, -1.4281362003810830456e-18, 4.3036207600764248414e-19, -1.4586462057847317455e-19, 5.4848233908258202493e-20, -2.262033328212963656e-20, 1.013360068865231996e-20, -4.8908661679054342143e-21, 2.525225008431604476e-21, -1.3862989846149158382e-21, 8.0491488113384216933e-22, -4.9199377870050137126e-22, 3.152879344926997132e-22, -2.1106630233863425826e-22, 1.4712759665839022026e-22, -1.0648350819876083367e-22, 7.9810808618529539835e-23, -6.1804797308270260789e-23, 4.9345910425900177001e-23, -4.054350874547959665e-23, 3.4220027529892056179e-23, -2.9623891598440920883e-23, 2.6265104216458453387e-23, -2.3818592952364526152e-23, 2.2066032992589462966e-23, -2.0860105162858248152e-23, 2.0102175898647558481e-23, -1.9728221087534365243e-23, 1.9699961999278910895e-23, -1.9999414248462406779e-23, 2.0625778784350794516e-23, -2.1594048666273226583e-23, 2.2934990043478153926e-23, -2.4696351127189228236e-23, 2.6945301631013956169e-23, -2.977223473378016173e-23, 3.3296194104278520874e-23, -3.7672336977993832366e-23, 4.3102027980761075546e-23, -4.984639753937532922e-23, 5.8244519110422296423e-23, -6.8737795780709969673e-23, 8.1902746897366217869e-23, -9.8495216213388477259e-23, 1.195101789859477742e-22, -1.4626294037579522542e-22, 1.8049978187997851479e-22, -2.2454929851473350148e-22, 2.8153016698235591312e-22, -3.5563745451379934439e-22, 4.5253862724449948852e-22, -5.7992331349212043344e-22},
    {1, 1, 1, 0.00008161482918594061337, -1.7975827099556500786e-7, 1.3366258469246285657e-9, -2.383227096621180084e-11, 8.1137010676943140703e-13, -4.528738074936278714e-14, 3.7296200492315063067e-15, -4.2026359047166445729e-16, 6.128821817071417159e-17, -1.1089689709372981933e-17, 2.4096923767118621829e-18, -6.1276798922632885999e-19, 1.7861563249545987278e-19, -5.8680346823461676692e-20, 2.1426728569846103834e-20, -8.5950502165190903652e-21, 3.7506039897484716261e-21, -1.7655504814185872867e-21, 8.9015648399612582284e-22, -4.7770902569093611961e-22, 2.7140917522153930103e-22, -1.6247795211915736668e-22, 1.020618358858928646e-22, -6.702385097394896481e-23, 4.5863552746696024505e-23, -3.2606607253789884794e-23, 2.4021580919865471478e-23, -1.8294806510594645551e-23, 1.4373272521200143308e-23, -1.1626345933163318181e-23, 9.6655152775140692158e-24, -8.2452301860545958557e-24, 7.2067268632803549357e-24, -6.4453339835325500955e-24, 5.8909700415186875988e-24, -5.496250882377597928e-24, 5.2291006512582282343e-24, -5.0680956073962240049e-24, 4.9995105276654299768e-24, -5.0154579442495504402e-24, 5.1127555420357754596e-24, -5.2923040257050448682e-24, 5.5588491658348100827e-24, -5.9210615727826890809e-24, 6.3919100872927547753e-24, -6.9893381796308434143e-24, 7.7372832085622664807e-24, -8.6671102090388256104e-24, 9.8195689991644504149e-24, -1.1247429999247988054e-23, 1.30190152526585472e-23, -1.5222923233531424961e-23, 1.7974357917358067766e-23, -2.1423626440702840406e-23, 2.5767582505554577022e-23, -3.1265088477769689453e-23, 3.8257981855998448518e-23, -4.7199609738282113786e-23, 5.8693806970141536985e-23, -7.3548338377351153338e-23, 9.2848444302311771045e-23, -1.1805842512425079508e-22},
    {1, 1, 1, 0.000070139839612951668823, -1.3728102528237569502e-7, 9.2651048777662497756e-10, -1.5224024890791602132e-11, 4.8304281681043369259e-13, -2.5343803511079819276e-14, 1.9752216694469575068e-15, -2.1178191960132679242e-16, 2.9518845638739746182e-17, -5.1241125252463603275e-18, 1.071552197146016539e-18, -2.6295507949758358791e-19, 7.414195553855884452e-20, -2.360992002030668306e-20, 8.3715959850536676124e-21, -3.2663072914409443732e-21, 1.3883464850741528248e-21, -6.3743184076628339784e-22, 3.1382816106775284532e-22, -1.6463769680441518049e-22, 9.1528934060532873276e-23, -5.3664685483811595835e-23, 3.3042944553418265639e-23, -2.1286252366974422955e-23, 1.429883360757033989e-23, -9.98592387553899046e-24, 7.2310414077160065575e-24, -5.4161770410142370145e-24, 4.1871600328955023647e-24, -3.3344527300360482755e-24, 2.7304177010185157401e-24, -2.2952124378632305345e-24, 1.9776818732252832189e-24, -1.744352263649150822e-24, 1.5729278932245626193e-24, -1.4483584513296767675e-24, 1.3604118114636143743e-24, -1.3021503653613162756e-24, 1.2689632832805449585e-24, -1.2579500160602358969e-24, 1.2675327230014439146e-24, -1.2972240083451235967e-24, 1.3475060124036145078e-24, -1.4197956941843241091e-24, 1.5164837468697538337e-24, -1.6410438004671149621e-24, 1.7982161533222148166e-24, -1.9942774607563134274e-24, 2.2374155700694591388e-24, -2.5382379450296751677e-24, 2.9104538802252504286e-24, -3.3717862084048130168e-24, 3.9451890823721149679e-24, -4.6604768832253932693e-24, 5.5565084376712593125e-24, -6.6841248129011728431e-24, 8.110114045031643368e-24, -9.922580800466682609e-24, 1.223824531398423177e-23, -1.5212401292967842302e-23, 1.9052551541700767698e-23, -2.4037148261969000587e-23},
    {1, 1, 1, 0.000060285304237560497102, -1.0485436495732005057e-7, 6.423123697414222651e-10, -9.7263425448900412128e-12, 2.8761306391837181408e-13, -1.4184791728575099771e-14, 1.0462219695516038592e-15, -1.0673647592391503273e-16, 1.4219315980148109551e-17, -2.3679632525644437212e-18, 4.7656501532987008552e-19, -1.1285588486905981207e-19, 3.0779800646666692447e-20, -9.5006559867026044294e-21, 3.2712815170465758848e-21, -1.2414320603333213142e-21, 5.1398664982278749347e-22, -2.3016787581116940637e-22, 1.1065596427454263393e-22, -5.674825320996312346e-23, 3.0870927746983038686e-23, -1.7727200327820694832e-23, 1.0699205596413341012e-23, -6.7612419073199785456e-24, 4.4585226897274035999e-24, -3.0586399043725241191e-24, 2.1769957510504267988e-24, -1.6036713213075490629e-24, 1.2199469566565966827e-24, -9.5645245326487790224e-25, 7.7141960771538694781e-25, -6.3899947796618264968e-25, 5.4279060311284236725e-25, -4.7215053166514959832e-25, 4.200377687906682126e-25, -3.8171840622497323587e-25, 3.5397396965351039134e-25, -3.3460699821547265843e-25, 3.2212777626509877097e-25, -3.1555402977236447216e-25, 3.1428299988892910965e-25, -3.1801143205969972116e-25, 3.2668868977314258315e-25, -3.4049418684686781056e-25, 3.5983419826884390299e-25, -3.8535577131377899585e-25, 4.1797747215393563063e-25, -4.58938429853850123e-25, 5.0986884004602270706e-25, -5.7288698244613755143e-25, 6.5073010451761191087e-25, -7.4692947260953389097e-25, 8.6604379838119972322e-25, -1.013970518544620731e-24, 1.1983615942983531267e-24, -1.4291803718736321892e-24, 1.7194496825330306292e-24, -2.0862602719269261641e-24, 2.5521349724986398226e-24, -3.1468807997279518609e-24, 3.9101126751546827941e-24, -4.8947049081284286036e-24},
    {1, 1, 1, 0.000051821215071259936534, -8.0096805136890024753e-8, 4.4534436049848005137e-10, -6.2147532975749451031e-12, 1.7127192102014476747e-13, -7.9401540816754017245e-15, 5.5422585289539367283e-16, -5.3801193057897116021e-17, 6.8503568969669015933e-18, -1.0944262206350207709e-18, 2.1197579583953405364e-19, -4.8442014594265379198e-20, 1.2779767443458553392e-20, -3.8235615710187541749e-21, 1.2784478828663220993e-21, -4.7189388605406443592e-22, 1.9030987055758851811e-22, -8.3121082161334233369e-23, 3.9022337798547615204e-23, -1.9562811632874276621e-23, 1.041349548296906112e-23, -5.856623501737819139e-24, 3.4648140096395054966e-24, -2.1478767491103049146e-24, 1.3903910871812926515e-24, -9.3696654599542107172e-25, 6.554958691566681246e-25, -4.7489054957913102144e-25, 3.5548231330058418875e-25, -2.7438342635594726345e-25, 2.1797562545681616578e-25, -1.7792369988494604453e-25, 1.489923249941205586e-25, -1.278151894688290973e-25, 1.1218210179533112394e-25, -1.006157229074095262e-25, 9.2114493742000377566e-26, -8.5993300574665853339e-26, 8.178299538564317392e-26, -7.9166200279286414694e-26, 7.7936038751804640206e-26, -7.7969763966588135071e-26, 7.9212409226948747729e-26, -8.1667505054629448677e-26, 8.5393114643547545845e-26, -9.0502235807209979403e-26, 9.7167181305025351321e-26, -1.0562799152640337013e-25, 1.1620532820034716834e-25, -1.2931870090179225546e-25, 1.4551133859050129949e-25, -1.65483586680738344e-25, 1.9013744534670090704e-25, -2.2063584228956761742e-25, 2.5848155365382824977e-25, -3.0562248791035062466e-25, 3.645925203888240026e-25, -4.3870047792202022891e-25, 5.3228459934737992544e-25, -6.5105636850626311261e-25, 8.025667826349328202e-25, -9.9684094717071196468e-25},
    {1, 1, 1, 0.000044550408783837460511, -6.1192041133474262945e-8, 3.0881451003439436512e-10, -3.9714656017046413073e-12, 1.0200386232107629415e-13, -4.4451660684066509358e-15, 2.9363171867561581226e-16, -2.712216203071201515e-17, 3.3006625950538468059e-18, -5.0588468263808372112e-19, 9.4298320461721929811e-20, -2.0795712077313559304e-20, 5.3068126192645266365e-21, -1.5389916792696018584e-21, 4.996913771200562869e-22, -1.7939878143441628963e-22, 7.0473291595390873278e-23, -3.0021435598350365626e-23, 1.3762760355651556467e-23, -6.7447191709062714794e-24, 3.5131542851495228156e-24, -1.9351218252370041019e-24, 1.1221791227790844838e-24, -6.8241115058265572559e-25, 4.3364750833711619063e-25, -2.8706068401242508329e-25, 1.9739504270449959891e-25, -1.406454233689559725e-25, 1.0359742680419275883e-25, -7.8723844137040786514e-26, 6.1599775095234598916e-26, -4.954741942555543605e-26, 4.0902459321431925562e-26, -3.4604964349184175444e-26, 2.9964895009340819956e-26, -2.6524216561246670328e-26, 2.3973894003495287661e-26, -2.2102845558133442554e-26, 2.0765951154726232918e-26, -1.9863686462524400919e-26, 1.932901624356901792e-26, -1.9118934914290954622e-26, 1.9209069929516660593e-26, -1.9590380420378666707e-26, 2.0267364173955311684e-26, -2.1257429646084325932e-26, 2.2591253990799605186e-26, -2.4314067915084695675e-26, 2.6487905694014803692e-26, -2.9194949249753539963e-26, 3.2542190851764215968e-26, -3.6667750808105021976e-26, 4.1749326358383345731e-26, -4.8015430125380290718e-26, 5.5760319374171048169e-26, -6.5363845733489911527e-26, 7.7317902930416182579e-26, -9.2261764628317504867e-26, 1.1102945154045486573e-25, -1.3471343917177180276e-25, 1.6475064530143110353e-25, -2.0303890421863495776e-25},
    {1, 1, 1, 0.000038303842929608376014, -4.6754585882577102293e-8, 2.1416570555603261865e-10, -2.5382166836668134124e-12, 6.0757281834759008927e-14, -2.4888488411662149675e-15, 1.5558608946305081284e-16, -1.3674404540214091484e-17, 1.5905263394082940541e-18, -2.3386672692916871126e-19, 4.1954018221372342278e-20, -8.9284765565175154714e-21, 2.2039237295705474044e-21, -6.1952165856647389759e-22, 1.9533170865013280237e-22, -6.8209791035735622431e-23, 2.609996430707113699e-23, -1.0844358268241445364e-23, 4.8545610030254870945e-24, -2.3256729277186298289e-24, 1.1853595262151648488e-24, -6.3947196726525569927e-25, 3.6349349999764524211e-25, -2.1683784862626280763e-25, 1.3526609823457532728e-25, -8.7958047949647992944e-26, 5.9450398438595137653e-26, -4.1659105009694666988e-26, 3.0194802941832518082e-26, -2.2589518718495279325e-26, 1.7410153275530699077e-26, -1.3799411298237541703e-26, 1.1230193291586610377e-26, -9.3701519867129040715e-27, 8.0048692834376399203e-27, -6.9931296404777965495e-27, 6.2402429404173800196e-27, -5.6817755001673070425e-27, 5.273427227455036147e-27, -4.9846218863793036508e-27, 4.7943915067305620314e-27, -4.6887113256051687469e-27, 4.6587754792410980469e-27, -4.699901784611911001e-27, 4.8108744420567209419e-27, -4.9936089521306587016e-27, 5.253072712027366089e-27, -5.5974287207665914111e-27, 6.0383953891592799038e-27, -6.5918370645733239853e-27, 7.2786207286249691875e-27, -8.125797025142611918e-27, 9.1681908445122775801e-27, -1.0450520925443012204e-26, 1.2030212773500161552e-26, -1.3981129125641458586e-26, 1.6398523269136889983e-26, -1.9405631074837935628e-26, 2.3162469180901991643e-26, -2.7877615499346488453e-26, 3.3824036753727604486e-26, -4.1360427972503511775e-26},
    {1, 1, 1, 0.00003293655737751807636, -3.572741116577701763e-8, 1.4854264018984714009e-10, -1.6223928775061318309e-12, 3.6193433616267207484e-14, -1.3936668509903224335e-15, 8.2449612380674772477e-17, -6.895133801133237719e-18, 7.665329631000880748e-19, -1.081273712727689905e-19, 1.8667816318800483876e-20, -3.8338165203859110949e-21, 9.1539766558745229699e-22, -2.4941764801173904044e-22, 7.636496088725994358e-23, -2.5937282005623532066e-23, 9.6673135790600359401e-24, -3.917660634108961951e-24, 1.7125566609653797827e-24, -8.0201776775474546621e-25, 3.9999408200251794019e-25, -2.1134176173190895566e-25, 1.1775564555318817272e-25, -6.8908801953655073713e-26, 4.2197985110152413117e-26, -2.6954302781779033046e-26, 1.7907046384165083312e-26, -1.2340846060410683763e-26, 8.8016907766746301761e-27, -6.4827361484463927059e-27, 4.9212649967698102151e-27, -3.8437114084168886443e-27, 3.0837257061219483724e-27, -2.5374971873196619064e-27, 2.1386831037736838604e-27, -1.8439592024699621173e-27, 1.6244828590816616704e-27, -1.4607324369558045825e-27, 1.3393214148966720906e-27, -1.2509942658465823698e-27, 1.1893453353619528162e-27, -1.1499899167616703339e-27, 1.1300247741286990565e-27, -1.1276788596769758034e-27, 1.1420931329616423608e-27, -1.1731919075010306128e-27, 1.2216231050839144814e-27, -1.2887547314438423657e-27, 1.3767219259731235946e-27, -1.4885244160271651783e-27, 1.6281790174041498626e-27, -1.800936632938540284e-27, 2.0135785833954939179e-27, -2.2748136342018333356e-27, 2.5958054169413380278e-27, -2.99087091063784143e-27, 3.4784053311050375892e-27, -4.0821086489637648804e-27, 4.8326160430445073142e-27, -5.769671705748707776e-27, 6.945036472631989579e-27, -8.426390247578144096e-27},
    {1, 1, 1, 0.000028324220528154684849, -2.7303953434066692789e-8, 1.0303858522798501519e-10, -1.0371255981337412553e-12, 2.1563015361820640872e-14, -7.8049092462912262853e-16, 4.3697343818383891144e-17, -3.477167707027747571e-18, 3.6946180411650249759e-19, -4.9997887390591653036e-20, 8.3073465391676582259e-21, -1.6463955352918631591e-21, 3.8025240430390233249e-22, -1.0042615236551830806e-22, 2.9858260529175886388e-23, -9.8639580520905184002e-24, 3.5811352765850599421e-24, -1.4154639545031349481e-24, 6.0421151508461614024e-25, -2.7661034507981948172e-25, 1.3499140799666351796e-25, -6.9855119087599938477e-26, 3.8151879574486331611e-26, -2.1900974090280036701e-26, 1.3165689176548145861e-26, -8.2609456083956600257e-27, 5.394389426760424943e-27, -3.6561928051029012301e-27, 2.5659558278950323733e-27, -1.8606251285596966684e-27, 1.3912336779201295947e-27, -1.070755135769248312e-27, 8.4686354068579693484e-28, -6.8724833079963434087e-28, 5.7146261863497569228e-28, -4.862730889424704571e-28, 4.2293923155548221525e-28, -3.7558348194414803025e-28, 3.4019339672007506755e-28, -3.1399854258937335462e-28, 2.9507449295369479363e-28, -2.8208746679058856652e-28, 2.7412798429768058206e-28, -2.7060219943350984588e-28, 2.7116163736444273037e-28, -2.7565940831764455311e-28, 2.8412555614454951228e-28, -2.9675716321583437486e-28, 3.1392084233933660501e-28, -3.3616670956423238709e-28, 3.6425411049955264395e-28, -3.9919044786505074693e-28, 4.4228556957109644173e-28, -4.9522545075629312894e-28, 5.6017047167760032217e-28, -6.3988561125452851686e-28, 7.3791254047690297486e-28, -8.587971722904624212e-28, 1.0083910602796449342e-27, -1.1942516266007919145e-27, 1.4261752188249184243e-27, -1.716909392154867852e-27},
    {1, 1, 1, 0.00002436017480054376191, -2.0868672805938890247e-8, 7.1481686971482072006e-11, -6.6306066389770125119e-13, 1.2848015334896212382e-14, -4.3714323932785197993e-16, 2.3161602559260496458e-17, -1.7537020767366064283e-18, 1.7809657587340771341e-19, -2.31214482617865816e-20, 3.6972477531114121048e-21, -7.0710595173099833683e-22, 1.5797254993561036809e-22, -4.0440264606566257612e-23, 1.1675686305936099608e-23, -3.7516778294399765993e-24, 1.3267320914027899061e-24, -5.1146794424718822546e-25, 2.1319679209826678443e-25, -9.5411446589950911032e-26, 4.556237290886981703e-26, -2.3091854298165380876e-26, 1.2362258116044836943e-26, -6.9614521226162543217e-27, 4.1081204347176891625e-27, -2.5320898766251205041e-27, 1.6252056890385237475e-27, -1.0833304565123348885e-27, 7.4813497800603606036e-28, -5.3408101361678567434e-28, 3.933427411052328445e-28, -2.9831651140978084423e-28, 2.3259419304912032807e-28, -1.8615279429344629208e-28, 1.5271333508689975541e-28, -1.2824986418567164832e-28, 1.1012567033342838059e-28, -9.6580634138420242538e-29, 8.6420088721217369708e-29, -7.882204726627471625e-29, 7.3215517059449287539e-29, -6.920242614988875678e-29, 6.6506876983474574846e-29, -6.4941898804124091858e-29, 6.4387676715811313428e-29, -6.4777526847877577975e-29, 6.6089298832041931056e-29, -6.8340780770819402352e-29, 7.158826916034413317e-29, -7.5927873234936918643e-29, 8.1499428376053106221e-29, -8.8493147986417122278e-29, 9.7159384280261066481e-29, -1.0782212547111470269e-28, 1.2089715775972440441e-28, -1.3691619558563500626e-28, 1.5655876924766802246e-28, -1.8069430209912791033e-28, 2.1043767265036070059e-28, -2.4722272526645661767e-28, 2.9289978407968311087e-28, -3.4986540099679146218e-28},
    {1, 1, 1, 0.000020952908527080935693, -1.5951743110227884624e-8, 4.959461663755046131e-11, -4.2395564785305369263e-13, 7.6561088584201136585e-15, -2.4486422622763533628e-16, 1.2278009265478897096e-17, -8.8456931497264948725e-19, 8.5859338389794674478e-20, -1.0693611787666291331e-20, 1.6456627029272897274e-21, -3.0372524450143771506e-22, 6.5635289477951934665e-23, -1.6286482534342379759e-23, 4.5661113885524162361e-24, -1.4270725305170001493e-24, 4.9157744692665502605e-25, -1.8483501631313449646e-25, 7.5234765737661022841e-26, -3.2913858310042375537e-26, 1.5379875010131946521e-26, -7.6342373028561891307e-27, 4.0061385156838850432e-27, -2.2130057616385920841e-27, 1.2820029939056702798e-27, -7.7620204411346451155e-28, 4.8968938571539950698e-28, -3.2102521176277719407e-28, 2.1815093328884124214e-28, -1.5332103183888049586e-28, 1.1122144675466272458e-28, -8.3121000982463063735e-29, 6.3889671054471812008e-29, -5.0428004440061189193e-29, 4.081430988248382748e-29, -3.3828283142621869209e-29, 2.8677779873533779667e-29, -2.4838192791969571071e-29, 2.1955832651605758237e-29, -1.978855891144759438e-29, 1.8168579775140265e-29, -1.6978731298710617098e-29, 1.613712357593476785e-29, -1.55870896588391029e-29, 1.529056928874305835e-29, -1.5223773078818703234e-29, 1.5374407490163282635e-29, -1.5740012087461106114e-29, 1.632713470360140282e-29, -1.7151186589678537687e-29, 1.8236902250709089217e-29, -1.9619392555031781634e-29, 2.1345834996244267626e-29, -2.3477899178222004355e-29, 2.6095065016091388575e-29, -2.9299062006878952851e-29, 3.321974717203245492e-29, -3.8022855448762468411e-29, 4.3920210505573505913e-29, -5.1183191059394232987e-29, 6.0160528074321766203e-29, -7.1301889686027674203e-29},
    {1, 1, 1, 0.000018023892150391096757, -1.2194507165830804372e-8, 3.4412629227016204846e-11, -2.711012725248543428e-13, 4.5627257147604642315e-15, -1.3717383624216993922e-16, 6.5092625885275345661e-18, -4.4622355595241136653e-19, 4.1396549347421419688e-20, -4.9462775658849302006e-21, 7.3256781252848999794e-22, -1.3047342352012382597e-22, 2.7273318400437482295e-23, -6.5597220060071903862e-24, 1.7858930039721934218e-24, -5.4288948130376637619e-25, 1.821568694385314256e-25, -6.680284903421883559e-26, 2.6552254305495317457e-26, -1.1355389650474908046e-26, 5.1921140146058834965e-27, -2.5241633440584452845e-27, 1.2983718123255573729e-27, -7.0357470989938368898e-28, 4.0011044741682731267e-28, -2.3796628707736838253e-28, 1.4756319207124508842e-28, -9.51398209500861504e-29, 6.3617880649420010133e-29, -4.4019118309126996913e-29, 3.1452194837475556106e-29, -2.3162704123394209177e-29, 1.7551226181991799897e-29, -1.3662149587406824863e-29, 1.0909202125018924337e-29, -8.9237632338775856325e-30, 7.4687432700309227222e-30, -6.3884423622112051836e-30, 5.5786634914901333263e-30, -4.9685043911942446811e-30, 4.5090378601418525091e-30, -4.1661433414961627717e-30, 3.9158924520978458473e-30, -3.7415378922800187764e-30, 3.6315298745964409594e-30, -3.5782057473412606779e-30, 3.5769321435263167697e-30, -3.6255613589563982935e-30, 3.724115690767698847e-30, -3.8746473793073278102e-30, 4.0812449052811196681e-30, -4.3501732095399485061e-30, 4.6901488304158025956e-30, -5.1127630003886620457e-30, 5.6330779505348988717e-30, -6.2704353956708978915e-30, 7.0495328011430838565e-30, -8.0018441598894395607e-30, 9.1674896298151264231e-30, -1.0597695110163757658e-29, 1.2358032188536165908e-29, -1.4532695670854404046e-29},
    {1, 1, 1, 0.000015505725761581584295, -9.323134870329400344e-9, 2.3880506561091940087e-11, -1.7337458310982272057e-13, 2.7194657718414042313e-15, -7.6852916595954033569e-17, 3.451269411156280112e-18, -2.2512122008877657231e-19, 1.9961084264940656202e-20, -2.288105196307825071e-21, 3.2613568065119819028e-22, -5.605401648797254736e-23, 1.1333971105864954777e-23, -2.6423305018305598223e-24, 6.9856682742719511818e-25, -2.0654772478553285847e-25, 6.7506063379202980133e-26, -2.4146231465224704284e-26, 9.3719050141523286102e-27, -3.918040789921746263e-27, 1.7529894890723401109e-27, -8.346664222455997977e-28, 4.208389142575583981e-28, -2.2370804391448971492e-28, 1.2488620133180543634e-28, -7.2962516778424117523e-29, 4.447122759062554727e-29, -2.8198710603342115554e-29, 1.8554320603330843038e-29, -1.2639348504946916711e-29, 8.8952287892439704952e-30, -6.4552268295293017995e-30, 4.8220091610750628495e-30, -3.7017752828003992968e-30, 2.9161998587241231555e-30, -2.3542896320809803512e-30, 1.9453301070059824296e-30, -1.6432882267652900469e-30, 1.4176016093551518611e-30, -1.2476160977260495063e-30, 1.1191558478676229376e-30, -1.0223673434214285712e-30, 9.5034036609043432222e-31, -8.9821241741698906592e-31, 8.6257999582374970875e-31, -8.4110867671906937533e-31, 8.3227493291040281709e-31, -8.3519764057893990496e-31, 8.4953280792040418169e-31, -8.7541503841417743438e-31, 9.1343584388272539218e-31, -9.6465351432063821307e-31, 1.0306326837810446947e-30, -1.11351456516705531e-30, 1.2161214606781610487e-30, -1.3421019073850773861e-30, 1.4961259819457030313e-30, -1.6841441749673365986e-30, 1.9137282201314668977e-30, -2.1945187898363939748e-30, 2.5388136599367992144e-30, -2.9623416273469182578e-30},
    {1, 1, 1, 0.000013340552800278257296, -7.1285328411443039676e-9, 1.657335430959645889e-11, -1.1088706023368440529e-13, 1.6210059868969630089e-15, -4.3061718511444438484e-17, 1.8300711428062034945e-18, -1.1358539540834401198e-19, 9.6260097104675558556e-21, -1.0585605906699624788e-21, 1.4520817534800711883e-22, -2.4084282921152764373e-23, 4.7105167524147560282e-24, -1.0644645406168522328e-24, 2.7327686064318341268e-25, -7.8590816521609256361e-26, 2.5019718863215605704e-26, -8.7286315631896013186e-27, 3.3082384041026589064e-27, -1.352005062320481016e-27, 5.919116216441977193e-28, -2.7602657895588460202e-28, 1.3641910205374819955e-28, -7.1136987627350036324e-29, 3.8984459532237635647e-29, -2.237312670986741909e-29, 1.3403638471202937494e-29, -8.3586990278633291177e-30, 5.4119457961181725227e-30, -3.6295314620167482563e-30, 2.5159720969590784432e-30, -1.7991869480872663145e-30, 1.324924527606258838e-30, -1.0030985695996371532e-30, 7.7962207781689856551e-31, -6.2117546407022623217e-31, 5.0673586184124329038e-31, -4.2274167636751119471e-31, 3.6026393758829453527e-31, -3.1331329640312296238e-31, 2.7780484900023175755e-31, -2.5091250074273810963e-31, 2.306588750224874859e-31, -2.1565053310207478285e-31, 2.0490460552877329451e-31, -1.977340657150762133e-31, 1.9367140734974198031e-31, -1.924180641486007832e-31, 1.9381158435121918868e-31, -1.9780552040310264964e-31, 2.0445890533007452355e-31, -2.1393347336498792294e-31, 2.2649769633422876711e-31, -2.4253741021255608138e-31, 2.625734066088866774e-31, -2.8728694149458789634e-31, 3.1755473448423850159e-31, -3.5449575995294635288e-31, 3.9953303573726718098e-31, -4.5447477996032080664e-31, 5.2162083910712015426e-31, -6.0390233343920243892e-31},
    {1, 1, 1, 0.000011478701359944563797, -5.45101942257891127e-9, 1.1503165832286579103e-11, -7.0927842245518147342e-14, 9.6633184302848841112e-16, -2.4130324016512917567e-17, 9.7050544650088337351e-19, -5.7315161973756745626e-20, 4.6424741225160315847e-21, -4.8977490815898460097e-22, 6.4658397430818048953e-23, -1.0349082547122573486e-23, 1.9579255733933653553e-24, -4.2886090511752102108e-25, 1.0691508360650577887e-25, -2.9906419637155071964e-26, 9.2739194651307217963e-27, -3.1556165944215740805e-27, 1.1679036028632613434e-27, -4.6658306576726280523e-28, 1.998829369840669902e-28, -9.1291464159254621065e-29, 4.4225810378993232189e-29, -2.2623024527933120303e-29, 1.2170541799310829795e-29, -6.8611181443104223008e-30, 4.0402444111339875148e-30, -2.477932525587624667e-30, 1.5787129674908466284e-30, -1.042360158437552619e-30, 7.11698226720137004e-31, -5.0151329039649096292e-31, 3.6407899459849981497e-31, -2.7184324708453719345e-31, 2.0844541924029472413e-31, -1.6391174781684636406e-31, 1.320113769888622188e-31, -1.0876214849368799106e-31, 9.1564848153492373463e-32, -7.8689734324334669195e-32, 6.8965275453333777021e-32, -6.1585579327914603453e-32, 5.5988981471811773381e-32, -5.1780164308192207981e-32, 4.867942800894231997e-32, -4.6489223283761339163e-32, 4.5071873986603191075e-32, -4.4334706543046405814e-32, 4.4220202272977899109e-32, -4.4699661172134672551e-32, 4.5769422528706835211e-32, -4.7449053702358540579e-32, 4.9781169009693773658e-32, -5.2832722956470589536e-32, 5.6697766517767304781e-32, -6.1501782852930157015e-32, 6.7407845771769906183e-32, -7.4624984309016072845e-32, 8.341930374479103665e-32, -9.4128622918027552041e-32, 1.0718165886417854819e-31, -1.2312314744548079624e-31},
    {1, 1, 1, 9.8775201822370731875e-6, -4.1686328539623102194e-9, 7.9847874977556359201e-12, -4.5372423925541324484e-14, 5.761127967064537909e-16, -1.3523049720672020418e-17, 5.1471612492853086943e-19, -2.8923871274790578552e-20, 2.2391985844131370097e-21, -2.266299541219132065e-22, 2.8793788513569985287e-23, -4.4474391920412788053e-24, 8.1388654467755800487e-25, -1.7279922998894071198e-25, 4.183263347938027475e-26, -1.1381438804132590998e-26, 3.4378296489993463595e-27, -1.1409390679471636962e-27, 4.1234179844434501152e-28, -1.6103482333833311102e-28, 6.7504816515147512397e-29, -3.0196015172044571875e-29, 1.4338925668749739764e-29, -7.1952530220888571991e-30, 3.7998681029334108815e-30, -2.104279065816670689e-30, 1.2179593632187844185e-30, -7.3465002489734708219e-31, 4.6056735532285392051e-31, -2.9938172316909669044e-31, 2.0133820585477161655e-31, -1.398069793364768591e-31, 1.0005536060341867097e-31, -7.367730740394011106e-32, 5.573664767956969503e-32, -4.3255977927300316328e-32, 3.4393894151689658353e-32, -2.7984706887796974048e-32, 2.3274322186364660931e-32, -1.9765035974687450675e-32, 1.7122272630323031149e-32, -1.5117363129147261352e-32, 1.359174596275359548e-32, -1.2434164213341683889e-32, 1.1565901826415475132e-32, -1.0931087601999110239e-32, 1.0490254261646040769e-32, -1.0216029289670100319e-32, 1.0090252022880759279e-32, -1.0102069499715208818e-32, 1.0246726251289708502e-32, -1.0524868283716171718e-32, 1.0942251557635341586e-32, -1.1509794155415749119e-32, 1.224394777940072171e-32, -1.3167394025320309389e-32, 1.4310098126886906606e-32, -1.5710780807204414044e-32, 1.741890041871034788e-32, -1.9497275725797455788e-32, 2.2025527977171252787e-32, -2.51045836868588053e-32},
    {1, 1, 1, 8.5003812295671840538e-6, -3.1882098871591100624e-9, 5.5430313264062192108e-12, -2.902722566285367691e-14, 3.4350041367081535452e-16, -7.5792246287857810283e-18, 2.7300858702668015247e-19, -1.459762230956943068e-20, 1.0801264265458193243e-21, -1.0487621260328629973e-22, 1.2823649107215171255e-23, -1.9114245422113094652e-24, 3.3835340784410907308e-25, -6.9631566591843868129e-26, 1.636931356945788378e-26, -4.3318056742546152479e-27, 1.2745135394063595696e-27, -4.125530029586218196e-28, 1.4559511905893283892e-28, -5.5583991823638097765e-29, 2.2799898441052643247e-29, -9.9886827154918006756e-30, 4.6493961969176840491e-30, -2.2886564762187836176e-30, 1.18649598898702104e-30, -6.454326340724379003e-31, 3.6719529932768328733e-31, -2.1782648246822054654e-31, 1.3437618048140283084e-31, -8.5994748978623246288e-32, 5.696337042519138575e-32, -3.8977540855498335734e-32, 2.7499467273697581824e-32, -1.9970462686151583764e-32, 1.4904881257630678953e-32, -1.1416195028814053799e-32, 8.9617027221426531369e-33, -7.2011680082175234919e-33, 5.9164950247677066987e-33, -4.964966711722665124e-33, 4.2513958546279432793e-33, -3.7111821090159923401e-33, 3.2997963597045078645e-33, -2.9861318785512918563e-33, 2.7482279327960118948e-33, -2.5704767139245093071e-33, 2.4417749301840027969e-33, -2.3542882705484634723e-33, 2.3026213255342021368e-33, -2.2832617142618700362e-33, 2.2942146986126929689e-33, -2.3347748581120248598e-33, 2.4054012079634777362e-33, -2.5076755540463345328e-33, 2.6443334267835865057e-33, -2.8193642462529856491e-33, 3.0381835417419416264e-33, -3.3078858831309994363e-33, 3.637593336786930399e-33, -4.0389213525383422044e-33, 4.5265926816206216519e-33, -5.1192410004559799407e-33},
    {1, 1, 1, 7.3158248209463342579e-6, -2.4385773588880032215e-9, 3.8482945406191464472e-12, -1.8571902409954475153e-14, 2.0482573011891825442e-16, -4.2482742337504470728e-18, 1.4481800549784562579e-19, -7.3679320942740732764e-21, 5.2106801367534744664e-22, -4.8537182612683543153e-23, 5.7116601678560422436e-24, -8.2156559586576046625e-25, 1.4067444581885301619e-25, -2.8061347505674407561e-26, 6.4059529028762908119e-27, -1.6488410475233888197e-27, 4.7254437141436973469e-28, -1.4918841734076866146e-28, 5.1413164422467605886e-29, -1.9187471096130994187e-29, 7.7013905357644252251e-30, -3.3044934310311540519e-30, 1.5076988287896638782e-30, -7.2803667311660086353e-31, 3.7051188093511410852e-31, -1.9798697657887291577e-31, 1.1071323750200705378e-31, -6.4592037991461354413e-32, 3.9209343051810619719e-32, -2.4703399391490589542e-32, 1.6117708720938581132e-32, -1.086771307245845199e-32, 7.5586867709130631134e-33, -5.4135314028000924754e-33, 3.9861565625886832501e-33, -3.0132473443316095442e-33, 2.335274181201055666e-33, -1.8532042541910520052e-33, 1.5041464810151301673e-33, -1.2473066622952259547e-33, 1.0556983749837100556e-33, -9.1114324569693149874e-34, 8.0119310879308316146e-34, -7.1719878774177566955e-34, 6.5307668097282325065e-34, -6.0450814190635390066e-34, 5.6841225080221586342e-34, -5.4259443269642822219e-34, 5.2551028223370513788e-34, -5.1610637676726714203e-34, 5.1371372794310957913e-34, -5.1797826616642612228e-34, 5.2881838693683805397e-34, -5.4640332020144327341e-34, 5.7114865210894221944e-34, -6.0372719675278034031e-34, 6.4509488626839854237e-34, -6.9653263016448530461e-34, 7.5970634692245751278e-34, -8.3674872714696777603e-34, 9.3036787909425950929e-34, -1.0439899790170311617e-33},
    {1, 1, 1, 6.2968268087812172843e-6, -1.8653555292251253568e-9, 2.6719315503981762487e-12, -1.1883479231150904851e-14, 1.2214573183632082106e-16, -2.3814253210137288952e-18, 7.6825529495318725335e-20, -3.719169034643965932e-21, 2.5139179482779221805e-22, -2.2465135085305786407e-23, 2.5441928595595110532e-24, -3.5315413359878979235e-25, 5.8492052023877286712e-26, -1.130961665731272298e-26, 2.5071137506780696276e-27, -6.2766173756707513803e-28, 1.7521762468518017318e-28, -5.3954481798589910171e-29, 1.8156784110353220589e-29, -6.6240382184564725449e-30, 2.6016116893269716442e-30, -1.0932982912572448874e-30, 4.8895596420413567479e-31, -2.3161304587078715794e-31, 1.1571112735668733095e-31, -6.0737855571460017567e-32, 3.3384057605520374587e-32, -1.9155100463214579663e-32, 1.1441790224413268446e-32, -7.0970628755067542951e-33, 4.5608740387527597068e-33, -3.0303936398388849292e-33, 2.0778087654691772847e-33, -1.4676089593591171808e-33, 1.0661476346317668154e-33, -7.9539957028729854301e-34, 6.0858660802746012896e-34, -4.7695873344374502209e-34, 3.8243086402623021159e-34, -3.1337714509703253198e-34, 2.6217141634370396688e-34, -2.2371661287629620323e-34, 1.9454696245269874598e-34, -1.7226906429481736564e-34, 1.5520753134532622306e-34, -1.4217650164242029387e-34, 1.3233003243348066586e-34, -1.250628207317301619e-34, 1.1994360750711830761e-34, -1.1667020221658403485e-34, 1.1503910163614743789e-34, -1.1492519914434358095e-34, 1.1626868713465561773e-34, -1.1906730216243173803e-34, 1.2337276528716328108e-34, -1.2929076220560796608e-34, 1.3698417409774664262e-34, -1.4667956702235969389e-34, 1.5867721517808041573e-34, -1.7336520260200951922e-34, 1.9123844654938981693e-34, -2.1292384185980201726e-34},
    {1, 1, 1, 5.42017025359118745e-6, -1.4269910750295924693e-9, 1.8553142933139939491e-12, -7.6044228985477913683e-15, 7.2846331648703181501e-17, -1.3350487186607864626e-18, 4.0759079737205584857e-20, -1.8775092714719692419e-21, 1.2129522435319127904e-22, -1.0398709919896204098e-23, 1.1333751578318897901e-24, -1.518176835235336398e-25, 2.432285449687781382e-26, -4.5585141390404897845e-27, 9.8129686644181191741e-28, -2.3895086625690578971e-28, 6.49754185746106418e-29, -1.9514438497596399326e-29, 6.412681193907961577e-30, -2.2869888899036994611e-30, 8.789252604954416042e-31, -3.6175004699233936913e-31, 1.5858461632484287266e-31, -7.3690065593912053273e-32, 3.6139678158365711661e-32, -1.8634530728569774005e-32, 1.0067343140313377592e-32, -5.6810166927508547019e-33, 3.3391396719933822892e-33, -2.0390918033572861712e-33, 1.2907111157851243815e-33, -8.4507670128402563058e-34, 5.7121684109901813993e-34, -3.9790212201518295811e-34, 2.8517835073016547453e-34, -2.0997719713052506183e-34, 1.5861458544252940781e-34, -1.2276499209903674798e-34, 9.7241569001023789129e-35, -7.874040059420408756e-35, 6.5112902062535851366e-35, -5.493460894923784105e-35, 4.7244138299044196144e-35, -4.1381980645532279248e-35, 3.6889064175758813648e-35, -3.3441806308223768558e-35, 3.0809854256355336956e-35, -2.8828187848705813258e-35, 2.7378473766960913928e-35, -2.6376484782943177832e-35, 2.5763569467599468659e-35, -2.5500883968922180327e-35, 2.556555552755542511e-35, -2.5948242176642912232e-35, 2.6651747409158871993e-35, -2.7690480933986851373e-35, 2.909065104478088096e-35, -3.0891146138463668802e-35, 3.3145123071963267971e-35, -3.5922376006020554016e-35, 3.9312617357357184192e-35, -4.3429868175650838205e-35},
    {1, 1, 1, 4.6659065960284458265e-6, -1.0917283613173598182e-9, 1.2883798157216884175e-12, -4.8665753495327134262e-15, 4.3448199953403225857e-17, -7.4850048391966724977e-19, 2.1626093797230334193e-20, -9.478797133488769655e-22, 5.8529028829711749171e-23, -4.8137656374752927058e-24, 5.0493145936008518375e-25, -6.5270288716522943124e-26, 1.0115034173118830809e-26, -1.8375276157591948284e-27, 3.8411560257955180044e-28, -9.0975979703148710646e-29, 2.4096589155510717604e-29, -7.0586198677578645869e-30, 2.2650384571806631374e-30, -7.8966065519667261055e-31, 2.9695910475679599108e-31, -1.1970540441548253096e-31, 5.1438418262835807316e-32, -2.3447152944184046352e-32, 1.1288303827465254088e-32, -5.7175862570966929599e-33, 3.0361686088288259637e-33, -1.6850118711306952326e-33, 9.7456421870568472563e-34, -5.859090163027297413e-34, 3.6529632435257775657e-34, -2.3568312842733786351e-34, 1.5704775006648140639e-34, -1.0788906906220767778e-34, 7.6287090333757164134e-35, -5.5436296847173495873e-35, 4.1342729618257207415e-35, -3.1601200042308134861e-35, 2.4727845745768599957e-35, -1.9786237816436169627e-35, 1.6172757351371594488e-35, -1.3490534219950349711e-35, 1.1473785337675408576e-35, -9.9414690313638931556e-36, 8.7683484174421590507e-36, -7.8665977950359390723e-36, 7.1739143303074398735e-36, -6.6457161149571773668e-36, 6.2499520023832389808e-36, -5.9636095761623338213e-36, 5.7703466974415155558e-36, -5.6588805860190315566e-36, 5.6218988199652329637e-36, -5.655339717972565867e-36, 5.7579435890611258245e-36, -5.9310123792642198807e-36, 6.1783401646467726949e-36, -6.5062950745281670945e-36, 6.9240474255237155284e-36, -7.4439510761033146082e-36, 8.0820967506150658317e-36, -8.8590685710062675755e-36},
    {1, 1, 1, 4.0168934959386761766e-6, -8.3529676467786593852e-10, 8.9475420445661189935e-13, -3.1146864137708625692e-15, 2.5916105488591963032e-17, -4.1968264893705422341e-19, 1.1475347089729777932e-20, -4.7858431472282098849e-22, 2.8244446814048716182e-23, -2.2285613654232933229e-24, 2.2497039774350858673e-25, -2.8063571288846252493e-26, 4.206824729058806196e-27, -7.4076202884286070135e-28, 1.5036881231786929561e-28, -3.4640102737015957533e-29, 8.9370949213314598173e-30, -2.5533940968160384641e-30, 8.0010289430876476131e-31, -2.7267874355766777809e-31, 1.0034036841910905356e-31, -3.9614417507111740782e-32, 1.6685856945241787208e-32, -7.4611480228875261304e-33, 3.5262043064149448546e-33, -1.7544514998428236337e-33, 9.1573809562437023344e-34, -4.9982071382955333538e-34, 2.8445968771683349093e-34, -1.6836738914964458555e-34, 1.0339414630214067575e-34, -6.5734787125260615004e-35, 4.3181411899398548393e-35, -2.9255871707858093566e-35, 2.0408916927740792509e-35, -1.4636953765470837829e-35, 1.0776794172808944389e-35, -8.1351773793852302009e-36, 6.2886152882224135463e-36, -4.9723678584614400984e-36, 4.0173115390975132908e-36, -3.3131926056889836092e-36, 2.7867623566655009039e-36, -2.3884947071044367579e-36, 2.084358363738014922e-36, -1.8506257326772575757e-36, 1.6705410867639332617e-36, -1.5321478022094065841e-36, 1.4268508093323814265e-36, -1.348453354702981998e-36, 1.2925049665622772834e-36, -1.2558572116111570221e-36, 1.2363608863457829961e-36, -1.2326616828368206789e-36, 1.2440664158839736669e-36, -1.2704617865594207846e-36, 1.312274338633145641e-36, -1.3704649588060485598e-36, 1.4465547616019694284e-36, -1.5426820072346539883e-36, 1.6616921954438507933e-36, -1.8072659462995072262e-36},
    {1, 1, 1, 3.4583983615864900357e-6, -6.3914453470203928295e-10, 6.2143567514657836989e-13, -1.9936001872189086129e-15, 1.5459687397411281491e-17, -2.3533312775344920716e-19, 6.0895721024426792013e-21, -2.4165564694337583991e-22, 1.3631012852873445114e-23, -1.0318048999499540811e-24, 1.0024245284123499184e-25, -1.206712265799488019e-26, 1.7497455136086606344e-27, -2.9864618701499768595e-28, 5.8869050294506840899e-29, -1.3190615279760363754e-29, 3.3149015286499215775e-30, -9.2373921278133028397e-31, 2.8265038802506069206e-31, -9.4166311647678090836e-32, 3.3906911783771548199e-32, -1.3110712574965192769e-32, 5.4130612136977990079e-33, -2.3744045385172655027e-33, 1.1015894702890917527e-33, -5.3839812259023421765e-34, 2.7621687654379301269e-34, -1.4827197503298592404e-34, 8.3035639776867790045e-35, -4.8385955052071717127e-35, 2.9267130805113084503e-35, -1.8335619084118183029e-35, 1.1873958081411502721e-35, -7.9338173211110857696e-36, 5.4603744291078058475e-36, -3.86492216337731446e-36, 2.8094000991270919005e-36, -2.0944212934736304435e-36, 1.5994009137635032998e-36, -1.2496743235509467047e-36, 9.979769854596747019e-37, -8.1376272554659880307e-37, 6.7690344816135297946e-37, -5.7389386220841868801e-37, 4.9551928252587496255e-37, -4.3539539156423348072e-37, 3.8903772691646008296e-37, -3.5325890873613244323e-37, 3.2577221037175454858e-37, -3.0492725121053420447e-37, 2.895317004752689582e-37, -2.7872991775718275969e-37, 2.719199545562659229e-37, -2.6869691793814361109e-37, 2.6881488843590241522e-37, -2.7216230586156271667e-37, 2.7874754471949058782e-37, -2.8869263896693546064e-37, 3.0223400089412279765e-37, -3.197296519897542088e-37, 3.4167304066103789043e-37, -3.6871403218821541323e-37},
    {1, 1, 1, 2.9777581753985254235e-6, -4.8908986746332078156e-10, 4.316387512913088825e-13, -1.2761269165510623457e-15, 9.222823452677213708e-18, -1.3197067535961561545e-19, 3.2317675502091818448e-21, -1.2203035702916792727e-22, 6.5789361308434729742e-24, -4.7775264190931793289e-25, 4.4669443706228725613e-26, -5.1891603440471596714e-27, 7.2782664729708996616e-28, -1.2041146645825646249e-28, 2.304883154456347071e-29, -5.0232358895136400317e-30, 1.2296387511457269343e-30, -3.3420548059472624074e-31, 9.9858720555835599639e-32, -3.2521645444606225421e-32, 1.1458650242230902393e-32, -4.339423206338333666e-33, 1.7561841589304666729e-33, -7.5567746519631366313e-34, 3.4416348101461075902e-34, -1.6523361954512461025e-34, 8.3322415674536439551e-35, -4.3988242734372177577e-35, 2.4240471633838339773e-35, -1.3906357291896486539e-35, 8.2850870185404663219e-36, -5.1147996534353732588e-36, 3.2653288027688514785e-36, -2.1517117301969997848e-36, 1.4610249430062071233e-36, -1.0206187563953738657e-36, 7.3243711472703254952e-37, -5.3925453891945129109e-37, 4.0681073183623613107e-37, -3.1409656986704771718e-37, 2.4793525908757252512e-37, -1.9988566858012782539e-37, 1.6443197510311445114e-37, -1.3790233964882379564e-37, 1.1780981888292389654e-37, -1.0244288956316182562e-37, 9.0606434104491737218e-38, -8.1455112331381503685e-38, 7.4384469144128285111e-38, -6.8958739680504685207e-38, 6.4862366329250761782e-38, -6.1867088423869852813e-38, 5.9809430925721727745e-38, -5.8575262586481228348e-38, 5.8089259801862878767e-38, -5.8307862018875672387e-38, 5.9214796147426237845e-38, -6.0818577720892905048e-38, 6.3151625933514794529e-38, -6.6270796723599365587e-38, 7.0259268216712678255e-38, -7.5229823721127772619e-38},
    {1, 1, 1, 2.5640875754529165281e-6, -3.7429051028044901939e-10, 2.9983044423678431115e-13, -8.1692271523861030058e-16, 5.5024807175057087566e-18, -7.401220725032890818e-20, 1.7152406857033108375e-21, -6.1626922836110636209e-23, 3.1755205037767858607e-24, -2.2122813054262576562e-25, 1.9906787042346123895e-26, -2.2316302696236226009e-27, 3.0276998815152924297e-28, -4.8552382013879440823e-29, 9.0249049504058940087e-30, -1.9130835747276208856e-30, 4.5615906738959000163e-31, -1.2092319863012626339e-31, 3.5282088783648037359e-32, -1.1232627294804601896e-32, 3.8726708367561270721e-33, -1.4363808906868507161e-33, 5.6980870432847599503e-34, -2.4051941628454147244e-34, 1.075329876666379704e-34, -5.0713690689500336367e-35, 2.5136535541263503882e-35, -1.3051068908175163562e-35, 7.0770054116757654168e-36, -3.9970482014090002304e-36, 2.3455567195424160692e-36, -1.4269002650851316943e-36, 8.9802881535734181767e-37, -5.8360355185775660065e-37, 3.9095321328749113306e-37, -2.69536938616169189e-37, 1.9096731758702465152e-37, -1.3885308426986502429e-37, 1.0348071522251170897e-37, -7.895170286820189599e-38, 6.1601036791640358926e-38, -4.9101808205945411315e-38, 3.994641287554558951e-38, -3.3139324868403182382e-38, 2.8011372345112539961e-38, -2.4105257902290214916e-38, 2.1103686667444868275e-38, -1.8783457821005023693e-38, 1.6985661919809396375e-38, -1.5596040985686701268e-38, 1.4531867623354274083e-38, -1.3733073575134458911e-38, 1.3156193132795211797e-38, -1.2770202478775098053e-38, 1.2553659571340599541e-38, -1.2492755513989236073e-38, 1.2580022245918135246e-38, -1.2813530074444447323e-38, 1.3196468864619204783e-38, -1.3737049209287366587e-38, 1.4448691477987981283e-38, -1.5350495877564259921e-38},
    {1, 1, 1, 2.2080283083317470879e-6, -2.8645658927105171426e-10, 2.0828661056308347074e-13, -5.2299632353305440596e-16, 3.2830981481457948702e-18, -4.1510701318390062674e-20, 9.1041811064120705559e-22, -3.1124615845127521356e-23, 1.5328693037472463578e-24, -1.0244920393888717986e-25, 8.8720260734562778411e-27, -9.597948834377940171e-28, 1.2595886766589028663e-28, -1.9578720438883669243e-29, 3.5340072868084426703e-30, -7.2864402165304741374e-31, 1.6923344153105842992e-31, -4.375591783282673663e-32, 1.2466762919873774489e-32, -3.8799069978348316293e-33, 1.3089372943449002827e-33, -4.7548669772601409413e-34, 1.848924808742039275e-34, -7.6558772735244817953e-35, 3.3600813771321854839e-35, -1.5566220555286937867e-35, 7.5836826978718830557e-36, -3.8724580464037880644e-36, 2.0662799537391519081e-36, -1.1489379201019512917e-36, 6.6408853609904585227e-37, -3.9809781297699670577e-37, 2.4699309348902252735e-37, -1.5830074254659870545e-37, 1.0462202764733223946e-37, -7.11875798855754071e-38, 4.9794221402932609151e-38, -3.5755960441132518639e-38, 2.6324348982556010633e-38, -1.9846823550739445576e-38, 1.5306255264431836482e-38, -1.2062699710815663993e-38, 9.7051109266615921699e-39, -7.9642866764809828477e-39, 6.6606788104892184813e-39, -5.6724797376446808561e-39, 4.915739173130756452e-39, -4.3317556647487453338e-39, 3.8789474466050159524e-39, -3.5275293304083189738e-39, 3.2559767261429846516e-39, -3.0486462320583443403e-39, 2.8941566018876330935e-39, -2.7842776421007672607e-39, 2.7131641542460728281e-39, -2.6768287227428640156e-39, 2.6727836132505440804e-39, -2.6998059280763485313e-39, 2.7577961455917057716e-39, -2.8477111715018297324e-39, 2.9715609113133253847e-39, -3.1324633530432551583e-39},
    {1, 1, 1, 1.9015341576771288381e-6, -2.1924917955255159672e-10, 1.4470268704416100003e-13, -3.3484679563984333971e-16, 1.9590214344829755893e-18, -2.3283424791360898145e-20, 4.8326667572442524349e-22, -1.5720548954880036707e-23, 7.3998950829502838853e-25, -4.744681735342801506e-26, 3.9543465978508838151e-27, -4.1282398535720502495e-28, 5.2405273770420484611e-29, -7.8956594119268946537e-30, 1.3839569950483201815e-30, -2.7754105245231706582e-31, 6.2789414384223485841e-32, -1.5834135263837346523e-32, 4.4053812755500230724e-33, -1.3402680916529367656e-33, 4.4244317601459032988e-34, -1.5741188982300353214e-34, 5.9998422840401112784e-35, -2.4370823450352134441e-35, 1.0499975902262956394e-35, -4.7782797161317989732e-36, 2.2881543584247076094e-36, -1.1490999220942912123e-36, 6.0333599423045487809e-37, -3.3028145016395897008e-37, 1.8803403763597713638e-37, -1.1107502307831938013e-37, 6.7937536812687199601e-38, -4.2941620554819759342e-38, 2.799960683378092818e-38, -1.8802717585484637967e-38, 1.2984621755477747258e-38, -9.2081381444739110152e-39, 6.6970928507823640765e-39, -4.9894305297483091418e-39, 3.803473091276116874e-39, -2.9636165769895057001e-39, 2.3580536941968113467e-39, -1.9141703397365358627e-39, 1.5839193667649902422e-39, -1.3349487636814697733e-39, 1.1451168552213665606e-39, -9.9903990530839378147e-40, 8.8588182693873269395e-40, -7.9791637390166110979e-40, 7.295778324579350259e-40, -6.7682562981465866736e-40, 6.3671384698618244753e-40, -6.0709656954184899116e-40, 5.8642473059239813088e-40, -5.7360563685831219304e-40, 5.6790628589628801308e-40, -5.6888801749477144653e-40, 5.7636429573860376069e-40, -5.9037629494695363137e-40, 6.1118296928495487244e-40, -6.3926374845215784953e-40},
    {1, 1, 1, 1.6376862983967245862e-6, -1.6782073954353660753e-10, 1.0053580882066885797e-13, -2.1439903069642528358e-16, 1.1690250067645217397e-18, -1.3060597114262642403e-20, 2.565442308323510527e-22, -7.9407387060900263603e-24, 3.5725270407696367783e-25, -2.1975316611815925064e-26, 1.7626100277717591849e-27, -1.775746694829706736e-28, 2.1804737289317021011e-29, -3.1843598535864806899e-30, 5.420101191678467346e-31, -1.0572281345883845174e-31, 2.3297876818561722196e-32, -5.7303556595112941694e-33, 1.5568364274996385028e-33, -4.6301144918840930396e-34, 1.4956358140934406302e-34, -5.211543841244576322e-35, 1.9471083779047719743e-35, -7.7584533909669774681e-36, 3.2813799853681631887e-36, -1.4668634214350569827e-36, 6.9043087927697881959e-37, -3.4100331118140631707e-37, 1.7618098713080252019e-37, -9.4951435641704259575e-38, 5.3244738242259926787e-38, -3.0993653354619744839e-38, 1.8688073147290302306e-38, -1.1649402523608557618e-38, 7.4939446205123185798e-39, -4.9666866984664358317e-39, 3.386175086580913882e-39, -2.3715099767209510021e-39, 1.7039025559860658576e-39, -1.2544134645442273966e-39, 9.4519518640889768112e-40, -7.2816411515857861436e-40, 5.7297627057143706919e-40, -4.6009132228815554672e-40, 3.7668420419470098691e-40, -3.1418540425679547595e-40, 2.6677218640959900662e-40, -2.3042600275192212266e-40, 2.0233334297015326992e-40, -1.804986215343686286e-40, 1.6349020465605731747e-40, -1.5027139960360566013e-40, 1.4008651595134418816e-40, -1.3238317421586565509e-40, 1.267588401514775498e-40, -1.2292380968478987523e-40, 1.2067556100338524721e-40, -1.1988112278906263061e-40, 1.2046524145837164403e-40, -1.2240288686069225286e-40, 1.257151532904846088e-40, -1.3046797868302156574e-40},
    {1, 1, 1, 1.4105347498611433178e-6, -1.284639093336097516e-10, 6.9854314683057431676e-14, -1.3728654702448186469e-16, 6.9764909846620904015e-19, -7.3266916124514675669e-21, 1.3619665007827291937e-22, -4.0112794513357311705e-24, 1.7248619331525783876e-25, -1.0178694382939572483e-26, 7.857178921347018654e-28, -7.6388156549345424224e-29, 9.0730992357342800152e-30, -1.2843542457642871231e-30, 2.1228590456210619341e-31, -4.0275335597094196615e-32, 8.6452039093965984669e-33, -2.0739476860232600563e-33, 5.5021378665143859023e-34, -1.5996345989886816087e-34, 5.0561881208798904736e-35, -1.7255370020746461473e-35, 6.3193067899215506489e-36, -2.4700692860622411769e-36, 1.0255427757577914168e-36, -4.5033616462494980373e-37, 2.0834550497876677595e-37, -1.0120184735689010329e-37, 5.1450297335311677112e-38, -2.7299075380316178798e-38, 1.5078077566266242608e-38, -8.6488464224719396749e-39, 5.14100837116668047e-39, -3.1605155608146258053e-39, 2.0058483116225814573e-39, -1.3120245794350160009e-39, 8.831176866386772273e-40, -6.1081142707617200151e-40, 4.3354307904574240119e-40, -3.1539841356538858468e-40, 2.3490471102637430602e-40, -1.7892276604130458025e-40, 1.3923507651008066167e-40, -1.1059527640062113366e-40, 8.9588203267507199342e-41, -7.3949713006331005558e-41, 6.2152752900526657977e-41, -5.3150728032047578244e-41, 4.6215558305070334928e-41, -4.0833770665018108198e-41, 3.6638773710870287584e-41, -3.3366061973448121698e-41, 3.0823181241548645458e-41, -2.8869342220344650717e-41, 2.7401436093744232569e-41, -2.6344362976368957046e-41, 2.5644312878870141639e-41, -2.5264104430983070017e-41, 2.5179988799374426995e-41, -2.5379525736531956774e-41, 2.5860273110152807651e-41, -2.6629124263230806018e-41},
    {1, 1, 1, 1.2149622204448843791e-6, -9.8343093448947753372e-11, 4.8539281166224756797e-14, -8.7914598693127513519e-17, 4.1636886763679179627e-19, -4.1103693646596372294e-21, 7.231005911243254459e-23, -2.0264369419490901615e-24, 8.3283956200107373827e-26, -4.7149515317494710614e-27, 3.5027189806368734747e-28, -3.2862402395272789616e-29, 3.7756240560923420211e-30, -5.1805492434020207675e-31, 8.3150186777725099223e-32, -1.5343977065509285023e-32, 3.2082076018149195582e-33, -7.5065840397596003084e-34, 1.944680649920626704e-34, -5.5268565879954610478e-35, 1.7094206477231523085e-35, -5.7136089965720289861e-36, 2.0510541124800154482e-36, -7.8645066470966998161e-37, 3.2053791445626256775e-37, -1.3826502827848838631e-37, 6.2874774869525539358e-38, -3.0036316188991258296e-38, 1.5026060215265346245e-38, -7.84915180917043978e-39, 4.2701550187438159692e-39, -2.4136372148198304498e-39, 1.4143617184900609929e-39, -8.5751280232749815662e-40, 5.3692558097652545916e-40, -3.466135845099568247e-40, 2.3033299762766314909e-40, -1.5733225113279590955e-40, 1.1031845444754799727e-40, -7.930612158513045378e-41, 5.8383531767811499194e-41, -4.3967355687997609765e-41, 3.3836785089064418627e-41, -2.6586277259801108288e-41, 2.1308489296675512343e-41, -1.7406658040591957926e-41, 1.448133529731513605e-41, -1.226070290341130289e-41, 1.0556923426895993035e-41, -9.2383306667055796009e-42, 8.2114252599785237738e-42, -7.4090409731391304104e-42, 6.7824563939831250571e-42, -6.2960680400021803769e-42, 5.9237513217688684263e-42, -5.6463503028176753964e-42, 5.4499339039393414535e-42, -5.3245810531430917636e-42, 5.263537713318541382e-42, -5.2626413494991542008e-42, 5.3199434480234461665e-42, -5.4354845531463587257e-42},
    {1, 1, 1, 1.0465671655789447489e-6, -7.5289293626401343858e-11, 3.3730321002680346016e-14, -5.6301665829228900961e-17, 2.4851168409723697027e-19, -2.3061163732877617637e-21, 3.839356993035150666e-23, -1.0237898051020845282e-24, 4.0215726720072051209e-26, -2.1841878237277462526e-27, 1.5616063326784819488e-28, -1.4138398553119248993e-29, 1.5712651365971692441e-30, -2.0897505405090426371e-31, 3.2571140837035060669e-32, -5.846075157370996615e-33, 1.1906315199248432999e-33, -2.7171561814519466567e-34, 6.8737355855565913328e-35, -1.9096919299286488875e-35, 5.7796612290807096699e-36, -1.8920145780399206675e-36, 6.6575206645729626618e-37, -2.5041570520676820399e-37, 1.001919360676084816e-37, -4.2453705918961880828e-38, 1.8975642632635733541e-38, -8.9152318001798331902e-39, 4.3886415313546844373e-39, -2.2569673183889349289e-39, 1.2093974898485639848e-39, -6.7361780358883297649e-40, 3.8913512419620928051e-40, -2.3267571236217302778e-40, 1.4373345535668079544e-40, -9.1575005177509244353e-41, 6.0078834959206023191e-41, -4.0528089113633945163e-41, 2.8073195419129943393e-41, -1.9942595889968546525e-41, 1.4511649198239509624e-41, -1.0804952534745429681e-41, 8.2235114762892202751e-42, -6.3915509708267846503e-42, 5.0685329230677454217e-42, -4.0975294200680540559e-42, 3.3743072757132186894e-42, -2.8284551754912992203e-42, 2.4116502282624675538e-42, -2.0902359012296119357e-42, 1.8404496570256956706e-42, -1.6453065442422775589e-42, 1.4925343865859033789e-42, -1.3731871242136065182e-42, 1.2807020344170605971e-42, -1.2102517604373376634e-42, 1.1582950270977347417e-42, -1.1222633130729848184e-42, 1.1003421135663484635e-42, -1.0913192994916694636e-42, 1.0944822311689611468e-42, -1.1095514341496568469e-42},
    {1, 1, 1, 9.0156333506679594963e-7, -5.7643271968476709814e-11, 2.344088713017018542e-14, -3.605853863502960159e-17, 1.4833448305863850188e-19, -1.2939229411423542259e-21, 2.0386616569944613659e-23, -5.1726778087974424952e-25, 1.9420367702325601831e-26, -1.0118816260875663221e-27, 6.9624938985196091416e-29, -6.0831447213607579914e-30, 6.5393899740970584223e-31, -8.4302437465982398297e-32, 1.275938577219358065e-32, -2.227501033300972204e-33, 4.4189527945311296954e-34, -9.8358961716478136103e-35, 2.4297659405438239476e-35, -6.5989604450466010218e-36, 1.9542626045066103375e-36, -6.265641546551070221e-37, 2.1611008690978589481e-37, -7.9740463580681109857e-38, 3.1319386686754222767e-38, -1.3036049151661828747e-38, 5.7272173810363034874e-39, -2.646340599283427911e-39, 1.2818648041635791436e-39, -6.4901531573446549626e-40, 3.4254814526973492236e-40, -1.8801056688239404781e-40, 1.0706992867794554404e-40, -6.3137682153660530244e-41, 3.8479443377319833821e-41, -2.4195548432329577151e-41, 1.5671623628544304968e-41, -1.0440508418043699308e-41, 7.1443493061485959987e-42, -5.0151487548518628588e-42, 3.6072009260182494102e-42, -2.6554769472718309834e-42, 1.9987232735276580934e-42, -1.5366754102952226367e-42, 1.2056993927310723257e-42, -9.6461924833640278688e-43, 7.862991549684085229e-43, -6.5254489231823383392e-43, 5.509579231729300864e-43, -4.7295984173281761596e-43, 4.1253092043485050544e-43, -3.6539181938719833151e-43, 3.2846480974529089218e-43, -2.9951405040612107054e-43, 2.7690229086376411606e-43, -2.5942438559899651229e-43, 2.4619221650659139699e-43, -2.3655452516817666039e-43, 2.3004081765845243393e-43, -2.2632215485415497529e-43, 2.2518403012255659527e-43, -2.2650812529023317507e-43},
    {1, 1, 1, 7.7669347382424637828e-7, -4.4135648071394106315e-11, 1.6291214955407323273e-14, -2.3095163892634059289e-17, 8.8544908795886801218e-20, -7.2604216686587547305e-22, 1.0825752777503947171e-23, -2.6136437274376391631e-25, 9.3787581876330879327e-27, -4.6880887623510507277e-28, 3.1044488000083187644e-29, -2.6174747588609201589e-30, 2.7217701019649905252e-31, -3.4010441434274580305e-32, 4.9986534460936000625e-33, -8.4878540974840405447e-34, 1.640166090203743759e-34, -3.5607358426562869863e-35, 8.5893945396173323732e-36, -2.2804168277654554161e-36, 6.608303263989441106e-37, -2.075071771950555504e-37, 7.015587846673205766e-38, -2.5393494369806471133e-38, 9.7908467065497231444e-39, -4.0031595235259536506e-39, 1.7286909991887603468e-39, -7.8557087092798677329e-40, 3.7443882387496358852e-40, -1.866427641589613716e-40, 9.70288127111726808e-41, -5.2478027453604354013e-41, 2.9461926213592119791e-41, -1.7133763280293289438e-41, 1.0302111446548715341e-41, -6.3932330895349829235e-42, 4.0882084807928756453e-42, -2.689761177819699344e-42, 1.8182770857192689898e-42, -1.2612828919674952079e-42, 8.9670677942254192453e-43, -6.526626589630741223e-43, 4.8581911175560161101e-43, -3.6947459995020009566e-43, 2.8682854790535994742e-43, -2.2709957936943527094e-43, 1.8323884825834658518e-43, -1.5055599507211970827e-43, 1.2587779047056786958e-43, -1.0702365978637102659e-43, 9.2473146824686566924e-44, -8.1151655345942059579e-44, 7.2290282119468882596e-44, -6.5332796478053518876e-44, 5.9873075743850040948e-44, -5.5612502999557284721e-44, 5.233063650418729389e-44, -4.9864828242357534653e-44, 4.8095965082491961871e-44, -4.6938464404906045294e-44, 4.6333278805964513404e-44, -4.6243075223033374414e-44},
    {1, 1, 1, 6.6915517287630340646e-7, -3.3795223933736913873e-11, 1.1322914008999674811e-14, -1.4793109190548181169e-17, 5.2857992443330114793e-20, -4.0741868760556213354e-22, 5.7490591528800432847e-24, -1.3206967446234234406e-25, 4.5295913181762755274e-27, -2.1721397178183109082e-28, 1.3842993359994435853e-29, -1.1263223635816164939e-30, 1.1328998163346360265e-31, -1.3721774663011613435e-32, 1.9584034973285108273e-33, -3.2344749619018927637e-34, 6.0881059230450227836e-35, -1.2891144851794829476e-35, 3.0365929648198894591e-36, -7.8809554536020250072e-37, 2.2347190816484289937e-37, -6.8726880696226732036e-38, 2.2776082850172822032e-38, -8.0870871083020664516e-39, 3.0609285935035288116e-39, -1.2293788943849511819e-39, 5.2181556076753647517e-40, -2.3321203450998433221e-40, 1.0938190397852917692e-40, -5.3677629643838080383e-41, 2.7485635073264601523e-41, -1.4648687359495505372e-41, 8.1073837086676866452e-42, -4.6498919210616865544e-42, 2.7583519335877618443e-42, -1.6893964542796888994e-42, 1.0665422148906666673e-42, -6.9299766136495274687e-43, 4.6278942613256783968e-43, -3.1722482653948462876e-43, 2.229238586079457352e-43, -1.6042091040547815964e-43, 1.1809255074931859344e-43, -8.8840913610557594583e-44, 6.8238848569653092643e-44, -5.3469081126071677321e-44, 4.2704464825291054466e-44, -3.4738550997070500529e-44, 2.876111887345698387e-44, -2.4219284991608305305e-44, 2.0730069453848679843e-44, -1.8024446473421979479e-44, 1.5910981921162776381e-44, -1.4251851465795281831e-44, 1.2946804245061650351e-44, -1.1922300960153500098e-44, 1.1124069917155494907e-44, -1.0511952485568429235e-44, 1.0056303683062763559e-44, -9.7354645874179307924e-45, 9.5339852965639135753e-45, -9.4413832332804381044e-45},
    {1, 1, 1, 5.7653715243951175338e-7, -2.5878884873250043549e-11, 7.8702351855158693248e-15, -9.4759513188858654851e-18, 3.1556060781925039444e-20, -2.2863629141605630471e-22, 3.0532377263357035562e-24, -6.6739820706790263357e-26, 2.1877512381160336548e-27, -1.0064795352911379192e-28, 6.173064286116210433e-30, -4.8469463030041147301e-31, 4.7158163439479334429e-32, -5.5364785111094480196e-33, 7.6732024161564781735e-34, -1.2326364941664010898e-34, 2.2599660682652783912e-35, -4.6673307166863962437e-36, 1.0735838553989562011e-36, -2.723760230706445114e-37, 7.5575543789319074718e-38, -2.2763839731881464101e-38, 7.3946799277966680094e-39, -2.575651843035505146e-39, 9.5699910784490004187e-40, -3.7756697427949348964e-40, 1.5752231085414529726e-40, -6.9237589603711698004e-41, 3.1954754388839685117e-41, -1.5438350712303863611e-41, 7.7863915380530949392e-42, -4.0892656406795173292e-42, 2.2311343671030251611e-42, -1.2619971209224288238e-42, 7.3858166125990569616e-43, -4.4644506943214419912e-43, 2.7825852554368899013e-43, -1.785563225867538418e-43, 1.1779645117715297282e-43, -7.978977691174131462e-44, 5.5422749878285561628e-44, -3.9432890115475520798e-44, 2.8707530390502881173e-44, -2.1363227625211283245e-44, 1.6235528504243272915e-44, -1.2589675525720261389e-44, 9.9530099986208624733e-45, -8.0158719785157169227e-45, 6.5718534658579379988e-45, -5.4811068242175260028e-45, 4.6474134175451764012e-45, -4.0036114322259241841e-45, 3.5021881677724813527e-45, -3.1091148076802965654e-45, 2.7997485418735012039e-45, -2.5560724115905423254e-45, 2.3648131777263946684e-45, -2.2161435375971996341e-45, 2.1027786145404334561e-45, -2.0193422006077110229e-45, 1.9619202865678709524e-45, -1.9277467624917149879e-45},
    {1, 1, 1, 4.9676450311961937719e-7, -1.9817993757989677093e-11, 5.4706818251312905466e-15, -6.0703061668123068763e-18, 1.8839936002077430203e-20, -1.2831398028019838712e-22, 1.6216200291153676986e-24, -3.372808141198301585e-26, 1.056723970507392609e-27, -4.6638746079251735904e-29, 2.7529373757747672709e-30, -2.0859239386931511656e-31, 1.9631208332495612933e-32, -2.2339926353154998654e-33, 3.0066018917370814324e-34, -4.6977624335144865729e-35, 8.3897002784236051463e-36, -1.6899368540296633803e-36, 3.7958599764526839916e-37, -9.4142061229890972817e-38, 2.5560212234578008879e-38, -7.5403106994747317498e-39, 2.4009580237135841836e-39, -8.2036483979671770377e-40, 2.9922282194046205493e-40, -1.1596504377901456214e-40, 4.7554538975408401891e-41, -2.0556907341006822039e-41, 9.3357743706219241338e-42, -4.4405149448556919963e-42, 2.2059291688117498417e-42, -1.1416074688829860073e-42, 6.1403848353266669362e-43, -3.4253009745793174779e-43, 1.9777538288647496985e-43, -1.1798569505384271335e-43, 7.2601197208843543505e-44, -4.6009082374707303852e-44, 2.9985124890763857028e-44, -2.0070223576535058686e-44, 1.3779849272542999842e-44, -9.6935112849578410736e-45, 6.9790135922688423034e-45, -5.1374262661947413764e-45, 3.8630119700786635579e-45, -2.9644985180855603661e-45, 2.3198530844151113071e-45, -1.849757221240249044e-45, 1.5017403580582378965e-45, -1.240509523136215455e-45, 1.0419496813908356037e-45, -8.8933788647286877285e-46, 7.7091564687506616579e-46, -6.7830828465762607533e-46, 6.0548078310908247334e-46, -5.4803856053625870474e-46, 5.0275321474164001331e-46, -4.6723705251772696154e-46, 4.3971736521708747336e-46, -4.1887848644442785766e-46, 4.0375055544889539014e-46, -3.936309459754540295e-46},
    {1, 1, 1, 4.2805162054268416481e-7, -1.5177394396202160581e-11, 3.802935357138067941e-15, -3.8888597251606812473e-18, 1.1248641326430194417e-20, -7.2015645836331964829e-23, 8.6131438387636195167e-25, -1.704599650526336575e-26, 5.1044548048113814e-28, -2.1612897185800734719e-29, 1.2277673717735061311e-30, -8.9774494055296585733e-32, 8.1726223588439186695e-33, -9.0147589256385669164e-34, 1.1781469655233207893e-34, -1.7904877708825129437e-35, 3.1146930939547244535e-36, -6.1192285217760576424e-37, 1.3421736057683228958e-37, -3.2540396113840179542e-38, 8.6451381777460454818e-39, -2.4977977867105332881e-39, 7.7960402156956390072e-40, -2.6130711767430740358e-40, 9.3562586589337286691e-41, -3.5619229411444649149e-41, 1.4357082853153990674e-41, -6.1037670243992816505e-42, 2.7276556651148498945e-42, -1.277291747912740771e-42, 6.2498735023064961056e-43, -3.1872242449441474314e-43, 1.6900120369386471525e-43, -9.2974415479587755958e-44, 5.2962721593619347984e-44, -3.1182794125886697495e-44, 1.8943640111197007585e-44, -1.1855947319117981985e-44, 7.6331512737395404312e-45, -5.0487227635899986031e-45, 3.4262986778886333242e-45, -2.3830216473928909953e-45, 1.6967451314675097675e-45, -1.2355169089324549386e-45, 9.1920004290712729046e-46, -6.9809141579345418672e-46, 5.4074297275778303759e-46, -4.2687729100090237549e-46, 3.4318333903813869197e-46, -2.8077362037197035777e-46, 2.3361813292799858821e-46, -1.9756319085491446891e-46, 1.6970655195940727165e-46, -1.4799323027359831969e-46, 1.309501702615753186e-46, -1.1750962715328176639e-46, 1.068900410027094665e-46, -9.8514698801999393501e-47, 9.1955572282808698603e-47, -8.6894154045174179845e-47, 8.3093925811941432996e-47, -8.038089843303078831e-47},
    {1, 1, 1, 3.6886174813885380708e-7, -1.1624055536557150118e-11, 2.643745710582846408e-15, -2.4914798119059020204e-18, 6.7165186397865922104e-21, -4.0420652165532934644e-23, 4.5750718579833856443e-25, -8.6154269930493321883e-27, 2.46581702860070238e-28, -1.0016195931971598991e-29, 5.4759511184975171651e-31, -3.8639473328555656545e-32, 3.4025113078383296806e-33, -3.6378957890272644095e-34, 4.6168607666492990154e-35, -6.8245723361573759921e-36, 1.1563993922722964655e-36, -2.2158821259059326675e-37, 4.7460356618766288159e-38, -1.1248270853385125295e-38, 2.9241740917970508623e-39, -8.2746401839336030212e-40, 2.5315551002389583361e-40, -8.3237543359377872227e-41, 2.9257252623429489087e-41, -1.0941220306835108892e-41, 4.3347520445803145784e-42, -1.8124329286187130956e-42, 7.9698944527157863726e-43, -3.6742669496201449974e-43, 1.770821241609185031e-43, -8.8988169403906167876e-44, 4.6516586815055730713e-44, -2.5237830095032328359e-44, 1.4183786741414117888e-44, -8.2418468412712571851e-45, 4.9431857261477249302e-45, -3.0552924942294670384e-45, 1.943236775430758107e-45, -1.2700905409960257854e-45, 8.5198081779450862555e-46, -5.858665497348261771e-46, 4.1253711595306352989e-46, -2.9714993173375960063e-46, 2.187347926613013764e-46, -1.6439825727473441528e-46, 1.2605066438117887713e-46, -9.851791810096849087e-47, 7.8429851529215306935e-47, -6.3553043927094668475e-47, 5.2382980698564162182e-47, -4.3890358298620248673e-47, 3.7360632166043669097e-47, -3.2290921274244256178e-47, 2.832276372522271467e-47, -2.5197626764931965335e-47, 2.2727072154028351674e-47, -2.0772494488884816703e-47, 1.9231198338214492527e-47, -1.802672878654727158e-47, 1.7102093831630134785e-47, -1.6414979346753679105e-47},
    {1, 1, 1, 3.1787219660860199137e-7, -8.9030864748179494499e-12, 1.8379902063969144398e-15, -1.5963032591062262388e-18, 4.0106189998673732037e-21, -2.2688347362110175592e-23, 2.4302853765381338854e-25, -4.3546613395350036462e-27, 1.191229734841314865e-28, -4.6421145129496302716e-30, 2.4424533009031485682e-31, -1.6631548069182140247e-32, 1.4166448269968512461e-33, -1.4681473933993529943e-34, 1.8093282418598008463e-35, -2.6013740737004715473e-36, 4.2936213356538231791e-37, -8.0245360660239877515e-38, 1.6783270822699045631e-38, -3.888409564181676259e-39, 9.8914013114083669156e-40, -2.7413486593283149043e-40, 8.2209880569528887688e-41, -2.6516157582564633299e-41, 9.1493067687396674885e-42, -3.3610141010587491733e-42, 1.3088372441016424399e-42, -5.3820691433085990304e-43, 2.3288361342853368562e-43, -1.0569991649206437975e-43, 5.0176637955266567402e-44, -2.4847074791776938475e-44, 1.2804105000868183327e-44, -6.8511568945784592409e-45, 3.7987213571618677817e-45, -2.1784991777095998009e-45, 1.2899526573796617038e-45, -7.8739504305758295958e-46, 4.9473305088077279163e-46, -3.1952967185391959299e-46, 2.1186429367129089399e-46, -1.440432040668021437e-46, 1.0030736019059453062e-46, -7.1470354108407288415e-47, 5.2053395016939543304e-47, -3.8717337180534970547e-47, 2.9384801019151495827e-47, -2.2737921387193214038e-47, 1.7925032519647497336e-47, -1.4385990637449050987e-47, 1.1746195605847046678e-47, -9.7511443312284837289e-48, 8.2253272830945087327e-48, -7.0459959194248396819e-48, 6.1261631326576815895e-48, -5.4034258464082601791e-48, 4.8325138240694306853e-48, -4.3802574681547165528e-48, 4.0221476633557024518e-48, -3.7399574302420935664e-48, 3.5200810336345629184e-48, -3.3523642067497621523e-48},
    {1, 1, 1, 2.739444404952330419e-7, -6.8193897889735533158e-12, 1.2778769066882128048e-15, -1.022812193658730209e-18, 2.3949759102798829523e-21, -1.2735764786987909759e-23, 1.2910387912466333763e-25, -2.2011754532450048874e-27, 5.7551006706625547772e-29, -2.1515509212549681852e-30, 1.0894711463506059464e-31, -7.1590751948581498551e-33, 5.8985486540524506784e-34, -5.9253206275826009975e-35, 7.0910532773817443497e-36, -9.9163759372931386286e-37, 1.5942720841299621701e-37, -2.9061366658269639969e-38, 5.9353322469620699131e-39, -1.3442531767545801884e-39, 3.3460717975430700211e-40, -9.0824344403819934005e-41, 2.6698292909341148759e-41, -8.4474333722396736038e-42, 2.8613150990022729319e-42, -1.0325183020081895325e-42, 3.9521178286997062487e-43, -1.5983042672800763234e-43, 6.8053137329981548718e-44, -3.0408951337106614823e-44, 1.4218417335858943098e-44, -6.938109133977853587e-45, 3.5246295685457819917e-45, -1.859938921231231045e-45, 1.0174323778791513253e-45, -5.7585498722742067216e-46, 3.3663827112341120184e-46, -2.0293428878526439489e-46, 1.2596184221906848552e-46, -8.039158321202333836e-47, 5.2687622263571635152e-47, -3.5416832899111980152e-47, 2.439076679128835346e-47, -1.7190919769601704291e-47, 1.2388055236008516857e-47, -9.118777820461978496e-48, 6.850515487257991563e-48, -5.2481855307297699271e-48, 4.0969569104660333674e-48, -3.2566124270562313539e-48, 2.6340689355687840422e-48, -2.1665311117562552776e-48, 1.810985804574451898e-48, -1.5375428461651178254e-48, 1.3251479653495630315e-48, -1.15878173785395173e-48, 1.0276035814346249804e-48, -9.237054700576438164e-49, 8.4126450472714380536e-49, -7.759599107930501215e-49, 7.2456763233437685235e-49, -6.846757430911711966e-49},
    {1, 1, 1, 2.3609840540752537827e-7, -5.2236252025067388387e-12, 8.8849854390322502943e-16, -6.5538788886838914568e-19, 1.4302534325819289712e-21, -7.1493961511170312152e-24, 6.8587270011330176944e-26, -1.1126976616689384248e-27, 2.7805620281608473579e-29, -9.9726296593570778502e-31, 4.8599018190940085776e-32, -3.0817931448304097519e-33, 2.4561318410032240566e-34, -2.3915330430539314991e-35, 2.779242714359078679e-36, -3.7802933876006794546e-37, 5.9200244153609303781e-38, -1.0525300285513623667e-38, 2.0991130703655484903e-39, -4.6474262924830340584e-40, 1.1319703807647164211e-40, -3.0092803323077049042e-41, 8.6709234340168294314e-42, -2.691295242388254537e-42, 8.9488158596769386949e-43, -3.1721051352991493792e-43, 1.1934288069982357217e-43, -4.746702798097675163e-44, 1.9887480618419033408e-44, -8.7488433904732896891e-45, 4.0292420501050265428e-45, -1.9374450307373200513e-45, 9.7028683373519456854e-46, -5.04958699054461998e-46, 2.72518555029307102e-46, -1.5222685452927183006e-46, 8.7856847219948892185e-47, -5.2304686430705168735e-47, 3.2072254951001419809e-47, -2.0227042353291969163e-47, 1.3103335021313119612e-47, -8.7086149295449673058e-48, 5.9311720675292961948e-48, -4.1351827831195878131e-48, 2.9483539867960099459e-48, -2.147782044869162931e-48, 1.5971517115537867801e-48, -1.2114068349378561275e-48, 9.3645144112230665812e-49, -7.3725003343464964427e-49, 5.9071699777575763295e-49, -4.8138959529068699979e-49, 3.9874874523506258092e-49, -3.3553242203498443225e-49, 2.8665704299978815873e-49, -2.4851723611753585244e-49, 2.1852470383178621728e-49, -1.9480038893351319276e-49, 1.7596631842214216902e-49, -1.6100313718230737431e-49, 1.4915150043673465003e-49, -1.3984312424429839685e-49},
    {1, 1, 1, 2.0349035647110216465e-7, -4.0014707608699506331e-12, 6.1779708638936809008e-16, -4.1997409929342534568e-19, 8.5417432876486045972e-22, -4.0136123000819219106e-24, 3.6439256513834651612e-26, -5.6249865088921292354e-28, 1.3434888597153893788e-29, -4.6226348808405734568e-31, 2.1680091558998240126e-32, -1.3266975037124442009e-33, 1.0227749010568648225e-34, -9.6530110330863374297e-36, 1.0893416230952382265e-36, -1.4411856602186338521e-37, 2.1983987180546232475e-38, -3.8121929148632731953e-39, 7.4241809664001210667e-40, -1.6068149834447551857e-40, 3.82962983590113493e-41, -9.9711431795613767963e-42, 2.8162366324104500888e-42, -8.5747180649044756325e-43, 2.7989000938371901923e-43, -9.7458411867453952917e-44, 3.6040025871520234064e-44, -1.409764496331791594e-44, 5.8121034620491946236e-45, -2.517223587210204093e-45, 1.1418719668230517914e-45, -5.4105274778599172155e-46, 2.6712144334219708992e-46, -1.3709922250321514841e-46, 7.2997596388033739634e-47, -4.0243095147409089077e-47, 2.293029821655842666e-47, -1.3481795492633197474e-47, 8.1666125497693256239e-48, -5.0895121018635333282e-48, 3.2589451920926798108e-48, -2.1414622368751989926e-48, 1.4423729241571233407e-48, -9.9474617453306025553e-49, 7.0174299572026523136e-49, -5.0590124882152055174e-49, 3.7238404196236976357e-49, -2.7963583966212299507e-49, 2.1405782161248050216e-49, -1.6691120782285345877e-49, 1.3248105363261427769e-49, -1.0696715464678299022e-49, 8.7802238540484916412e-50, -7.3225733883811783351e-50, 6.2013013465384356608e-50, -5.3300755531991598126e-50, 4.6472651227852630917e-50, -4.1083563116950698499e-50, 3.6808533024278477488e-50, -3.3408068838338787696e-50, 3.0704236580924680221e-50, -2.856401645010351278e-50},
    {1, 1, 1, 1.7539388183733840273e-7, -3.0654061227161827508e-12, 4.2959182441192810754e-16, -2.6913350026714005536e-19, 5.1015397121685072415e-22, -2.2533194821084877129e-24, 1.9360511769091248096e-26, -2.8437221452037372276e-28, 6.491677548268084299e-30, -2.1428456955122294013e-31, 9.671996926136296503e-33, -5.7116524518336588205e-34, 4.2592182342873926576e-35, -3.8964640446398491139e-36, 4.2699538565914521715e-37, -5.4945967666498776588e-38, 8.1641485985069043028e-39, -1.3808189053013209454e-39, 2.6259276162579006697e-40, -5.5557247868269213092e-41, 1.2956867740602121116e-41, -3.3040663103530906484e-42, 9.1473318527302262079e-43, -2.7321205497980605992e-43, 8.7544875035430580738e-44, -2.994419174358321557e-44, 1.088416662600461192e-44, -4.1871894413811104588e-45, 1.698667659093271119e-45, -7.2429328538143237402e-46, 3.236182332523286144e-46, -1.5110239367528184093e-46, 7.3542584545639129953e-47, -3.7225080731094189186e-47, 1.9554313206829052325e-47, -1.0639298721291208317e-47, 5.9850160336103830032e-48, -3.4751727251241183873e-48, 2.0795814600457569777e-48, -1.2806824229122462717e-48, 8.1057610948009178037e-49, -5.2661507236532003135e-49, 3.5078105681180463138e-49, -2.3930478254825097159e-49, 1.6703138523890749655e-49, -1.1916887265973760581e-49, 8.6827538271532078252e-50, -6.4553111684872797673e-50, 4.8932620942160735331e-50, -3.7790073537615917252e-50, 2.9713213121471279738e-50, -2.3769810943319149092e-50, 1.9334519685876271535e-50, -1.5981388709491426743e-50, 1.3416047790504525719e-50, -1.1432251160366422679e-50, 9.88361683972562459e-51, -8.6649872438565968875e-51, 7.6999696631337749464e-51, -6.9325011296903396067e-51, 6.3210689112379015246e-51, -5.8347059731661170067e-51},
    {1, 1, 1, 1.5118353651379287699e-7, -2.34842516987598679e-12, 2.9873545419531477731e-16, -1.7247801958581432434e-19, 3.0470308718879717916e-22, -1.2651179110045277768e-24, 1.028691325121516727e-26, -1.4377182388327515416e-28, 3.136900759068190285e-30, -9.9337492565363782324e-32, 4.3151131663884613449e-33, -2.4590799494589677188e-34, 1.7737840641999128387e-35, -1.5728944388229570774e-36, 1.6737990097821168125e-37, -2.0949456864905244059e-38, 3.0320502402376471431e-39, -5.0017230172228744528e-40, 9.2883380937425765001e-41, -1.9210410394866026026e-41, 4.3839372056941239381e-42, -1.0948979229392013212e-42, 2.9712610150761071776e-43, -8.7056448769110674874e-44, 2.7383889977077740837e-44, -9.2008287220848036403e-45, 3.2872017381206302992e-45, -1.2437117626357881193e-45, 4.9648324011450419803e-46, -2.0841463870557203032e-46, 9.1721183545906014793e-47, -4.2201139628005046664e-47, 2.0248372574081167542e-47, -1.0107817808633846978e-47, 5.2383881574242064316e-48, -2.8129093285572515593e-48, 1.5622194401764039031e-48, -8.9583115172559228232e-49, 5.2957933350077149583e-49, -3.2227591799506015359e-49, 2.0161909020839686153e-49, -1.2950819756963199304e-49, 8.5313117738195701293e-50, -5.7572036207034465952e-50, 3.9759342257331320735e-50, -2.8072495647758920846e-50, 2.0246270184331042629e-50, -1.4902620606748789193e-50, 1.1186313072599818743e-50, -8.5564002114838372574e-51, 6.6644847532342117871e-51, -5.2822886343832702253e-51, 4.2577710447558919306e-51, -3.4880798197864413237e-51, 2.9026019035156626483e-51, -2.4521739134248990631e-51, 2.1021101609839582165e-51, -1.8276325031954254815e-51, 1.6108334495906657511e-51, -1.4386317707975408708e-51, 1.3013792191968105802e-51, -1.1918999416028739452e-51}};

  /** Encode a graph, return the first word */
  long encode1c(long c[]) {
    int ib = 0;
    long code = 0;

    for (int i = 0; i < n - 1; i++) {
      code |= (c[i] >> (i + 1)) << ib;
      ib += n - 1 - i;
    }
    return code;
  }

  /** From code to the connectivity matrix */
  void decode1c(long c[], long code)
  {
    int ib = 0;

    // empty the connectivity matrix
    for (int i = 0; i < n; i++) c[i] = 0;

    // make links
    for (int i = 0; i < n - 1; i++) {
      for (int j = i + 1; j < n; j++) {
        if ( ((code >> ib) & 1) != 0 )
          linkc(c, i, j);
        ++ib;
      }
    }
  }

  private static DiagramMap [] dgmap = new DiagramMap [ DiagramMap.NMAX + 1 ];

  /** Get the diagram index */
  short getMapIdx() {
    if (dgmap[n] == null)
      dgmap[n] = new DiagramMap(n);
    long code = encode1c(c);
    return dgmap[n].map[(int) code];
  }

  private static boolean fbnrMapInitialized = false;
  private static double [][][] fbnrMap = new double [DiagramMap.NMAX + 1][][];

  /** Get Fb and Nr from the lookup table */
  double [] getFbNrMap() {
    if (n < 1 || n > DiagramMap.NMAX) {
      System.out.printf("cannot use lookup table with n %d\n", n);
      return null;
    }
    short mapID = getMapIdx();
    short nUnique = dgmap[n].nUnique;
    if ( fbnrMap[n] == null ) {
      long cc [] = new long [n];
      fbnrMap[n] = new double [nUnique][2];
      for (int ig = 0; ig < nUnique; ig++) {
        double fb, nr;
        decode1c(cc, dgmap[n].first[ig]);
        if ( biconnectedc(n, cc) ) {
          fb = getFbc(cc);
          /*
          // checking code
          double fb1 = getFbDirectc(cc);
          double fb2 = getFbRJWc(cc);
          if (abs(fb - fb1) > 0.01 || abs(fb1 - fb2) > 0.01) {
            System.out.printf("corruption %f %f %f, ig %d, code %x\n",
                fb, fb1, fb2, ig, dgmap[n].first[ig]);
            printc(cc);
            System.exit(1);
          }
          */
          nr = getNrc(cc);
        } else {
          fb = nr = 0;
        }
        fbnrMap[n][ig][0] = fb;
        fbnrMap[n][ig][1] = nr;
      }
    }
    return fbnrMap[n][mapID];
  }
}



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



/** Bitwise operations.
 *  All functions are `static' for no instance is needed */
class Bits {
  /** Return a long integer with only ith lowest bit being 1 */
  public static long makeBit(int i) {
    return (long) 1 << i;
  }

  /** Return a long integer with the lowest i bits being 1 */
  public static long makeBitsMask(int i) {
    return ((long) 1 << i) - 1;
  }

  private static final int [] bitrem = {
    -1,  0,  1, 39,  2, 15, 40, 23,  3, 12, 16, 59, 41, 19, 24, 54,
     4, -1, 13, 10, 17, 62, 60, 28, 42, 30, 20, 51, 25, 44, 55, 47,
     5, 32, -1, 38, 14, 22, 11, 58, 18, 53, 63,  9, 61, 27, 29, 50,
    43, 46, 31, 37, 21, 57, 52,  8, 26, 49, 45, 36, 56,  7, 48, 35,
     6, 34, 33};

  /** From bit to id */
  public static int bit2id(long bit) {
    return bitrem[(int) (bit % 67)];
  }

  /** Seek the least significant bit */
  public static int bitFirst(long x) {
    long bit = x & (-x);
    return bit2id(bit);
  }

  private static final int [] byteBitCount = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8 };

  /** Count the number of one-bits in x */
  public static int bitCount(long x) {
    return byteBitCount[(int) (x         & 0xffl)] +
           byteBitCount[(int) ((x >>  8) & 0xffl)] +
           byteBitCount[(int) ((x >> 16) & 0xffl)] +
           byteBitCount[(int) ((x >> 24) & 0xffl)] +
           byteBitCount[(int) ((x >> 32) & 0xffl)] +
           byteBitCount[(int) ((x >> 40) & 0xffl)] +
           byteBitCount[(int) ((x >> 48) & 0xffl)] +
           byteBitCount[(int) ((x >> 56) & 0xffl)];
  }
}



/** Average and standard deviation */
class Ave {
  double cnt, xsm, x2sm;

  public void clear() {
    cnt = xsm = x2sm = 0;
  }

  public void add(double x) {
    cnt += 1;
    xsm += x;
    x2sm += x*x;
  }

  public double getAve() {
    return (cnt > 0) ? xsm/cnt : 0;
  }

  public double getStd() {
    if (cnt <= 1) return 0;
    double av = xsm/cnt;
    return sqrt(x2sm/cnt - av*av);
  }
}




class MCSamp {
  public int D = 3; // the dimension
  public int N = 7; // the order
  public Diagram g, ng;

  int step; // the current simulation step
  double mcAmp = 1.5; // amplitude of randomly moving a particle
  int sampFreq = 1; // sampling frequency
  long mcacc = 0, mctot = 0;
  double x[][], xi[], xo[];
  Ave avfb = new Ave(), avnr = new Ave();
  double xyz[]; // 3D coordinates to output

  double fb, nr; // star and ring contents
  private boolean fbdirty = true;

  Random rng = new Random();

  /** Constructor */
  MCSamp(int dim, int order) {
    init(dim, order);
  }

  /** Initialize system */
  public void init(int dim, int order) {
    D = dim;
    N = order;
    mcAmp = 1.5/D;
    x = new double [N][D];
    xyz = new double [(N + 1) * 3];
    double rr = 0.7/(2.*sin(PI/N));
    for (int i = 0; i < N; i++) {
      x[i][0] = rr * cos(2*i*PI/N);
      x[i][1] = rr * sin(2*i*PI/N);
      for (int j = 2; j < D; j++) x[i][j] = 0;
    }
    xi = new double [D];
    xo = new double [D];
    g  = new Diagram(N, D, x);
    ng = new Diagram(N, D, x);
    clearData();
  }

  /** Clear trajectory data */
  void clearData() {
    step = 0;
    mcacc = mctot = 0;
    avfb.clear();
    avnr.clear();
  }

  public int moveAtom = -1;
  public boolean moveAcc = false;

  /** Step of the Metropolis algorithm */
  public void mcStep() {
    int i, j, k;

    moveAtom = i = (int) (N * rng.nextDouble());
    // displace particle i
    for (k = 0; k < D; k++) {
      xo[k] = x[i][k];
      xi[k] = xo[k] + (2 * rng.nextDouble() - 1) * mcAmp;
    }
    ng.copy(g);
    for (j = 0; j < N; j++) {
      if (j == i) continue;
      // distance between i and j
      double dx, x2 = 0;
      for (k = 0; k < D; k++) {
        dx = xi[k] - x[j][k];
        x2 += dx * dx;
      }
      if (x2 < 1) {
        ng.link(i, j);
      } else {
        ng.unlink(i, j);
      }
    }
    mctot += 1;
    if ( ng.biconnected() ) {
      mcacc += 1;
      for (k = 0; k < D; k++)
        x[i][k] = xi[k];
      g.copy(ng);
      //for (int ii = 0; ii < N; ii++) System.out.print("x" + ii + ": " + x[ii][0] + " " + x[ii][1] + " " + x[ii][2] + "; ");
      //System.out.println("moving " + i);
      moveAcc = true;
      fbdirty = true;
    } else {
      moveAcc = false;
    }
    if ((step + 1) % sampFreq == 0) {
      if ( fbdirty ) {
        if (N <= DiagramMap.NMAX) {
          double [] arr2 = g.getFbNrMap();
          fb = arr2[0];
          nr = arr2[1];
        } else {
          fb = g.getFb();
          nr = g.getNr();
        }
        fbdirty = false;
      }
      avfb.add(fb);
      avnr.add(nr);
    }
    step++;
  }

  /** Remove the center of mass motion */
  void rmcom() {
    for (int k = 0; k < D; k++) {
      double xc = 0;
      for (int i = 0; i < N; i++)
        xc += x[i][k];
      xc /= N;
      for (int i = 0; i < N; i++)
        x[i][k] -= xc;
      xi[k] -= xc;
      xo[k] -= xc;
    }
  }

  /** Get the virial coefficients */
  public double getVir() {
    double fb = avfb.getAve();
    double sc = avnr.getAve();
    if (sc == 0) return 0;
    double Bring = abs( Diagram.BringArr[D][N] );
    return -fb/sc*Bring;
  }
}



public class VirSampApp extends JApplet implements ActionListener
{
  MCSamp mc = new MCSamp(3, 6);
  int delay = 100; // interval of repainting
  int speed = 10000; // mc steps per frame
  Timer timer;

  MyCanvas canvas; // 3D animation
  MyScheme scheme; // schematic plot

  JPanel cpnl, spnl, mpnl;
  JTextField      tDim     = new JTextField("    " + mc.D);
  JTextField      tOrder   = new JTextField("    " + mc.N);
  JTextField      tDelay   = new JTextField("   " + delay);
  JTextField      tSpeed   = new JTextField("   " + speed);
  JTextField      tMCAmp   = new JTextField("   " + mc.mcAmp*mc.D);
  JTextField      tSampFreq = new JTextField("   " + mc.sampFreq);
  JToggleButton   bStart   = new JToggleButton("Start", false);
  JButton         bReset   = new JButton("Reset");
  JButton         bMDS     = new JButton("MDS mode 0");
  JButton         bClear   = new JButton("Clear");
  JLabel          lStatus  = new JLabel("Status");
  JTextField      tVir     = new JTextField("   0");
  JTextField      tAvStar  = new JTextField("   0");
  JTextField      tAvRing  = new JTextField("   0");
  JTextField      tAcc     = new JTextField("   0");
  public static final long serialVersionUID = 1L;

  /** initialize the handler */
  public void init() {
    Container box = getContentPane();
    box.setLayout(new BorderLayout());

    // cpnl is the panel for controls, such as buttons, textfields, etc.
    cpnl = new JPanel();
    cpnl.setLayout( new GridLayout(4, 6) ); // four rows, six columns
    box.add(cpnl, BorderLayout.NORTH);

    // add controls
    cpnl.add(bStart);
    bStart.addActionListener(this);

    cpnl.add(bReset);
    bReset.addActionListener(this);

    cpnl.add(new JLabel(" Dimension (D):"));
    tDim.addActionListener(this);
    cpnl.add(tDim);

    cpnl.add(new JLabel(" Order (n):"));
    tOrder.addActionListener(this);
    cpnl.add(tOrder);

    cpnl.add(new JLabel(" Delay (ms):"));
    tDelay.addActionListener(this);
    cpnl.add(tDelay);

    cpnl.add(new JLabel(" Steps/frame:"));
    tSpeed.addActionListener(this);
    cpnl.add(tSpeed);

    cpnl.add(new JLabel(" Samp. freq.:"));
    tSampFreq.addActionListener(this);
    cpnl.add(tSampFreq);

    cpnl.add(new JLabel(" Bn/B2^(n-1):"));
    tVir.setEditable(false);
    cpnl.add(tVir);

    cpnl.add(new JLabel(" < HSFb >:"));
    tAvStar.setEditable(false);
    cpnl.add(tAvStar);

    cpnl.add(new JLabel(" < Ring >:"));
    tAvRing.setEditable(false);
    cpnl.add(tAvRing);

    cpnl.add(new JLabel(" Move size D*\u0394X:"));
    tMCAmp.addActionListener(this);
    cpnl.add(tMCAmp);

    cpnl.add(new JLabel(" Acc. ratio:"));
    tAcc.setEditable(false);
    cpnl.add(tAcc);

    cpnl.add(bMDS);
    bMDS.addActionListener(this);

    cpnl.add(bClear);
    bClear.addActionListener(this);

    // spnl is the panel at the bottom
    spnl = new JPanel(); // create a panel for status
    box.add(spnl, BorderLayout.SOUTH);
    lStatus.setFont(new Font("Courier", Font.BOLD, 14));
    spnl.add(lStatus);

    // create a panel for animation and schematic plot
    mpnl = new JPanel();
    mpnl.setLayout(new GridLayout(1, 2));
    box.add(mpnl, BorderLayout.CENTER);

    // scheme is the schematic diagram
    scheme = new MyScheme();
    scheme.setMC(mc);
    mpnl.add(scheme);

    // canvas is the place for 3D animation
    canvas = new MyCanvas();
    canvas.setMC(mc);
    mpnl.add(canvas);

    Font font = new Font("Arial", Font.BOLD, 14);
    tVir.setFont(font);

    timer = new Timer(delay, this);
    timer.start(); timer.stop();
  }

  public void start() {
    scheme.start();
    canvas.start();
  }

  public void update(Graphics g) {
    paint(g);
  }

  //int frame = 0;

  /** paint the applet */
  public void paint(Graphics g) {
    //System.out.println("frame: " + (frame++));
    double vir = mc.getVir();
    lStatus.setText("steps = " + mc.step + ", "
         + "D = " + mc.D + ", "
         + "n = " + mc.N + ", "
         + "vir = " + String.format("%+.5e", vir) + ", "
         + "avfb = " + String.format("%+8.5f", mc.avfb.getAve()) + ", "
         + "avnr = " + String.format("%+8.5f", mc.avnr.getAve())
         + "");
    tVir.setText(" " + String.format("%.5e", vir));
    tAvStar.setText(" " + String.format("%g", mc.avfb.getAve()));
    tAvRing.setText(" " + String.format("%g", mc.avnr.getAve()));
    tAcc.setText(" " + String.format("%g", mc.mctot > 0 ? 1.0*mc.mcacc/mc.mctot : 0.0));
    scheme.repaint();
    if (canvas.model != null) {
      canvas.update();
      canvas.repaint();
    }
    cpnl.repaint();
    spnl.repaint();
  }

  /** Event handler */
  public void actionPerformed(ActionEvent e) {
    Object src = e.getSource();

    // the regular time for the regular animation
    if (src == timer) {
      for (int i = 0; i < speed; i++) // sample a few steps
        mc.mcStep();
      //mc.g.print();
      repaint(); // calls the paint() function above
      return;
    }

    //System.out.println( "trigger: " + src.toString() );

    // the user may change the value of some parameters
    // but forget to press Enter to trigger an event
    // so we always refresh all parameters when an effect occurs
    int dim = mc.D;
    try {
      dim = Integer.parseInt(tDim.getText().trim());
    } catch (NumberFormatException err) {}
    if (dim < 2) {
      dim = 2;
      tDim.setText(" " + dim);
    } else if (dim > 100) {
      dim = 100;
      tDim.setText(" " + dim);
    }

    int n = mc.N;
    try {
      n = Integer.parseInt(tOrder.getText().trim());
    } catch (NumberFormatException err) {}
    if (n < 2) {
      n = 2;
      tOrder.setText(" " + n);
    } else if (n > Diagram.NMAX) {
      n = Diagram.NMAX;
      tOrder.setText(" " + n);
    }

    try {
      delay = Integer.parseInt(tDelay.getText().trim());
    } catch (NumberFormatException err) {}
    if (delay < 1) {
      delay = 1;
      tDelay.setText(" " + delay);
    }

    try {
      speed = Integer.parseInt(tSpeed.getText().trim());
    } catch (NumberFormatException err) {}
    if (speed < 1) {
      speed = 1;
      tSpeed.setText(" " + speed);
    }

    int sampFreq = mc.sampFreq;
    try {
      sampFreq = Integer.parseInt(tSampFreq.getText().trim());
    } catch (NumberFormatException err) {}
    if (sampFreq < 1) {
      sampFreq = 1;
      tSpeed.setText(" " + sampFreq);
    }
    mc.sampFreq = sampFreq;

    double amp = mc.mcAmp * mc.D;
    try {
      amp = Double.parseDouble(tMCAmp.getText().trim());
    } catch (NumberFormatException err) {}
    if (amp < 0.) {
      amp = -amp;
      tMCAmp.setText(" " + amp);
    }
    mc.mcAmp = amp / dim;

    if (src == tDelay) {
      if ( timer.isRunning() ) {
        timer.stop();
        timer = new Timer(delay, this);
        timer.start();
      }
    }

    if (src == tSampFreq || src == tMCAmp
        || src == bReset || src == bClear) {
      mc.clearData();
    }

    if (src == bStart) {
      boolean on = bStart.isSelected();
      if (on) {
        timer.restart();
        bStart.setText("Pause");
      } else {
        timer.stop();
        bStart.setText("Resume");
      }
    }

    if (src == tDim || src == tOrder || src == bReset
     || dim != mc.D || n != mc.N) {
      if ( src == bReset ) {
        if ( timer.isRunning() ) timer.stop();
        bStart.setSelected(false);
        bStart.setText("Start");
      }
      mc.init(dim, n);
      canvas.amat.unit(); // reset the view
      canvas.bStarted = false;
    }

    // change the mode of multidimensional scaling
    if (src == bMDS) {
      int mdsmode = 0;
      mdsmode = Integer.parseInt( bMDS.getText().substring(9) );
      mdsmode = (mdsmode + 1) % 4;
      //System.out.println("" + mdsmode);
      bMDS.setText( String.format("MDS mode %d", mdsmode) );
      scheme.useMDS = ((mdsmode & 0x1) != 0);
      canvas.useMDS = ((mdsmode & 0x2) != 0);
    }

    // in case the schematic diagram has not been initialized
    // initialize it
    if ( !scheme.bStarted ) {
      scheme.newImgBuf();
    }

    // in case the canvas has not been initialized, initialize it
    if ( !canvas.bStarted ) {
      canvas.gauge();
    }
    repaint();
  }
}



/** Panel for drawing the schematic diagram */
class MyScheme extends JPanel
    implements MouseListener, MouseMotionListener, MouseWheelListener {
  Image img;
  Graphics imgG;
  Dimension imgSize;
  MCSamp mc;
  public boolean bStarted = false;
  public static final long serialVersionUID = 2L;

  public MyScheme() {
    super();
    setBackground( Color.WHITE );
  }

  void setMC(MCSamp m) { mc = m; }

  public void start() {
    newImgBuf();
    addMouseListener(this);
    addMouseMotionListener(this);
    addMouseWheelListener(this);
  }

  public void newImgBuf() {
    if (getSize().width == 0 || getSize().height == 0) {
      System.out.println("scheme not ready.\n");
      return;
    }
    //System.out.println("scheme " + getSize().toString());
    img = createImage(getSize().width, getSize().height);
    if (imgG != null) { imgG.dispose(); }
    imgG = img.getGraphics();
    imgSize = getSize();
    bStarted = true;
  }

  public void update(Graphics g) {
    if (img == null)
      g.clearRect(0, 0, getSize().width, getSize().height);
    paintComponent(g);
  }

  private static final Color colorLine         = Color.BLACK;
  private static final Color colorLargeCircle  = new Color(230, 230, 230);
  private static final Atom atomNormal = new Atom(  0,   0, 255);
  private static final Atom atomFailed = new Atom(255, 100, 100);
  private static final Atom atomMoved  = new Atom(100, 255, 100);

  boolean useMDS = false;

  protected void paintComponent(Graphics g) {
    super.paintComponent(g);
    if ( img == null ) return;
    if ( !imgSize.equals(getSize()) )
      newImgBuf();
    imgG.setColor( Color.WHITE );
    int w = getSize().width;
    int h = getSize().height;
    imgG.fillRect(0, 0, w, h);

    if ( mc == null) return;
    int np = mc.N;

    // make lines thicker
    Graphics2D imgG2 = (Graphics2D) imgG;
    int lineWidth = (int) (40./np);
    if (lineWidth > 4) lineWidth = 4;
    else if (lineWidth < 1) lineWidth = 1;
    imgG2.setStroke(new BasicStroke(lineWidth));

    int rad = 400 / np; // radius of the balls
    if (rad > 40) rad = 40;
    int [][] xy;

    if (useMDS) {
      xy = getMDS(np, rad);
    } else {
      int radius = (int) (w/2 - rad - 5); // the big circle
      xy = new int [np][2];
      // compute the 2D coordinates around the circle
      for (int i = 0; i < np; i++) {
        double theta = 2*PI*i/np - PI*.5;
        xy[i][0] = (int) (w/2 + radius * cos(theta) + .5);
        xy[i][1] = (int) (h/2 + radius * sin(theta) + .5);
      }
      // draw the big circle
      imgG.setColor( colorLargeCircle );
      imgG.drawOval(w/2 - radius, h/2 - radius, 2*radius, 2*radius);
    }

    // draw connections
    imgG.setColor( colorLine );
    for (int i = 0; i < np; i++) {
      for (int j = i + 1; j < np; j++) {
        if ( mc.g.isLinked(i, j) ) {
          // draw a line between i and j
          imgG2.drawLine(xy[i][0], xy[i][1], xy[j][0], xy[j][1]);
        }
      }
    }

    // draw the small circles
    for (int i = 0; i < np; i++) {
      Atom atom = atomNormal;
      if (i == mc.moveAtom)
        atom = mc.moveAcc ? atomMoved : atomFailed;
      atom.paint(imgG, xy[i][0], xy[i][1], 15, rad);
    }

    g.drawImage(img, 0, 0, this);
  }

  /** Get the best coordinates from multidimensional scaling */
  int [][] getMDS(int np, int rad) {
    double [][] xy0 = new double [np][2];
    int [][] xy = new int [np][2];
    for (int i = 0; i < np; i++) {
      double theta = 2*PI*i/np - PI*.5;
      xy0[i][0] = cos(theta);
      xy0[i][1] = sin(theta);
    }
    MDS mds = new MDS(mc.x);
    mds.min0(xy0, 0);
    int w = getSize().width, h = getSize().height;
    int margin = 5 + rad;

    double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
    for (int i = 0; i < np; i++) {
      if (xy0[i][0] > xmax) xmax = xy0[i][0];
      else if (xy0[i][0] < xmin) xmin = xy0[i][0];
      if (xy0[i][1] > ymax) ymax = xy0[i][1];
      else if (xy0[i][1] < ymin) ymin = xy0[i][1];
    }
    xmax = max(xmax, -xmin);
    ymax = max(ymax, -ymin);

    double scl = min((.5 * w - margin) / xmax, (.5 * h - margin) / ymax);
    for (int i = 0; i < np; i++) {
      xy[i][0] = (int) (w / 2 + scl * xy0[i][0]);
      xy[i][1] = (int) (h / 2 + scl * xy0[i][1]);
    }
    xy0 = null;
    return xy;
  }

  /** Event handling */
  public void mouseClicked(MouseEvent e) { }
  public void mousePressed(MouseEvent e) {
    //prevx = e.getX();
    //prevy = e.getY();
    // consume this event so that it will not be processed
    // in the default manner
    e.consume();
  }
  public void mouseReleased(MouseEvent e) { }
  public void mouseEntered(MouseEvent e) { }
  public void mouseExited(MouseEvent e) { }

  public void mouseDragged(MouseEvent e) {
    //int x = e.getX();
    //int y = e.getY();
    //double xtheta = 360.0 * (prevy - y) / getSize().height;
    //double ytheta = 360.0 * (x - prevx) / getSize().width;
    repaint();
    //prevx = x;
    //prevy = y;
    e.consume();
  }

  public void mouseMoved(MouseEvent e) { }

  public void mouseWheelMoved(MouseWheelEvent e) {
    //int notches = e.getWheelRotation();
    //repaint();
    e.consume();
  }
}



/** Panel for drawing particles */
class MyCanvas extends JPanel
    implements MouseListener, MouseMotionListener, MouseWheelListener {
  Image img;
  Graphics imgG;
  Dimension imgSize;

  double real2scrn;
  double xtheta, ytheta;
  double zoomscale = 1.0f;
  public Matrix3D amat = new Matrix3D(); // view matrix
  private Matrix3D tmat = new Matrix3D(); // temporary matrix
  int prevx, prevy;
  String message;
  XYZModel model;
  MCSamp mc;
  public boolean bStarted = false;

  public static final long serialVersionUID = 2L;

  public MyCanvas() {
    super();
  }

  void setMC(MCSamp m) { mc = m; }

  boolean useMDS = false; // use multidimensional scaling to transform coordinates

  /** Get 3D position */
  double [][] getPos3D() {
    int n = mc.N, D = mc.D;
    double r[][] = new double[n + 1][3];

    mc.rmcom();

    // copy the direct coordinates
    for (int i = 0; i < n; i++) {
      r[i][0] = mc.x[i][0];
      r[i][1] = mc.x[i][1];
      r[i][2] = (D == 2) ? 0.0 : mc.x[i][2];
    }
    // the last is the trial position
    double [] xlast0 = mc.moveAcc ? mc.xo : mc.xi;
    r[n][0] = xlast0[0];
    r[n][1] = xlast0[1];
    r[n][2] = (D == 2) ? 0 : xlast0[2];

    if ( D > 3 && useMDS ) {
      double r0[][] = new double[n + 1][D];
      double [] xlast = mc.moveAcc ? mc.xo : mc.xi;

      for (int i = 0; i < n; i++)
        for (int k = 0; k < D; k++)
          r0[i][k] = mc.x[i][k];
      for (int k = 0; k < D; k++)
        r0[n][k] = xlast[k];

      MDS mds = new MDS(r0);
      mds.min0(r, 0);
    }
    return r;
  }

  /** Copy coordinates from MCSamp to XYZModel */
  public void update() {
    model.updateXYZ( getPos3D() );
    model.moveAtom = mc.moveAtom;
    model.moveAcc = mc.moveAcc;
  }

  /** Determine scaling factors based on the current coordinates */
  void gauge() {
    int dim = useMDS ? 3 : mc.D;
    double w = 1.5*Math.pow(mc.N, 1./dim);
    if (model == null)
      model = new XYZModel( getPos3D() );
    double f1 = getSize().width  / w;
    double f2 = getSize().height / w;
    real2scrn = 0.7f * (f1 < f2 ? f1 : f2);
    System.out.println("real2scrn " + real2scrn + " w " + w);
    newImgBuf();
  }

  public void start() {
    gauge();
    addMouseListener(this);
    addMouseMotionListener(this);
    addMouseWheelListener(this);
  }

  public void newImgBuf() {
    if (getSize().width == 0 || getSize().height == 0) {
      System.out.println("canvas not ready.\n");
      return;
    }
    //System.out.println("canvas " + getSize().toString());
    img = createImage(getSize().width, getSize().height);
    if (imgG != null) { imgG.dispose(); }
    imgG = img.getGraphics();
    imgSize = getSize();
    bStarted = true;
  }

  public void update(Graphics g) {
    if (img == null)
      g.clearRect(0, 0, getSize().width, getSize().height);
    paintComponent(g);
  }

  protected void paintComponent(Graphics g) {
    super.paintComponent(g);
    if (model != null) {
      model.mat.unit();
      model.mat.mult(amat);
      // real2scrn is the default scaling from real coordinate
      // to the screen coordinats, it must be multiplied by
      // the zoom scale to get the actual scaling
      double real2scrn1 = real2scrn * zoomscale;
      model.mat.scale(real2scrn1, real2scrn1, real2scrn1);
      model.mat.translate(getSize().width / 2, getSize().height / 2, 8);
      model.real2scrn = real2scrn1; // tell the XYZModel the scaling factor
      model.transformed = false;
      if ( img != null ) {
        if ( !imgSize.equals(getSize()) )
          newImgBuf();
        imgG.setColor(Color.BLACK);
        imgG.fillRect(0, 0, getSize().width, getSize().height);
        model.paint(imgG);
        g.drawImage(img, 0, 0, this);
      } else
        model.paint(g);
    } else if ( message != null ) {
      g.drawString("Error in model:", 3, 20);
      g.drawString(message, 10, 40);
    }
  }

  /** Event handling */
  public void mouseClicked(MouseEvent e) { }
  public void mousePressed(MouseEvent e) {
    prevx = e.getX();
    prevy = e.getY();
    // consume this event so that it will not be processed
    // in the default manner
    e.consume();
  }
  public void mouseReleased(MouseEvent e) { }
  public void mouseEntered(MouseEvent e) { }
  public void mouseExited(MouseEvent e) { }

  public void mouseDragged(MouseEvent e) {
    int x = e.getX();
    int y = e.getY();
    tmat.unit();
    double xtheta = 360.0 * (prevy - y) / getSize().height;
    double ytheta = 360.0 * (x - prevx) / getSize().width;
    tmat.xrot(xtheta);
    tmat.yrot(ytheta);
    amat.mult(tmat);
    repaint();
    prevx = x;
    prevy = y;
    e.consume();
  }

  public void mouseMoved(MouseEvent e) { }

  public void mouseWheelMoved(MouseWheelEvent e) {
    int notches = e.getWheelRotation();
    if ((zoomscale -= 0.05f*notches) < 0.09999f)
      zoomscale = 0.1f;
    repaint();
  }
}



/** A set of atoms with 3D coordinates */
class XYZModel {
  double realxyz[][]; // 3D real coordinates 3*np
  int scrnxyz[][];  // 3D screen coordinates in pixels
                    // only the first two dimensions are used in drawing
  int zorder[]; // z-order, larger is closer to the viewer
  int np = -1;

  boolean transformed;
  Matrix3D mat; // rotation/scaling/translation matrix for the conversion
                // from realxyz[] to scrnxyz[]
  double real2scrn = 100.0; // real size --> screen size, to be set by MyCanvas
  double ballSize = 0.5; // ball size (radius) in terms of the real coordinates
                          // it is 0.5 because we are simulating hard spheres

  int moveAtom = -1;
  boolean moveAcc = false;

  // various atom types
  static Atom atomNormal = new Atom(  0,   0, 255);
  // failed move
  static Atom atomTrial  = new Atom(255,   0,   0);  // failed trial postion
  static Atom atomFailed = new Atom( 50,   0, 100);  // new and old position
  // successful move
  static Atom atomMoved  = new Atom(  0, 255,   0);  // trial and new position
  static Atom atomOldpos = new Atom(  0, 100, 100);  // old position

  /** Constructor from the real coordinates */
  XYZModel(double r[][]) {
    mat = new Matrix3D();
    updateXYZ(r);
  }

  /** Refresh coordinates */
  void updateXYZ(double r[][]) {
    int n = r.length;
    if (n != np) {
      np = n;
      realxyz = new double[np][3];
      scrnxyz = new int[np][3];
      zorder = new int[np];
    }
    for (int i = 0; i < np; i++) {
      realxyz[i][0] = r[i][0];
      realxyz[i][1] = r[i][1];
      realxyz[i][2] = r[i][2];
    }
    transformed = false;
    //for (int i = 0; i < np; i++)
    //  System.out.println("i " + i + ": " + r[i][0] + " " + r[i][1] + " " + r[i][2]);
  }

  /** Paint this model to the graphics context `g' */
  void paint(Graphics g) {
    if (realxyz == null || np <= 0) return;

    // transform the coordinates
    if ( !transformed ) {
      mat.transform(realxyz, scrnxyz, np);
      transformed = true;
    }

    // bubble sort z-order
    // zorder[0] is the fartherest from the viewer
    // zorder[np - 1] is the nearest to the viewer
    for (int i = 0; i < np; i++)
      zorder[i] = i;
    for (int i = 0; i < np; i++) {
      // find the particle with the smallest z
      int jm = i, k;
      for (int j = i + 1; j < np; j++)
        if (scrnxyz[zorder[j]][2] < scrnxyz[zorder[jm]][2])
          jm = j;
      if (jm != i) {
        k = zorder[i];
        zorder[i] = zorder[jm];
        zorder[jm] = k;
      }
    }

    double zmin = scrnxyz[zorder[0]][2];
    double zmax = scrnxyz[zorder[np - 1]][2];
    //System.out.println("" + zmin + " " + zmax);

    // radius in pixels
    double radius = real2scrn * ballSize;

    //System.out.println("real2scrn " + real2scrn + ", ball " + ballSize);
    for (int iz = 0; iz < np; iz++) {
      int i = zorder[iz], greyscale;
      if (zmin == zmax) { // two-dimensional case
        greyscale = 15;
      } else {
        greyscale = (int) (15.999 * (scrnxyz[i][2] - zmin) / (zmax - zmin));
      }
      // scrnxyz[3*i] and scrnxyz[3*i + 1] are the (x, y) coordinates on screen
      Atom atom = atomNormal;
      if ( i == moveAtom )
        atom = moveAcc ? atomMoved : atomFailed;
      if ( i == np - 1 ) {
        if ( moveAtom < 0 ) {
          atom = null;
        } else if ( moveAcc ) { // the original position before pos
          atom = atomOldpos;
        } else { // the failed trial position
          atom = atomTrial;
        }
      }
      if ( atom != null )
        atom.paint(g, scrnxyz[i][0], scrnxyz[i][1], greyscale, radius);
    }
  }
}



class Atom {
  private final static int R = 120;
  private final static int hx = 45;  // (hx, hy) is the offset of the spot light from the center
  private final static int hy = 45;
  private final static int bgGrey = 0; // further atoms are to be blended with this color
  private final static int nBalls = 16;
  private final static double spotlightmag = .9f; // spotlight brightness, 1.f for full
  private static int maxr; // maximal intensity

  private static byte data[];
  static { // data[] is a bitmap image of the ball of radius R
    data = new byte[R * 2 * R * 2];
    for (int Y = -R; Y < R; Y++) {
      int x0 = (int) (Math.sqrt(R * R - Y * Y) + 0.5);
      for (int X = -x0; X < x0; X++) {
        // sqrt(x^2 + y^2) gives distance from the spot light
        int x = X + hx, y = Y + hy;
        int r = (int) (Math.sqrt(x * x + y * y) + 0.5);
        // set the maximal intensity to the maximal distance
        // (in pixels) from the spot light
        if (r > maxr) maxr = r;
        data[(Y + R) * (R * 2) + (X + R)] = (r <= 0) ? 1 : (byte) r;
      }
    }
  }

  // the following variables are atom dependent
  private int Rl = 30, Gl = 10, Bl = 255;
  private Image balls[]; // 0..nBalls-1, at different z distances
  public double ballSize; // ball size

  /** Constructor */
  Atom(int r, int g, int b) {
    setRGB(r, g, b);
  }

  /** Set color */
  void setRGB(int r, int g, int b) {
    Rl = r; Gl = g; Bl = b;
    makeBalls();
  }

  /** Linearly interpolate colors */
  private int blend(int fg, int bg, double fgfactor) {
    return (int) (bg + (fg - bg) * fgfactor);
  }

  // need a component instance to call createImage()
  private static Component component = new Applet();

  /** Prepare ball images with different sizes */
  private void makeBalls() {
    balls = new Image[nBalls];
    byte red[] = new byte[256];
    red[0] = (byte) bgGrey;
    byte green[] = new byte[256];
    green[0] = (byte) bgGrey;
    byte blue[] = new byte[256];
    blue[0] = (byte) bgGrey;
    for (int r = 0; r < nBalls; r++) {
      // smaller b means greyer
      double b = (r + 1. + 7) / (nBalls + 7);
      for (int i = maxr; i >= 1; --i) {
        // contrast of the spotlight
        double d = 1 - (1 - 1. * i / maxr) * spotlightmag;
        red[i]    = (byte) blend(blend(Rl, 255, d), bgGrey, b);
        green[i]  = (byte) blend(blend(Gl, 255, d), bgGrey, b);
        blue[i]   = (byte) blend(blend(Bl, 255, d), bgGrey, b);
      }
      // 256 color model
      IndexColorModel model = new IndexColorModel(
          8, maxr + 1, red, green, blue, 0);
      balls[r] = component.createImage(
          new MemoryImageSource(R * 2, R * 2, model, data, 0, R * 2) );
    }
  }

  /** Draw a ball at screen coordinate (x, y) with a ball index `r'
   *  the `radius' is measured in pixels
   *  (0, 0) represents the top-left corner
   *  x, y, radius are given in terms of pixels
   *  the ball index (gray code) can be 0 to 15 */
  void paint(Graphics gc, int x, int y, int r, double radius) {
    if (balls == null) makeBalls();
    Image img = balls[r]; // r = [0..15]

    int size = (int) (radius * 2 + .5);
    gc.drawImage(img, x - size/2, y - size/2, size, size, null);
    //System.out.println("" + x + " " + y + " " + r + " " + radius);
  }
}



class Matrix3D {
  double xx, xy, xz, xo;
  double yx, yy, yz, yo;
  double zx, zy, zz, zo;
  static final double pi = 3.14159265;

  /** Create a new unit matrix */
  Matrix3D() {
    xx = 1.0;
    yy = 1.0;
    zz = 1.0;
  }

  /** Scale along each axis independently */
  void scale(double xf, double yf, double zf) {
    xx *= xf; xy *= xf; xz *= xf; xo *= xf;
    yx *= yf; yy *= yf; yz *= yf; yo *= yf;
    zx *= zf; zy *= zf; zz *= zf; zo *= zf;
  }

  /** Translate the origin */
  void translate(double x, double y, double z) {
    xo += x;
    yo += y;
    zo += z;
  }

  /** Rotate theta degrees around the y axis */
  void yrot(double theta) {
    theta *= (pi / 180);
    double ct = Math.cos(theta);
    double st = Math.sin(theta);

    double Nxx = (xx * ct + zx * st);
    double Nxy = (xy * ct + zy * st);
    double Nxz = (xz * ct + zz * st);
    double Nxo = (xo * ct + zo * st);

    double Nzx = (zx * ct - xx * st);
    double Nzy = (zy * ct - xy * st);
    double Nzz = (zz * ct - xz * st);
    double Nzo = (zo * ct - xo * st);

    xo = Nxo; xx = Nxx; xy = Nxy; xz = Nxz;
    zo = Nzo; zx = Nzx; zy = Nzy; zz = Nzz;
  }

  /** Rotate theta degrees about the x axis */
  void xrot(double theta) {
    theta *= (pi / 180);
    double ct = Math.cos(theta);
    double st = Math.sin(theta);

    double Nyx = (yx * ct + zx * st);
    double Nyy = (yy * ct + zy * st);
    double Nyz = (yz * ct + zz * st);
    double Nyo = (yo * ct + zo * st);

    double Nzx = (zx * ct - yx * st);
    double Nzy = (zy * ct - yy * st);
    double Nzz = (zz * ct - yz * st);
    double Nzo = (zo * ct - yo * st);

    yo = Nyo; yx = Nyx; yy = Nyy; yz = Nyz;
    zo = Nzo; zx = Nzx; zy = Nzy; zz = Nzz;
  }

  /** Rotate theta degrees about the z axis */
  void zrot(double theta) {
    theta *= pi / 180;
    double ct = Math.cos(theta);
    double st = Math.sin(theta);

    double Nyx = (yx * ct + xx * st);
    double Nyy = (yy * ct + xy * st);
    double Nyz = (yz * ct + xz * st);
    double Nyo = (yo * ct + xo * st);

    double Nxx = (xx * ct - yx * st);
    double Nxy = (xy * ct - yy * st);
    double Nxz = (xz * ct - yz * st);
    double Nxo = (xo * ct - yo * st);

    yo = Nyo; yx = Nyx; yy = Nyy; yz = Nyz;
    xo = Nxo; xx = Nxx; xy = Nxy; xz = Nxz;
  }

  /** Multiply this matrix by a second: M = M*R */
  void mult(Matrix3D rhs) {
    double lxx = xx * rhs.xx + yx * rhs.xy + zx * rhs.xz;
    double lxy = xy * rhs.xx + yy * rhs.xy + zy * rhs.xz;
    double lxz = xz * rhs.xx + yz * rhs.xy + zz * rhs.xz;
    double lxo = xo * rhs.xx + yo * rhs.xy + zo * rhs.xz + rhs.xo;

    double lyx = xx * rhs.yx + yx * rhs.yy + zx * rhs.yz;
    double lyy = xy * rhs.yx + yy * rhs.yy + zy * rhs.yz;
    double lyz = xz * rhs.yx + yz * rhs.yy + zz * rhs.yz;
    double lyo = xo * rhs.yx + yo * rhs.yy + zo * rhs.yz + rhs.yo;

    double lzx = xx * rhs.zx + yx * rhs.zy + zx * rhs.zz;
    double lzy = xy * rhs.zx + yy * rhs.zy + zy * rhs.zz;
    double lzz = xz * rhs.zx + yz * rhs.zy + zz * rhs.zz;
    double lzo = xo * rhs.zx + yo * rhs.zy + zo * rhs.zz + rhs.zo;

    xx = lxx; xy = lxy; xz = lxz; xo = lxo;
    yx = lyx; yy = lyy; yz = lyz; yo = lyo;
    zx = lzx; zy = lzy; zz = lzz; zo = lzo;
  }

  /** Recover the unit matrix */
  void unit() {
    xo = 0; xx = 1; xy = 0; xz = 0;
    yo = 0; yx = 0; yy = 1; yz = 0;
    zo = 0; zx = 0; zy = 0; zz = 1;
  }

  /** Transform np points from v into tv.
   *  v contains the input coordinates in floating point.
   *  Three successive entries in the array constitute a point.
   *  tv ends up holding the transformed points as integers;
   *  three successive entries per point */
  void transform(double v[][], int tv[][], int np) {
    // np can be different from v.length
    for (int i = 0; i < np; i++) {
      double x = v[i][0], y = v[i][1], z = v[i][2];
      tv[i][0] = (int) (xx * x + xy * y + xz * z + xo);
      tv[i][1] = (int) (yx * x + yy * y + yz * z + yo);
      tv[i][2] = (int) (zx * x + zy * y + zz * z + zo);
    }
  }
}



/** Multidimensional scaling */
class MDS {
  int N = -1;
  double [][] mat;

  /** Constructor from the coordinates */
  MDS(double x[][]) {
    N = x.length;
    mat = new double[N][N];
    getDisMat(mat, x);
  }

  /** Compute the distance matrix */
  void getDisMat(double [][] m, double x[][]) {
    int D = x[0].length;
    for (int i = 0; i < N - 1; i++) {
      for (int j = i + 1; j < N; j++) {
        double d2 = 0;
        for (int k = 0; k < D; k++) {
          double dx = x[i][k] - x[j][k];
          d2 += dx * dx;
        }
        m[i][j] = m[j][i] = sqrt(d2);
      }
    }
  }

  int ITER_MAX = 1000;

  /** Get coordinates */
  double min0(double y[][], double tol) {
    if (y.length != N) {
      System.out.printf("dimension mismatch %d vs %d\n", y.length, N);
      System.exit(1);
    }
    int D = y[0].length;
    int npr = N * (N - 1) / 2;

    if (tol <= 0) tol = 1e-4;

    double [][] f  = new double [N][D];
    double [][] yp = new double [N][D];
    double [][] fp = new double [N][D];
    double ene = getForce(y, f);
    double dt = 0.1;
    int iter;

    for (iter = 0; iter < ITER_MAX; iter++) {
      double enep = ene;
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < D; j++) {
          // backup the force
          yp[i][j] = y[i][j];
          fp[i][j] = f[i][j];

          // move along the force
          y[i][j] += f[i][j] * dt;
        }
      }
      ene = getForce(y, f);
      if (ene > enep) {
        dt *= 0.7;
        // recover the force
        for (int i = 0; i < N; i++) {
          for (int j = 0; j < D; j++) {
            y[i][j] = yp[i][j];
            f[i][j] = fp[i][j];
          }
        }
      } else {
        if (abs(ene - enep) < tol * npr * dt)
          break;
        dt *= 1.03f; // attempt to increase the step size
      }
    }
    if (iter >= ITER_MAX) {
      System.out.println("failed to reach convergence in " + iter + " steps");
    }
    f = null;
    fp = null;
    yp = null;
    rmcom(y);
    //print(y);
    return ene;
  }

  final double DMIN = 1e-6;

  /** Compute the force and energy */
  double getForce(double y[][], double f[][]) {
    int D = y[0].length;

    // clear the force
    for (int i = 0; i < N; i++)
      for (int k = 0; k < D; k++)
        f[i][k] = 0;

    double ene = 0;
    double [] dy = new double [D];
    // compute the pair difference
    for (int i = 0; i < N; i++) {
      for (int j = i + 1; j < N; j++) {
        double dy2 = 0;
        for (int k = 0; k < D; k++) {
          dy[k] = y[i][k] - y[j][k];
          dy2 += dy[k] * dy[k];
        }
        double dij = sqrt(dy2);
        double dref = mat[i][j];
        double dsc = dij/dref - 1;
        ene += dsc * dsc;
        if (dij < DMIN) dij = DMIN;
        for (int k = 0; k < D; k++) {
          dy[k] *= -(dij - dref) / (dref * dref * dij);
          f[i][k] += dy[k];
          f[j][k] -= dy[k];
        }
      }
    }
    dy = null;
    return ene;
  }

  /** Remove the center of mass */
  void rmcom(double x[][]) {
    int D = x[0].length;
    int N = x.length;

    for (int k = 0; k < D; k++) {
      double xc = 0;
      for (int i = 0; i < N; i++)
        xc += x[i][k];
      xc /= N;
      for (int i = 0; i < N; i++)
        x[i][k] -= xc;
    }
  }

  /** Print coordinates */
  void print(double x[][]) {
    int D = x[0].length;
    int N = x.length;

    for (int i = 0; i < N; i++) {
      System.out.printf("%2d: ", i);
      for (int j = 0; j < D; j++)
        System.out.printf("%+8.3f ", x[i][j]);
      System.out.println("");
    }
  }

  /** Free memory */
  void finish() {
    mat = null;
  }
}

