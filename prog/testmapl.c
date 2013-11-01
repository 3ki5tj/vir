#include "dgmapl.h"
#include "testutil.h"


static void testspeed(int k, int n, int nsamp, int nedmax)
{
  dg_t *g;
  int t, ned, eql = 1, nequil = 1000, isamp = 0, good = 0, tot = 0;
  int err = 0, st[DGMAPL_NMAX + 1];
  clock_t t0;
  double tsum = 0, rnp = 0.1;

  printf("getchain: speed test for n %d, k %d, nsamp %d, nedmax %d\n",
      n, k, nsamp, nedmax);
  g = dg_open(n);
  dg_full(g);
  for (t = 0; ; t++) {
    dg_rndswitchedge(g, &ned, nedmax, rnp);
    if (eql && t >= nequil) {
      t = 0;
      eql = 0;
      continue;
    }
    if (t % 10 != 0) continue;
    adjustrnp(ned, nedmax, t, 1000000, &good, &tot, &rnp);
    if ( ned > nedmax ) continue;
    t0 = clock();
    err += (dgmapl_getchain(g, k, st) == 0);
    tsum += clock() - t0;
    if (++isamp >= nsamp) break;
  }
  tsum /= CLOCKS_PER_SEC;
  printf("getchain() err %d(%g%%), time %g/%d = %gmcs\n",
    err, 100.*err/nsamp, tsum, nsamp, tsum/nsamp*1e6);
  dg_close(g);
}



static void testbyte3(void)
{
  int x0, x1, i;
  byte3_t b3;
  int arr[] = {17, -129, -654321, -123456, 0x808080, 0x800000, -8388608,
    0x7fffff, 0x8080, 0x7f7f7f, 0x808080, 0x7f7f7f7f, 0x80808080, 0};

  for (i = 0; ; i++) {
    x0 = arr[i];
    if (x0 == 0) break;
    x1 = b3toi( itob3(x0, b3) );
    printf("x0 %12d(0x%8x) --> byte3 0x%02x%02x%02x --> x1 %12d(0x%8x)\n",
        x0, x0, b3[2], b3[1], b3[0], x1, x1);
  }
}



int main(int argc, char **argv)
{
  int k = 7, n = 9, nsamp = 1000000, nedmax = 1000000;

  testbyte3();
  printf("short %d, %x \n", (int16_t) 0x8080, 32639);
  if (argc >= 2) k = atoi(argv[1]);
  if (argc >= 3) n = atoi(argv[2]);
  if (argc >= 4) nsamp = atoi(argv[3]);
  if (argc >= 5) nedmax = atoi(argv[4]);
  testspeed(k, n, nsamp, nedmax);
  return 0;
}

