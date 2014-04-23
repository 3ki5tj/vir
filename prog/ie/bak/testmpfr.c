/* compute e using MPFR
 *   gcc testmpfr.c -lmpfr -lgmp && ./a.out
 * -lgmp can be omitted
 * */
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

int main(void)
{
  int i, ndigits = 40;
  mpfr_t x, s;
  mpfr_set_default_prec(256);
  mpfr_init_set_d(x, 1, MPFR_RNDD);
  mpfr_init_set_d(s, 2, MPFR_RNDD);
  for (i = 2; i <= 40; i++) {
    mpfr_div_d(x, x, i, MPFR_RNDD);
    mpfr_add(s, s, x, MPFR_RNDD);
    mpfr_printf("%4d: %.*Rf %.*Rf\n", i, ndigits, x, ndigits, s);
  }
  mpfr_clear(s);
  mpfr_clear(x);
  return 0;
}
