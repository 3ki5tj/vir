/* compute e using GMP
 *   gcc testgmp.c -lgmp && ./a.out
 * */
#include <stdio.h>
#include "gmp.h"

int main(void)
{
  int i, ndigits = 40;
  mpf_t x, s;
  mpf_set_default_prec(256);
  mpf_init_set_si(x, 1);
  mpf_init_set_si(s, 2);
  for (i = 2; i <= 40; i++) {
    mpf_div_ui(x, x, i);
    mpf_add(s, s, x);
    gmp_printf("%4d: %.*Ff %.*Ff\n", i, ndigits, x, ndigits, s);
  }
  mpf_clear(s);
  mpf_clear(x);
  return 0;
}
