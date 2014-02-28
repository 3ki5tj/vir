#include "mcutil.h"

int main(void)
{
  int n;
  for (n = 1; n <= 12; n++)
    printf("n %2d, Z3/V^2 %.15f Z4/V^3 %.15f B3/B2^2 %.15f B4/B2^3 %+.15f\n",
        n, Z3rat(n), Z4rat(n), B3rat(n), B4rat(n));
  return 0;
}
