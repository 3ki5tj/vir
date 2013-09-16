#include "mcutil.h"

int main(void)
{
  int n;
  for (n = 2; n <= 12; n++)
    printf("n %2d, Z3 %.14f Z4 %.14f B3 %.14f B4 %+.14f\n",
        n, Z3rat(n), Z4rat(n), B3rat(n), B4rat(n));
  return 0;
}
