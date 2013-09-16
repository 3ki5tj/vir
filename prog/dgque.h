#ifndef DGQUE_H__
#define DGQUE_H__
/* quick history lookup list for recent values */
#include "dg.h"



#ifndef DGQUETYPE
#define DGQUETYPE double 
#endif

typedef struct {
  struct {
    DGQUETYPE val; /* e.g., the value of fb */
    code_t c[DG_NMAX/2]; /* code */
  } *arr;
  int nmax;
  int cnt; /* number of entries */
  int i; /* current point from 0 to cnt - 1 */
  int k; /* number of code_t */
} dgque_t;



INLINE dgque_t *dgque_open(int cnt, int k)
{
  dgque_t *q;

  xnew(q, 1);
  q->nmax = cnt;
  q->cnt = 0;
  q->k = k;
  xnew(q->arr, q->nmax);
  return q;
}



INLINE void dgque_close(dgque_t *q)
{
  free(q->arr);
  free(q);
}



INLINE int dgque_find(const dgque_t *q, const code_t *c, DGQUETYPE *x)
{
  int i, j, cnt = q->cnt, k = q->k;

  for (i = 0; i < cnt; i++) {
    for (j = 0; j < k; j++)
      if (q->arr[i].c[j] != c[j])
        break;
    if (j == k) {
      *x = q->arr[i].val;
      return i;
    }
  }
  return -1;
}



INLINE void dgque_add(dgque_t *q, code_t *c, DGQUETYPE x)
{
  int i, j, k = q->k;

  if (q->cnt < q->nmax) {
    i = q->cnt;
    q->arr[i].val = x;
    for (j = 0; j < k; j++) q->arr[i].c[j] = c[j];
    if ((q->cnt = i + 1) == q->nmax)
      q->i = 0; /* set point */
  } else { /* replace the q->i th point */
    i = q->i;
    q->arr[i].val = x;
    for (j = 0; j < k; j++) q->arr[i].c[j] = c[j];
    q->i = (i + 1) % q->nmax;
  }
}



#endif

