#ifndef DGDB_H__
#define DGDB_H__
/* database of the diagrams */



typedef struct {
  int n; /* number of entries in the list */
  code_t *code;
} dgls_t;



/* compute the hash value from the degree sequence */
INLINE int dg_gethash(dg_t *d, int *degseq)
{
  int sum, par, i, n = d->n;

  sum = dg_degseq(d, degseq);
  /* compute the number of partitions in the degrees */
  for (par = 0, i = 0; i < n - 1; i++)
    if (degseq[i] != degseq[i + 1]) par++;
  return sum * (n + 1) + par;
}



/* find the sequence `code' in the list `ls' */
INLINE int dgls_find(dgls_t *ls, code_t *code, int n)
{
  int i, j;

  for (i = 0; i < ls->n; i++) {
    for (j = 0; j < n; j++) {
      if (code[j] != ls->code[i * n + j])
        break;
    }
    if (j >= n) return i;
  }
  return -1;
}



/* add a code `code' into the list `ls' and the degree sequence */
INLINE void dgls_add(dgls_t *ls, unsigned *code, int n)
{
  int j;

  if (ls->n == 0) {
    xnew(ls->code, n);
  } else {
    xrenew(ls->code, (ls->n + 1) * n);
  }
  for (j = 0; j < n; j++)
    ls->code[ls->n * n + j] = code[j];
  ls->n++;
}



typedef struct {
  int n; /* number of particles */
  int npr; /* number of pairs */
  unsigned tot; /* number of entries */
  int sz; /* byte size of each entry */
  int nhash; /* number of hash lists */
  dgls_t *hash;

  /* local storage space */
  code_t *code; /* code */
  int *degseq; /* degree sequence */
  dg_t *d;
} dgdb_t;



/* add the graph to the database */
INLINE int dgdb_add(dgdb_t *db, dg_t *d)
{
  int add = 0;
  dgls_t *ls = db->hash + dg_gethash(d, db->degseq);

  if (dgls_find(ls, d->c, d->n) < 0) {
    dgls_add(ls, d->c, d->n);
    add = 1;
    db->tot++;
  }
  return add;
}



/* open the database */
INLINE dgdb_t *dgdb_open(int n)
{
  dgdb_t *db;
  int i, npr;

  xnew(db, 1);
  db->n = n;
  db->tot = 0;
  npr = db->n  * (db->n - 1) / 2;
  db->nhash = (npr + 1) * (db->n + 1);
  xnew(db->hash, db->nhash);
  for (i = 0; i < db->nhash; i++) {
    db->hash[i].n = 0;
    db->hash[i].code = NULL;
  }
  db->d = dg_open(n);
  xnew(db->degseq, n);
  return db;
}



/* close the database */
INLINE void dgdb_close(dgdb_t *db)
{
  int i;

  if (db->hash != NULL) {
    for (i = 0; i < db->nhash; i++)
      if (db->hash[i].code != NULL)
        free(db->hash[i].code);
    free(db->hash);
  }
  if (db->degseq) free(db->degseq);
  if (db->d) dg_close(db->d);
  free(db);
}



#endif

