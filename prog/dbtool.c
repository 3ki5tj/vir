#define ZCOM_PICK
#define ZCOM_ARGOPT
#include "zcom.h"

/* utilities to handle the database */
#include "dghash.h"
#include "dgrjw.h"
#include "dgring.h"



char *fninp = NULL;
char *fninp2 = NULL;
char *fninp3 = NULL;
char *fninp4 = NULL;
char *fnout = NULL;
int binary = -1;
int check = 1; /* default checking level */
int verbose = 0;



int hash_bits = 20;
int hash_blkmem = 0;
double hash_memmax = 3e10;
int auto_level = -1;
int hash_isoenum = 0;
int hash_isomax = 0;



/* handle arguments */
static void doargs(int argc, char **argv)
{
  argopt_t *ao;

  ao = argopt_open(0);
  ao->desc = "convert/check/merge database(s) for star and ring contents";
  argopt_add(ao, "-i", NULL, &fninp,    "input database");
  argopt_add(ao, "-j", NULL, &fninp2,   "second database");
  argopt_add(ao, "-k", NULL, &fninp3,   "third database");
  argopt_add(ao, "-l", NULL, &fninp4,   "fourth database");
  argopt_add(ao, "-o", NULL, &fnout,    "output database");
  argopt_add(ao, "-b", "%d", &binary,   "output format 1: binary, 0: text, -1: same as the input");
  argopt_add(ao, "-c", "%d", &check,    "checking, 0: disable, 1: biconnectivity, 2: fb and nr are correct");
  argopt_add(ao, "-v", "%d", &verbose,  "verbose");
  argopt_addhelp(ao, "-h");
  argopt_addhelp(ao, "--help");

  /* hash table parameters */
  argopt_add(ao, "--hash-bits",   "%d",   &hash_bits,   "number of bits in the key of the hash table");
  argopt_add(ao, "--hash-blkmem", "%d",   &hash_blkmem, "number of list links to cache each time (pass to the pooled memory allocator blkmem_new() to avoid memory fragmentation), 0: default");
  argopt_add(ao, "--hash-memmax", "%lf",  &hash_memmax, "maximal memory for the hash table");
  argopt_add(ao, "--auto-level",  "%d",   &auto_level,  "automorphism level, -1: canonical label, 0: no transformation, 1: degree sequence, 2 or 3: first automorphism in the searching tree");
  argopt_add(ao, "--hash-isoenum","%d",   &hash_isoenum,"enumerate isomorphic graphs after a new graph is found, 1: yes, 0: no, -1: default");
  argopt_add(ao, "--hash-isomax", "%d",   &hash_isomax, "maximal number of items to used in the above enumeration");
  argopt_parse(ao, argc, argv);
  argopt_close(ao);
}



dghash_t *hash;
dgdb_t *db;



/* check if an item is corrupted */
static int checkitem(dgdbitem_t *it, dg_t *g, int level,
    size_t id, size_t cnt)
{
  dg_decode(g, it->c);
  if ( !dg_biconnected(g) ) {
    fprintf(stderr, "id %.0f: graph is not biconnected\n", 1.*id);
    dg_printcode(it->c, hash->cwords);
    dg_print(g);
    return 1;
  }

  /* check if the fb and nr values are valid */
  if (level >= 2) {
    double fb;
    fb = dg_hsfb_mixed(g);
    if (fabs(fb - it->fb) > 0.1) {
      fprintf(stderr, "id %.0f: fb %g mismatch %g\n",
          1.*id, fb, 1.*it->fb);
      dg_printcode(it->c, hash->cwords);
      dg_print(g);
      return 2;
    }
#ifndef DG_NORING
    if (db->hasnr) {
      double nr = dg_nring_mixed(g);
      if (fabs(nr - it->nr) > 0.1) {
        fprintf(stderr, "id %.0f: nr %g mismatch %g\n",
            1.*id, nr, 1.*it->nr);
        dg_printcode(it->c, hash->cwords);
        dg_print(g);
        return 3;
      }
    }
#endif
    if ((cnt + 1) % 100 == 0)
      fprintf(stderr, "checked %.0f items   \r", 1.*cnt + 1);
  } else {
    if ((cnt + 1) % 1000000 == 0)
      fprintf(stderr, "checked %.0f items   \r", 1.*cnt + 1);
  }
  return 0;
}



/* load the database */
static int load(const char *fninp)
{
  int binary;
  FILE *fp;
  dgdbitem_t it[1] = {{{0}}};
  dgword_t hashid;
  dg_t *g;
  dgls_t *ls, *ls1;
  size_t cnt = 0, err = 0, id = 0, dup = 0;
  int ret, pos = 0;

  binary = dgdb_detectbinary(fninp);
  die_if ((db = dgdb_open(0, 0)) == NULL,
      "cannot open the database\n");
  die_if ((fp = fopen(fninp, binary ? "rb" : "r")) == NULL,
      "cannot open %s\n", fninp);
  die_if (dgdb_loadhead(db, fp, fninp, 0, 0, binary, 0) != 0,
      "cannot load header from %s\n", fninp);
  printf("checking level %d, D %d, n %d, binary %d\n", check, db->dim, db->n, binary);
  g = dg_open(db->n);
  hash = dghash_open(db->n, hash_bits, hash_blkmem, (size_t) hash_memmax,
      auto_level, hash_isoenum, hash_isomax);
  /* loop over entries of the hash table */
  id = 0;
  while (dgdbitem_load(it, fp, fninp, binary, db->fbtype,
                       hash->cwords, db->hasnr) == 0) {
    if ( check ) {
      if (checkitem(it, g, check, id, cnt) != 0) {
        fprintf(stderr, "database %s corrupted at id %.0f\n", fninp, 1.*id);
        break;
      }
    }
    DGHASH_GETID(hashid, it->c, hash->cwords, hash->bits);
    ls = hash->ls + hashid;
    if ( check ) {
      ls1 = dgls_find(ls, it->c, hash->cwords, &pos);
      if (ls1 != NULL) {
        if (dup == 0)
          fprintf(stderr, "duplication found!\n");
        dup++;
        continue;
      }
    }
    ret = dgls_additem(ls, it, hash, db->hasnr);
    if (ret != 0) {
      if (verbose)
        fprintf(stderr, "failed to add an item %.0f\n", 1.*cnt);
      err++;
    } else {
      cnt++;
    }
    id++;
  }
  printf("loaded %.0f entries, %.0f erors, %.0f dups\n", 1.*cnt, 1.*err, 1.*dup);
  fclose(fp);
  dg_close(g);
  return 0;
}



/* append the database `fninp' to the existing `db' */
static int append(const char *fninp)
{
  int binary;
  FILE *fp;
  dgdbitem_t it[1] = {{{0}}};
  dgword_t hashid;
  dgdb_t *db2;
  dg_t *g;
  dgls_t *ls, *ls1;
  size_t cnt = 0, err = 0, id = 0;
  int ret, pos = 0;

  binary = dgdb_detectbinary(fninp);
  die_if ((db2 = dgdb_open(0, 0)) == NULL,
      "cannot open the database\n");
  die_if ((fp = fopen(fninp, binary ? "rb" : "r")) == NULL,
      "cannot open %s\n", fninp);
  /* skip the header part of file */
  die_if (dgdb_loadhead(db2, fp, fninp, 0, 0, binary, 0) != 0,
      "cannot load header from %s\n", fninp);
  die_if (db->dim != db2->dim, "dimension mismatch %d vs %d\n", db->dim, db2->dim);
  die_if (db->n != db2->n, "order mismatch %d vs %d\n", db->n, db2->n);
  printf("checking level %d, D %d, n %d, binary %d\n",
      check, db->dim, db->n, binary);
  g = dg_open(db->n);
  /* loop over entries of the hash table */
  id = 0;
  while (dgdbitem_load(it, fp, fninp, binary, db2->fbtype,
                       hash->cwords, db->hasnr) == 0) {
    if ( check ) {
      if (checkitem(it, g, check, id, cnt) != 0) {
        fprintf(stderr, "database %s corrupted at item %.0f\n", fninp, 1.*id);
        break;
      }
    }
    DGHASH_GETID(hashid, it->c, hash->cwords, hash->bits);
    ls = hash->ls + hashid;
    ls1 = dgls_find(ls, it->c, hash->cwords, &pos);
    if (ls1 != NULL) { /* compare fb and nr values */
      //fprintf(stderr, "existing entry id %.0f, fb %g vs %g, nr %g vs %g\n",
      //    1.*id, 1.*ls1->arr[pos].fb, 1.*it->fb, 1.*ls1->arr[pos].nr, 1.*it->nr); getchar();
      if (ls1->arr[pos].fb != it->fb) {
        fprintf(stderr, "id %.0f, hashid %.0f, fb mismatch %g vs %g\n",
            1.*id, 1.*hashid, 1.*ls1->arr[pos].fb, 1.*it->fb);
        exit(1);
      }
#ifndef DG_NORING
      if (db2->hasnr && ls1->arr[pos].nr != it->nr) {
        fprintf(stderr, "id %.0f, hashid %.0f, nr mismatch %g vs %g\n",
            1.*id, 1.*hashid, 1.*ls1->arr[pos].nr, 1.*it->nr);
        exit(1);
      }
#endif
      continue;
    }
    ret = dgls_additem(ls, it, hash, db2->hasnr);
    if (ret != 0) {
      if (verbose)
        fprintf(stderr, "failed to add an item %.0f\n", 1.*cnt);
      err++;
    } else {
      cnt++;
    }
    id++;
  }
  printf("added %.0f entries, %.0f erors\n", 1.*cnt, 1.*err);
  fclose(fp);
  dg_close(g);
  dgdb_close(db2);
  return 0;
}



/* save the database */
static void save(const char *fnout, int binary)
{
  FILE *fp;
  int j, k;
  size_t i, cnt = 0;

  die_if ((fp = fopen(fnout, binary ? "wb" : "w")) == NULL,
      "cannot write %s\n", fnout);
  die_if (dgdb_savehead(db, fp, fnout, binary) != 0,
      "cannot write header to %s\n", fnout);
  /* loop over entries of the hash table */
  for (i = 0; i < hash->lsn; i++) {
    dgls_t *ls = hash->ls + i;
    if (ls->cnt == 0) continue;
    /* loop over links in the linked list */
    for (k = 0; ls != NULL; ls = ls->next, k++) {
      /* loop over items in the list */
      for (j = 0; j < ls->cnt; j++) {
        dgdbitem_save(ls->arr + j, fp, fnout, binary,
                      db->fbtype, hash->cwords, db->hasnr);
        //fprintf(stderr, "hashid %d, link %d, item %d/%d, %x\n",
        //    i, k, j, ls->cnt, ls->arr[j].c[0]); getchar();
      }
      cnt += ls->cnt;
    }
  }
  fprintf(stderr, "saved %.0f items to %s database %s\n",
      1.*cnt, binary ? "binary" : "text", fnout);
  fclose(fp);
}



int main(int argc, char **argv)
{
  doargs(argc, argv);
  load(fninp);
  if (fninp2) append(fninp2);
  if (fninp3) append(fninp3);
  if (fninp4) append(fninp4);
  dghash_printstat(hash, stdout);
  if (binary < 0) {
    if (strstr(fnout, ".txt") || strstr(fnout, ".tdb"))
      binary = 0;
    else if (strstr(fnout, ".bin") || strstr(fnout, ".bdb"))
      binary = 1;
    else
      binary = db->binary;
  }
  if (fnout != NULL)
    save(fnout, binary);
  return 0;
}

