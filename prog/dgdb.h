#ifndef DGDB_H__
#define DGDB_H__


/* compact data structure to read and write fb and nr */

enum { DGDB_FBTYPE_INT32 = 0x00, DGDB_FBTYPE_DBL = 0x10 };

/* the maximal fb is (n - 2)!, and the maximal nr is n!/(2 n) > (n - 2)!
 * a 32-bit signed integer can hold a number up to 2.1e9
 * a 32-bit unsigned integer can hold a number up to 4.2e9
 * the cutoff N is 14 because 13!/(2*13) = 2.4e8, 14!/(2*14) = 3.1e9
 * using 32-bit integers can save memory */
#if defined(N) && (N <= 14)
  typedef int32_t dgdb_fb_t;
  typedef uint32_t dgdb_nr_t;
  /* convert double to integer and save */
  #define DGDB_SAVEFB(a, b) \
    (a) = (dgdb_fb_t) (((b) < 0 ? ((b) - .5) : ((b) + .5)))
  #define DGDB_SAVENR(a, b) \
    (a) = (dgdb_nr_t) ((b) + .5)
  #define DGDB_FBTYPE DGDB_FBTYPE_INT32
  #define DGDB_FBPRI PRId32
  #define DGDB_FBSCN SCNd32
  #define DGDB_NRPRI PRIu32
  #define DGDB_NRSCN SCNu32
#else
  typedef double dgdb_fb_t;
  typedef double dgdb_nr_t;
  #define DGDB_SAVEFB(a, b) (a) = (b)
  #define DGDB_SAVENR(a, b) (a) = (b)
  #define DGDB_FBTYPE DGDB_FBTYPE_DBL
  #define DGDB_FBPRI  ".0f"
  #define DGDB_NRPRI  ".0f"
  #define DGDB_FBSCN  "lf"
  #define DGDB_NRSCN  "lf"
#endif


typedef struct {
  dgword_t c[DG_CWORDS];
  dgdb_fb_t fb;
#ifndef DG_NORING
  dgdb_nr_t nr;
#endif
} dgdbitem_t;

/* convenient formula to compute the size of the dgdbitem_t
 * this does not include padding introduced by sizeof() */
#define DGDBITEM_SIZE_NONR \
  (DG_CWORDS * sizeof(dgword_t) + sizeof(dgdb_fb_t))
#define DGDBITEM_SIZE_WITHNR \
  (DGDBITEM_SIZE_NONR + sizeof(dgdb_nr_t))

#ifdef DG_NORING
#define DGDBITEM_SIZE DGDBITEM_SIZE_NONR
#else
#define DGDBITEM_SIZE DGDBITEM_SIZE_WITHNR
#endif



#define dgdbitem_load(it, fp, fn, binary, fbtype, cwords, nr) \
  ((binary) ? dgdbitem_fread(it, fp, fn, fbtype, cwords, nr) \
            : dgdbitem_fscan(it, fp, fn, fbtype, cwords, nr))

#define dgdbitem_save(it, fp, fn, binary, fbtype, cwords, nr) \
  ((binary) ? dgdbitem_fwrite(it, fp, fn, fbtype, cwords, nr) \
            : dgdbitem_fprint(it, fp, fn, fbtype, cwords, nr))



/* read an item from a binary file */
INLINE int dgdbitem_fread(dgdbitem_t *it, FILE *fp, const char *fn,
    int fbtype, int cwords, int hasnr)
{
  if (cwords != fread(it->c, sizeof(dgword_t), cwords, fp)) {
    /* probably the end of file */
    //fprintf(stderr, "cannot read item c from %s\n", fn);
    return -1;
  }
  /* input type */
  if (fbtype == DGDB_FBTYPE_DBL) {
    double fb;
    if (1 != fread(&fb, sizeof(fb), 1, fp)) {
      fprintf(stderr, "cannot read double fb from %s\n", fn);
      return -2;
    }
    it->fb = (dgdb_fb_t) fb;
  } else if (fbtype == DGDB_FBTYPE_INT32) {
    int32_t fb;
    if (1 != fread(&fb, sizeof(fb), 1, fp)) {
      fprintf(stderr, "cannot read int32_t fb from %s\n", fn);
      return -2;
    }
    it->fb = (dgdb_fb_t) fb;
  } else {
    fprintf(stderr, "bad fbtype %#x\n", fbtype);
    exit(1);
  }
#ifndef DG_NORING
  if (hasnr) {
    if (fbtype == DGDB_FBTYPE_DBL) {
      double nr;
      if (1 != fread(&nr, sizeof(nr), 1, fp)) {
        fprintf(stderr, "cannot read double nr from %s\n", fn);
        return -3;
      }
      it->nr = (dgdb_nr_t) nr;
    } else if (fbtype == DGDB_FBTYPE_INT32) {
      uint32_t nr;
      if (1 != fread(&nr, sizeof(nr), 1, fp)) {
        fprintf(stderr, "cannot read uint32_t nr from %s\n", fn);
        return -3;
      }
      it->nr = (dgdb_nr_t) nr;
    } else {
      fprintf(stderr, "bad fbtype %#x\n", fbtype);
      exit(1);
    }
  } /* hasnr */
#endif
  return 0;
}



/* write an item to a binary file */
INLINE int dgdbitem_fwrite(dgdbitem_t *it, FILE *fp, const char *fn,
    int fbtype, int cwords, int hasnr)
{
  if (cwords != fwrite(it->c, sizeof(dgword_t), cwords, fp)) {
    fprintf(stderr, "cannot write item c to %s\n", fn);
    return -1;
  }
  if (fbtype == DGDB_FBTYPE_DBL) {
    double fb = (double) it->fb;
    if (1 != fwrite(&fb, sizeof(fb), 1, fp)) {
      fprintf(stderr, "cannot write double fb to %s\n", fn);
      return -2;
    }
  } else if (fbtype == DGDB_FBTYPE_INT32) {
    int32_t fb;
    if (DGDB_FBTYPE == fbtype) fb = (int32_t) it->fb;
    else fb = (int32_t) (it->fb >= 0 ? it->fb + .5 : it->fb - 0.5);
    if (1 != fwrite(&fb, sizeof(fb), 1, fp)) {
      fprintf(stderr, "cannot write int32_t fb to %s\n", fn);
      return -2;
    }
  } else {
    fprintf(stderr, "bad fbtype %x\n", fbtype);
    exit(1);
  }
#ifndef DG_NORING
  if (hasnr) {
    if (fbtype == DGDB_FBTYPE_DBL) {
      double nr = (double) it->nr;
      if (1 != fwrite(&nr, sizeof(nr), 1, fp)) {
        fprintf(stderr, "cannot write double nr to %s\n", fn);
        return -2;
      }
    } else if (fbtype == DGDB_FBTYPE_INT32) {
      int32_t nr;
      if (DGDB_FBTYPE == fbtype) nr = (int32_t) it->nr;
      else nr = (int32_t) (it->nr >= 0 ? it->nr + .5 : it->nr - 0.5);
      if (1 != fwrite(&nr, sizeof(nr), 1, fp)) {
        fprintf(stderr, "cannot write int32_t nr to %s\n", fn);
        return -2;
      }
    } else {
      fprintf(stderr, "bad fbtype %x\n", fbtype);
      exit(1);
    }
  } /* hasnr */
#endif
  return 0;
}



/* read an item from a text file */
INLINE int dgdbitem_fscan(dgdbitem_t *it, FILE *fp, const char *fn,
    int fbtype, int cwords, int hasnr)
{
  int i;

  if (feof(fp)) return 1;

  for (i = 0; i < cwords; i++) {
    if (1 != fscanf(fp, "%" DG_CSCN, &it->c[i])) {
      if (i > 0) /* i == 0 means the end of file */
        fprintf(stderr, "cannot read item c[%d] from %s\n", i, fn);
      return -1;
    }
  }
  if (fbtype == DGDB_FBTYPE_DBL) {
    double fb;
    if (1 != fscanf(fp, "%lf", &fb)) {
      fprintf(stderr, "cannot read fb from %s\n", fn);
      return -2;
    }
    it->fb = (dgdb_fb_t) fb;
  } else if (fbtype == DGDB_FBTYPE_INT32) {
    int32_t fb;
    if (1 != fscanf(fp, "%" SCNd32, &fb)) {
      fprintf(stderr, "cannot read fb from %s\n", fn);
      return -2;
    }
    it->fb = (dgdb_fb_t) fb;
  } else {
    fprintf(stderr, "bad fbtype %#x\n", fbtype);
    exit(1);
  }
#ifndef DG_NORING
  if (hasnr) {
    if (fbtype == DGDB_FBTYPE_DBL) {
      double nr;
      if (1 != fscanf(fp, "%lf", &nr)) {
        fprintf(stderr, "cannot read nr from %s\n", fn);
        return -3;
      }
      it->nr = (dgdb_nr_t) nr;
    } else if (fbtype == DGDB_FBTYPE_INT32) {
      uint32_t nr;
      if (1 != fscanf(fp, "%" SCNu32, &nr)) {
        fprintf(stderr, "cannot read uint32_t nr from %s\n", fn);
        return -3;
      }
      it->nr = (dgdb_nr_t) nr;
    } else {
      fprintf(stderr, "bad fbtype %#x\n", fbtype);
      exit(1);
    }
  } /* hasnr */
#endif
  return 0;
}



/* print an item to a text file
 * `fbtype' can be different from DGDB_FBTYPE
 * dgword_t must match the input file */
INLINE int dgdbitem_fprint(dgdbitem_t *it, FILE *fp, const char *fn,
    int fbtype, int cwords, int hasnr)
{
  int i;

  for (i = 0; i < cwords; i++)
    fprintf(fp, "%" DG_CPRI " ", it->c[i]);
  fprintf(fp, "%.0f", 1.*it->fb);
#ifndef DG_NORING
  if (hasnr)
    fprintf(fp, " %.0f", 1.*it->nr);
#endif
  fprintf(fp, "\n");
  return 0;
}



typedef struct {
  int version;
  int dim, n; /* dimension and order */
  int wordsz; /* sizeof(dgword_t) */
  int cwords; /* value of cwords */
  int fbtype; /* (DGDB_FBTYPE) 0x0: int32, 0x10: double */
  int hasnr;  /* has ring content */
  int itemsz; /* DGDBITEM_SIZE */
  int binary; /* input file is binary */
} dgdb_t;



/* create a database according to the default setting */
INLINE dgdb_t *dgdb_open(int dim, int n)
{
  dgdb_t *db;

  xnew(db, 1);
  db->version = 0;
  db->dim = dim;
  db->n = n;
  db->wordsz = (int) sizeof(dgword_t);
  db->cwords = (int) DG_CWORDS;
  db->fbtype = DGDB_FBTYPE;
#ifdef DG_NORING
  db->hasnr = 0;
#else
  db->hasnr = 1;
#endif
  db->itemsz = DGDBITEM_SIZE;
  return db;
}


#define dgdb_close(db) free(db)



#define dgdb_loadhead(db, fp, fn, dim, n, binary, check) \
  ((binary) ? dgdb_freadhead(db, fp, fn, dim, n, check) \
            : dgdb_fscanhead(db, fp, fn, dim, n, check))

#define dgdb_savehead(db, fp, fn, binary) \
  ((binary) ? dgdb_fwritehead(db, fp, fn) : \
              dgdb_fprinthead(db, fp, fn))



/* read database information from a binary file
 * if `check' is 0, just read in the information
 * if `check' is 1, check the information against `db'
 * */
INLINE int dgdb_freadhead(dgdb_t *db, FILE *fp, const char *fn,
    int dim, int n, int check)
{
  db->binary = 1;
  if (1 != fread(&db->version, sizeof(int), 1, fp)) {
    fprintf(stderr, "cannot read version from %s\n", fn);
    return -1;
  }
  if (1 != fread(&db->dim, sizeof(int), 1, fp)) {
    fprintf(stderr, "cannot read dimension from %s\n", fn);
    return -1;
  }
  if (check && db->dim != dim) {
    fprintf(stderr, "%s: bad dim %d vs %d\n", fn, db->dim, dim);
    return -1;
  }
  if (1 != fread(&db->n, sizeof(int), 1, fp)) {
    fprintf(stderr, "cannot read order from %s\n", fn);
    return -1;
  }
  if (1 != fread(&db->wordsz, sizeof(int), 1, fp)) {
    fprintf(stderr, "cannot read wordsz from %s\n", fn);
    return -1;
  }
  if (1 != fread(&db->cwords, sizeof(int), 1, fp)) {
    fprintf(stderr, "cannot read cwords from %s\n", fn);
    return -1;
  }
  if (1 != fread(&db->fbtype, sizeof(int), 1, fp)) {
    fprintf(stderr, "cannot read fbtype from %s", fn);
    return -1;
  }
  if (1 != fread(&db->hasnr, sizeof(int), 1, fp)) {
    fprintf(stderr, "cannot read hasnr from %s\n", fn);
    return -1;
  }
  if (1 != fread(&db->itemsz, sizeof(int), 1, fp)) {
    fprintf(stderr, "cannot read itemsz from %s\n", fn);
    return -1;
  }

  if (!check) return 0;

  if (db->n != n) {
    fprintf(stderr, "%s: bad order %d vs %d\n", fn, db->n, n);
    return -2;
  }
  if (db->wordsz != (int) sizeof(dgword_t)) {
    fprintf(stderr, "%s: bad wordsz %d vs %d\n",
        fn, db->wordsz, (int) sizeof(dgword_t));
    return -2;
  }
  if (db->cwords != DG_CWORDS) {
    fprintf(stderr, "%s: bad cwords %d vs %d\n",
        fn, db->cwords, DG_CWORDS);
    return -2;
  }
  if (db->fbtype != DGDB_FBTYPE) {
    fprintf(stderr, "%s: bad fbtype %d vs %d\n",
        fn, db->fbtype, DGDB_FBTYPE);
    return -2;
  }
  if (db->itemsz != DGDBITEM_SIZE) {
    fprintf(stderr, "%s: bad itemsz %d vs %d\n",
        fn, db->itemsz, (int) DGDBITEM_SIZE);
    return -2;
  }
  return 0;
}



/* write the header part of a binary database */
INLINE int dgdb_fwritehead(dgdb_t *db, FILE *fp, const char *fn)
{
  if (1 != fwrite(&db->version, sizeof(int), 1, fp)) goto ERR;
  if (1 != fwrite(&db->dim, sizeof(int), 1, fp)) goto ERR;
  if (1 != fwrite(&db->n, sizeof(int), 1, fp)) goto ERR;
  if (1 != fwrite(&db->wordsz, sizeof(int), 1, fp)) goto ERR;
  if (1 != fwrite(&db->cwords, sizeof(int), 1, fp)) goto ERR;
  if (1 != fwrite(&db->fbtype, sizeof(int), 1, fp)) goto ERR;
  if (1 != fwrite(&db->hasnr, sizeof(int), 1, fp)) goto ERR;
  if (1 != fwrite(&db->itemsz, sizeof(int), 1, fp)) goto ERR;
  return 0;
ERR:
  fprintf(stderr, "error occuried in writing %s\n", fn);
  return -1;
}



/* scan the header of a text database */
INLINE int dgdb_fscanhead(dgdb_t *db, FILE *fp, const char *fn,
    int dim, int n, int check)
{
  if (8 != fscanf(fp, "# V%d%d%d%d%d%d%d%d",
      &db->version, &db->dim, &db->n, &db->wordsz,
      &db->cwords, &db->fbtype, &db->hasnr, &db->itemsz)) {
    return -1;
  }

  if (!check) goto EXIT;

  if (db->dim != dim) {
    fprintf(stderr, "%s: dimension mismatch %d vs %d\n",
        fn, db->dim, dim);
    return -2;
  }
  if (db->n != n) {
    fprintf(stderr, "%s: order mismatch %d vs %d\n",
        fn, db->n, n);
    return -2;
  }
  if (db->wordsz != sizeof(dgword_t)) {
    fprintf(stderr, "%s: word size mismatch %d vs %d\n",
        fn, db->wordsz, (int) sizeof(dgword_t));
    return -2;
  }
  if (db->cwords != DG_CWORDS) {
    fprintf(stderr, "%s: cwords mismatch %d vs %d\n",
        fn, db->cwords, DG_CWORDS);
    return -2;
  }
  if (db->fbtype != DGDB_FBTYPE) {
    fprintf(stderr, "%s: fbtype mismatch %x vs %x\n",
        fn, db->fbtype, DGDB_FBTYPE);
    return -2;
  }
  if (db->itemsz != DGDBITEM_SIZE) {
    fprintf(stderr, "%s: item size mismatch %x vs %x\n",
        fn, db->itemsz, (unsigned) DGDBITEM_SIZE);
    return -2;
  }
EXIT:
  db->binary = 0;
  return 0;
}



/* print the header of a text database */
INLINE int dgdb_fprinthead(dgdb_t *db, FILE *fp, const char *fn)
{
  fprintf(fp, "# V%d %d %d %d %d %d %d %d\n",
      db->version, db->dim, db->n, db->wordsz,
      db->cwords, db->fbtype, db->hasnr, db->itemsz);
  return 0;
}



/* detect if the input file is binary */
INLINE int dgdb_detectbinary(const char *fn)
{
  FILE *fp;
  int bin;

  if ((fp = fopen(fn, "rb")) == NULL) {
    fprintf(stderr, "cannot open database %s\n", fn);
    return 0; /* assume text file */
  }
  /* text file starts with '#' */
  bin = (fgetc(fp) != '#');
  fclose(fp);
  return bin;
}


#endif /* DGDB_H__ */
