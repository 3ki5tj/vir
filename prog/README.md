Overview
=========


Directories
--------------

  o ./intg    Monte-Carlo integration
  o ./samp    Mayer sampling
  o ./java    Java code for Mayer sampling
  o ./ie      integral equation (numerical solution, C with MPFR, GSL), This is a separate project.
  o ./pyhs    analytical solution for the PY closure (Mathematica)
  o ./gaussf  Gaussian model (C with nauty, GMP)
  o .         basic routines of manipulating diagrams (graphs)
              serving ./samp

Files
------
dg.h        basic graph operations
dgmap.h     basic lookup table for n <= 8
dgutil.h    nonessential utility routines
dgcsep.h    detection of clique separators
dgcs.h      Ree-Hoover star content, which can be used to compute
            the total weight of a configuration of a hard-sphere system
dgrjw.h     the total weight of a configuration of a hard-sphere system
            alternative algorithm by R. J. Wheatley
            consider #define RJW32 1 on a 32-bit machine
dgrjwb.h    theoretically ``improved'' version of RJW, failed in practice (unused).
dgring.h    ring content
dgdb.h      database for graphs, used in the advanced lookup tables
            such as dgmapl.h and dghash.h
dgmapl.h    larger lookup table for n = 9, 10
dgaut.h     graph automorphism (used by dghash.h)
dghash.h    hash table
dgcryr.h    for computation of direct correlation function, c(r)
            and the cavity distribution function, y(r)


Utility programs
-----------------
bcstat.c    compute the probability of making a biconnected diagram
            from a chain configuration (as in the MC integration)

testrjw.c   test the speed


