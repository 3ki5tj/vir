(* compute the ring contributions to the virial coefficients
   Usage
    math < ring.ma [dmin] [dmax] [nmin] [nmax] [prec] [output]
 *)
dmin = 2;
dmax = 50;
nmin = 1;
nmax = 64;
prec = 20;
fnls = "Bring.dat";



(* the contribution of Mayer ring diagrams to the virial series
   Marshall Luban, Asher Baram J. Chem. Phys. 76(6) 3233 1982 *)
Bnr[d_, n_, prec_: 20] := (-1)^(n - 1) Gamma[d/2]^(n - 2) *
  (n - 1)/n d^(n - 1) 2^(n d/2 - d) *
  NIntegrate[k^(d - 1) (BesselJ[d/2, k]/k^(d/2))^n, {k, 0, Infinity},
             PrecisionGoal -> prec, WorkingPrecision -> prec*5];



(* construct a list of virial series *)
mkBls[dmin_, dmax_, nmin_, nmax_, prec_] := Module[{n, d, x},
  (* the first column is reserved for the dimension *)
  ls = Table[ d,
              {d, dmin, dmax},
              {n, nmin, nmax + 1} ];
  For [d = dmin, d <= dmax, d++,
    For [ n = nmin, n <= nmax, n++,
      If [ n <= 2,
        x = 1,
        x = N[Bnr[d, n, prec], prec]
      ];
      ls[[d - dmin + 1, n - nmin + 2]] = x;
      Print["d ", d, ", n ", n, ", B ", (x // InputForm)];
    ];
  ];
  ls
];



(* main function starts here *)
argc = Length[$CommandLine];
argid = 2;
If [ argc >= argid, dmin = ToExpression[ $CommandLine[[argid]] ] ];
argid += 1;
If [ argc >= argid, dmax = ToExpression[ $CommandLine[[argid]] ] ];
argid += 1;
If [ argc >= argid, nmin = ToExpression[ $CommandLine[[argid]] ] ];
argid += 1;
If [ argc >= argid, nmax = ToExpression[ $CommandLine[[argid]] ] ];
argid += 1;
If [ argc >= argid, prec = ToExpression[ $CommandLine[[argid]] ] ];
argid += 1;
If [ argc >= argid, fnls = $CommandLine[[argid]] ];

Print["dim (", dmin, ", ", dmax, "), n (", nmin, ", ", nmax, ", ",
      "prec ", prec, ", file ", fnls];
ls = mkBls[dmin, dmax, nmin, nmax, prec];
Print["output saved to ", fnls];
Export[fnls, ls];

