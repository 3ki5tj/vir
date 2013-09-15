(* compute the ring contributions to the virial coefficients
   Copyright 2013 Cheng Zhang *)
nmax = 64
dmax = 50
prec = 20
fnls = "Bring.dat"



(* the contribution of Mayer ring diagrams to the virial series
   Marshall Luban, Asher Baram J. Chem. Phys. 76(6) 3233 1982 *)
Bnr[d_, n_, prec_: 20] := (-1)^(n - 1) Gamma[d/2]^(n - 2) *
  (n - 1)/n d^(n - 1) 2^(n d/2 - d) *
  NIntegrate[k^(d - 1) (BesselJ[d/2, k]/k^(d/2))^n, {k, 0, Infinity},
             PrecisionGoal -> prec, WorkingPrecision -> prec*5];



(* construct a list of virial series *)
mkBls[dmax_, nmax_, prec_] := Module[{n, d, x, d0 = 2},
  ls = Table[If[n == 1, d, 1], {d, d0, dmax}, {n, 1, nmax + 1}];
  For [d = d0, d <= dmax, d++,
    For [ n = 3, n <= nmax, n++,
      ls[[d - d0 + 1, n + 1]] = x = N[Bnr[d, n, prec], prec];
      Print["d ", d, ", n ", n, ", B ", (x // InputForm)];
    ];
  ];
  ls
];



(* main function starts here *)
argc = Length[$CommandLine];
If [ argc >= 2, dmax = ToExpression[ $CommandLine[[2]] ] ];
If [ argc >= 3, nmax = ToExpression[ $CommandLine[[3]] ] ];
If [ argc >= 4, prec = ToExpression[ $CommandLine[[4]] ] ];
If [ argc >= 5, fnls = ToExpression[ $CommandLine[[5]] ] ];

Print["dmax ", dmax, ", nmax ", nmax];
ls = mkBls[dmax, nmax, prec];
Print["output saved to ", fnls];
Export[fnls, ls];

