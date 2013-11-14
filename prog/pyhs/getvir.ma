(* compute the virial coefficients from the polynomials
   specifialized program to handle the output of pyhs2.ma
   Copyright (c) Cheng Zhang 2013

  Usage:

    math < getvir.ma [dim] [opt] [nmax] [prec]

  `dim' dimension, must be an odd integer
  `opt' can be any of `c' (compressibility), `p' (pressure),
  `nmax' is the maximal number of terms in the virial series
  `prec' is the precision
*)

Clear[xsave, xload];

xsave[fn_, xp_, append_: False, verbose_: False] := Module[{fp, s},
  If [ verbose,
    Print[ If [ append, "appending ", "writing " ], fn ];
  ];
  fp = If [ append, OpenAppend[fn], OpenWrite[fn] ];
  Write[fp, xp];
  Close[fp];
];

xload[fn_, verbose_: False] := Module[{fp, xp},
  If [ verbose, Print["reading ", fn] ];
  fp = OpenRead[fn];
  xp = Read[fp, Expression];
  Close[fp];
  xp
];

(* compute virial coefficients from the exact polynomial *)
vircoef[d_, poly_, nmax_, X_, X0_, t_, ch_, prec_ : 20] := Module[
  {l, xx, fn, app, ser, bk, bkrat},
  fn = "vir" <> ToString[d] <> ch <> ".txt";
  (* here we start with the contant term app = X0,
     then try app = X0 + X1 t, and determine X1,
     then try app = X0 + X1 t + X2 t^2, and determine X2, etc. *)
  For [ app = X0; l = 1, l <= nmax, l++,
    app += xx t^l;
    (* Print[Collect[(poly /. {X -> app}), t] // InputForm]; *)
    sol = Solve[
      Limit[(poly /. {X -> app}) / t^l, t -> 0] == 0,
    xx][[1]];
    app = app /. sol;
    Print["t^", l, ": ", (N[xx /. sol, 16] // InputForm)];
  ];
  (* in the compressibility route, the above procedure gives
     the expansion for sqrt(beta/(rho*kappa))
     so we square it to get beta/(rho*kappa) = d(beta*p)/d(rho) *)
  If [ ch == "c", app = Normal[ Series[ app^2, {t, 0, nmax} ] ]; ];
  ser = CoefficientList[app, t];
  (* in the compressibility route
     we now have the expansion of d(beta*p)/d(rho), integrate it
     over rho to get the virial expansion of beta*p *)
  bk = Table[ser[[k]] If[ch == "p", 1, 1/k] , {k, 1, Length[ser]}];
  xsave[fn, bk, False, True];

  (* Bk/ B2^(k-1) *)
  bkrat = Table[bk[[k]]/bk[[2]]^(k - 1), {k, Length[bk]}] ;
  Print["Bk/B2^(k-1) ", N[bkrat] // InputForm];
  xsave[fn, bkrat, True, False];
  xsave[fn, N[bkrat, prec], True];
];



nmax = 64; (* number of terms in the virial expansion *)
prec = 40; (* number of digits for virial coefficients *)


(* 1. handle input arguments *)
d = 11; (* must be an odd number *)
If [ Length[ $CommandLine ] >= 2,
  d = ToExpression[ $CommandLine[[2]] ];
];
If [ Mod[d, 2] == 0, Print["d ", d, " must be odd"]; Exit[]; ];
k = (d - 1)/2; (* an integer *)

op = "cp";
If [ Length[ $CommandLine ] >= 3,
  op = $CommandLine[[3]];
];
Print["compressibility ", StringCount[op, "c"] > 0,
      ", pressure ", StringCount[op, "p"] > 0];

If [ Length[ $CommandLine ] >= 4,
  nmax = ToExpression[ $CommandLine[[4]] ];
];

If [ Length[ $CommandLine ] >= 5,
  prec = ToExpression[ $CommandLine[[5]] ];
];

Print["number of terms ", nmax, ", internal precsion ", prec];

(* handle the compressibility-route equation *)
If [ StringCount[op, "c"] > 0,
  fn = "Xt" <> ToString[d] <> "c.txt";
  Xt = xload[fn];
  (* Print[ Xt // InputForm ]; *)
  Print["computing the compressibility-route virial coefficients..."];
  vircoef[d, Xt, nmax, X, 1, t, "c", prec];
];

(* handle the pressure-route equation *)
If [ StringCount[op, "p"] > 0,
  fn = "Zt" <> ToString[d] <> "p.txt";
  Zt = xload[fn];
  (* Print[ Zt // InputForm ]; *)
  Print["computing the pressure-route virial coefficients..."];
  vircoef[d, Zt, nmax, Z, 1, t, "p", prec];
];

