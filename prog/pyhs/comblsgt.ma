(* Exact solution of the PY equation for the d-dimensional
   hard-sphere system in by Groebner basis
   Copyright (c) Cheng Zhang 2013
   utility script for pyhs2.ma *)

Clear[Q, X, t, xsave, geteqsQ, numXt2, vircoef];


exprmax = 10000; (* maximal size for expression *)

xsave[fn_, xp_, append_: False, verbose_: False] := Module[{fp, s},
  If [ verbose,
    Print[ If [ append, "appending ", "writing " ], fn ];
  ];
  fp = If [ append, OpenAppend[fn], OpenWrite[fn] ];
  Write[fp, xp];
  Close[fp];
];



xload[fn_, verbose_: False] := Module[{fp, xp},
  If[verbose, Print["reading ", fn]];
  fp = OpenRead[fn];
  xp = Read[fp, Expression];
  Close[fp];
xp];



(* compute the polynomial of X and t *)
numXt2[ils_, X_, t_, pref_, nXmax_] := Module[
  {ls, tv, val, ls1, v0, vr, i, j, s, sr, sref = 1, ex, pr},

  ex = Exponent[pref, X];
  pr = Coefficient[pref, X^ex];
  (* make a copy of the input array, so it can be modified *)
  ls = Table[ils[[i, j]], {i, Length[ils]}, {j, 2}];

  (* compute the reference scaling sref *)
  s = Table[1, {i, Length[ls]} ];
  For [ i = 1, i <= Length[ls], i++,
    v0 = Coefficient[ ls[[i, 2]], X^ex ];
    vr = sref pr /. {t -> ls[[i, 1]]};
    s[[i]] = vr / v0;
    sr = LCM[vr, v0] / Abs[vr]; (* compute the reference scaling *)
    If [ sr != 1, (* need to scale the reference *)
      Print["multiplying ", sr, " to ", sref];
      sref *= sr;
      For [ j = 1, j <= i, j++, s[[j]] *= sr; ];
    ];
  ];
  Print["computed scaling, length ", Length[ls]];

  (* apply the scaling to list items *)
  ls = Table[{ls[[i, 1]], ls[[i, 2]] * s[[i]]}, {i, Length[ls]}];

  For [ p = 0; i = 0, i <= nXmax, i++,
    (* list of coefficients of X^i at different t's *)
    Print["computing coefficient of X^", i];
    ls1 = Table[{ls[[j, 1]], Coefficient[ls[[j, 2]], X, i]}, {j, Length[ls]}];
    p += Factor[ InterpolatingPolynomial[ls1, t] ] X^i;
  ];
  p
];



(* compute virial coefficients from the exact polynomial *)
vircoef[d_, poly_, nmax_, X_, X0_, t_, ch_, prec_ : 20] := Module[
  {l, xx, fn, app, ser, bk, bkrat},
  fn = "vir" <> ToString[d] <> ch <> ".txt";
  For [ app = X0; l = 1, l <= nmax, l++,
    app += xx t^l;
    (* Print[Collect[(poly /. {X -> app}), t] // InputForm]; *)
    app = app /. Solve[
      Limit[(poly /. {X -> app}) / t^l, t -> 0] == 0,
    xx][[1]];
  ];
  If [ ch == "c", app = Normal[ Series[ app^2, {t, 0, nmax} ] ]; ];
  ser = CoefficientList[app, t];
  bk = Table[ser[[k]] If[ch == "p", 1, 1/k] , {k, 1, Length[ser]}];
  xsave[fn, bk, False, True];

  (* Bk/ B2^(k-1) *)
  bkrat = Table[bk[[k]]/bk[[2]]^(k - 1), {k, Length[bk]}] ;
  Print["Bk/B2^(k-1) ", N[bkrat] // InputForm];
  xsave[fn, bkrat, True, False];
  xsave[fn, N[bkrat, prec], True];
];



(* main function starts here *)

nmax = 20; (* number of terms in the virial expansion *)
prec = 20; (* number of digits for virial coefficients *)

(* 1. handle input arguments *)
d = 19;
k = (d - 1)/2; (* an integer *)
op = "p";




(* 3. solve the compressibility equation *)

(* 3c. solve the equation of state f(g, t) == 0 *)
tm = Timing[
  nXmax = If [k >= 2, 2^(k - 1), 1];
  dt = 1;
  nt = nXmax * 3;
  prc = (1 - t)^(2^k) t^(2^(k - 1) - 1) g^(2^(k - 1));
  (* Remember to add "{" before the combined list
     replace all "}" by "},"
     and add "}" after the combined list *)
  ls = xload["lsgt768.txt", True];
  (* ls = xload["lsgt96.txt", True]; *)
  Print["list has ", Length[ls], " items, nXmax ", nXmax];
  gt = numXt2[ls, g, t, prc, nXmax];
][[1]];
Print[ (g^Exponent[gt, g] // InputForm), ", time ", tm];
(* convert to Zt *)
Zt = FactorTermsList[ Numerator[ Together[ gt /. {g -> (Z - 1) / 2^(2k) / t} ] ] ][[2]];
Zt = Collect[Zt, Z, Factor];
xsave["Zt" <> ToString[d] <> "p.txt", Zt, False, True];

(* 3d. compute the virial coefficients *)
vircoef[d, Zt, nmax, Z, 1, t, "p", prec];


