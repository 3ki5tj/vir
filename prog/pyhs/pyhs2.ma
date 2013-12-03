(* Exact solution of the PY equation for the d-dimensional
   hard-sphere system in by Groebner basis
   Copyright (c) Cheng Zhang 2013

  Example:

    math < pyhs.ma 11 c

  This computes the polynomial from the compressibility route for the
  11-dimensional system. The last letter can be `p' for the pressure
  route.  If the last letter is missing, it means both routes.
  The script works best for a dimensional <= 15.

  The analytical solution is based on a paper by E. Leutheusser,
  ``Exact solution of the Percus-Yevick equation for a hard-core fluid
    in odd dimensions'' Physica 127A (1984) 667-676 *)
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



(* get equations *)
geteqsQ[k_, t_] := Module[{Q, Qr, dQ0, i, n, r, r1, eqs, eq1, eq2, eq, lam},
  (* define the function Q(r) and the coefficients *)
  Q = Table[ ToExpression["Q" <> ToString[i]], {i, 0, k} ];
  Qr[r_] := Sum[ Q[[i + 1]] (r - 1)^(i + k), {i, 0, k} ];
  (* derivative of Q(r) at r = 0 *)
  dQ0[l_] := (D[Qr[r], {r, l}] /. {r-> 0});

  lam = (2 k + 1)!! 2^(2 k) t; (* packing fraction *)

  (* the first linear equation Eq. (21a) *)
  eq1 = Numerator[ Together[
    -(-1)^k - k! 2^k Q[[k + 1]] +  lam Sum[(-1)^i Q[[i + 1]]/(k + i + 1), {i, 0, k}]
  ] ];
  Qq0 = (-1)^(k + 1) k! 2^k Q[[-1]];
  eqs = { eq1 };

  (* the second linear equation Eq. (21b) *)
  If [ k >= 1,
    eq2 = Numerator[ Together[
      -(-1)^k - (k - 1)! 2^(k - 1) Q[[k]] + lam Sum[(-1)^i Q[[i + 1]]/(k + i + 2), {i, 0, k}]
    ] ];
    eqs = Append[eqs, eq2];
  ];

  (* the k-1 quadratic equations Eq. (22) *)
  For [ n = 0, n < k - 1, n++,
    eq = -dQ0[2 n + 1] + (lam/2)*(-1)^(n + 1) * (dQ0[n]^2)
       -lam * Sum[(-1)^nu dQ0[nu] * dQ0[2*n - nu], {nu, 0, n - 1}];
    eqs = Append[eqs, eq];
  ];
  {eqs, Q}
];



(* compute the polynomial of X = 1/sqrt(kappa) and t *)
numXt2[eqs_, X_, t_, Q_, pref_, nt_, dt_, nXmax_, i0_ : 0, i1_ : 100000] := Module[
  {tv, ls = {}, ls1, p1, val, v0, vr, i, j, s, sr, sref = 1, ex, pr, fn},

  ex = Exponent[pref, X];
  pr = Coefficient[pref, X^ex];
  fn = "ls" <> ToString[X] <> ToString[t] <> ToString[nt] <> ".txt";

  Print["total ", nt, " values, i ", i0, " - ", i1];

  (* compute the Groebner basis for each density t *)
  For [ tv = (-Ceiling[nt/2] + i0) * dt; ii = i0,
        Length[ls] < nt && ii < i1,
        tv += dt; ii += 1,
    If [ tv == 0 || tv == 1, Continue[]; ];
    Print["trying t = ", tv, ", i ", ii, " in [", i0, ", ", i1, ")"];
    MemoryConstrained[
      val = GroebnerBasis[eqs /. {t -> tv}, {X}, Q],
    2000000000];
    If [ Length[val] > 0, val = val[[1]], Continue[] ];
    v0 = Coefficient[ val, X^ex ];
    If [ v0 == 0, Continue[]; ];
    ls = Append[ls, {tv, val, v0}];
    If [ nt >= 64,
      Print["t ", (tv // InputForm), ", i ", ii, "/", nt, ", ", fn];
      xsave[fn, {tv, val}, True];
    ];
  ];
  (* Print[ls // InputForm]; *)
  If [ Length[ls] < nt && ii >= i0,
    Print[ "Collected partial list of ", i1 - i0, " items, exiting" ];
    Exit[];
  ];

  (* compute the reference scaling sref *)
  s = Table[1, {i, Length[ls]} ];
  For [ i = 1, i <= Length[ls], i++,
    v0 = ls[[i, 3]];
    vr = sref pr /. {t -> ls[[i, 1]]};
    s[[i]] = vr / v0;
    sr = LCM[vr, v0] / Abs[vr]; (* compute the reference scaling *)
    If [ sr != 1, (* need to scale the reference *)
      Print["multiplying ", sr, " to ", sref];
      sref *= sr;
      For [ j = 1, j <= i, j++, s[[j]] *= sr; ];
    ];
  ];

  (* apply the scaling to list items *)
  ls = Table[{ls[[i, 1]], ls[[i, 2]] * s[[i]]}, {i, Length[ls]}];
  (* Print[ls // InputForm]; *)

  For[p = 0; i = 0, i <= nXmax, i++,
    (* list of coefficients of X^i at different t's *)
    ls1 = Table[{ls[[j, 1]], Coefficient[ls[[j, 2]], X, i]}, {j, Length[ls]}];
    p += Factor[ InterpolatingPolynomial[ls1, t] ] X^i;
  ];
  If [ByteCount[p] < exprmax, Print[p // InputForm]];
  p
];



(* solve the first two Q's from the first two linear equations
   this proves to be time-saving in constructing Grobner bases *)
trimc[eqs_, Q_, i0_ : 1, i1_ : 2] := Module[{sub, eqs1},
  (* Print[Table[eqs[[j]] == 0, {j, 1, i1 - i0 + 1}] // InputForm]; *)
  sub = Solve[Table[eqs[[j]] == 0, {j, 1, i1 - i0 + 1}], Take[Q, {i0, i1}] ][[1]];
  (* Print[sub // InputForm]; *)
  (* replace Q[[i0..i1]] by the solution *)
  eqs1 = Drop[eqs, {1, i1 - i0 + 1}] /. sub;
  eqs1 = Table[ FactorTermsList[ Numerator[ Together[ eqs1[[j]] ] ] ][[2]],
  {j, Length[eqs1]} ];
  {eqs1, Drop[Q, {i0, i1}]}
];



(* compute virial coefficients from the exact polynomial *)
vircoef[d_, poly_, nmax_, X_, X0_, t_, ch_, prec_ : 20] := Module[
  {l, xx, fn, app, ser, bk, bkrat},
  fn = "vir" <> ToString[d] <> ch <> ".txt";
  (* Given poly(X, t) = 0, we need X = X0 + X1 t + X2 t^2 + ...
     We start with the contant term X = X0,
      then if poly(X0, t) = p0 + p1 t + ..., p0 must be zero
      i.e., poly(X0, 0) = 0
     Next, we try X = X0 + X1 t, and X1 is determined such that
      in the expansion poly(X0 + X1 t, t) = p0 + p1 t + ...,
      we have p1 == 0.
      Note, p0 == 0 is guaranteed by the previous step
     We can repeat the previous step and try X = X0 + X1 t + X2 t^2,
      and determine X2, etc. *)
  If [ (poly /. {X -> X0, t -> 0}) != 0,
    Print[ "bad initial value X = ", X0];
    Exit[];
  ];
  For [ app = X0; l = 1, l <= nmax, l++,
    app += xx t^l;
    (* Print[Collect[(poly /. {X -> app}), t] // InputForm]; *)
    app = app /. Solve[
      Limit[(poly /. {X -> app}) / t^l, t -> 0] == 0,
    xx][[1]];
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



(* main function starts here *)

nmax = 20; (* number of terms in the virial expansion *)
prec = 20; (* number of digits for virial coefficients *)

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

i0 = 0;
If [ Length[ $CommandLine ] >= 4,
  i0 = ToExpression[ $CommandLine[[4]] ];
];

i1 = 1000000;
If [ Length[ $CommandLine ] >= 5,
  i1 = ToExpression[ $CommandLine[[5]] ];
];
Print["start from i0 ", i0, ", ends at i1 ", i1];



(* 2. list basic equations *)
{eqs, Q} = geteqsQ[k, t];
If [ ByteCount[eqs] < exprmax, Print[(eqs // InputForm), " ", Q] ];
fneqs = "eqs" <> ToString[d] <> ".txt";
xsave[fneqs, eqs, False];
xsave[fneqs, Q, True];



(* 3. solve the compressibility equation *)
If [ StringCount[op, "c"] > 0,
  (* 3a. replace the last Q by X *)
  eqsc = eqs /. { Q[[-1]] -> X / ((-1)^(k + 1) k! 2^k) };
  Qc = Drop[Q, -1];
  xsave[fneqs, eqsc, True]; xsave[fneqs, Qc, True];

  (* 3b. remove the two linear equations *)
  {eqsc, Qc} = If [k >= 3, trimc[eqsc, Qc, 1, 2], {eqsc, Qc}];
  If [ ByteCount[eqsc] < exprmax, Print[(eqsc // InputForm), " ", Qc]; ];
  xsave[fneqs, eqsc, True]; xsave[fneqs, Qc, True];

  (* 3c. solve the equation of state f(X, t) == 0 *)
  tm = Timing[
    nXmax = If [k >= 2, 2^(k - 1), 1];
    dt = 1;
    nt = nXmax * 2 + 1;
    prc = (1 - t)^(2^k) X^(2^(k-1));
    Xt = numXt2[eqsc, X, t, Qc, prc, nt, dt, nXmax, i0, i1];
  ][[1]];
  Print[ (If [d <= 9, Xt, X^Exponent[Xt, X]] // InputForm), ", time ", tm];
  xsave["Xt" <> ToString[d] <> "c.txt", Xt, False, True];

  (* 3d. compute the virial coefficients *)
  vircoef[d, Xt, nmax, X, 1, t, "c", prec];
];


(* 4. pressure route *)
If [ StringCount[op, "p"] > 0,
  (* 4a. replace the first Q by g(1)
    g(1) = (-)^(k+1) k! Q[[1]], Z = 1 + 2^(d-1) t g(1) *)
  (* may leave a few fractions *)
  eqsp = eqs /. { Q[[1]] -> g / ((-1)^(k+1) k!) };
  Qp = Drop[Q, 1];
  (* Print[{eqsp, Qp} // InputForm]; *)
  xsave[fneqs, eqsp, True]; xsave[fneqs, Qp, True];

  (* 4b. remove the two linear equations *)
  {eqsp, Qp} = If [k >= 3, trimc[eqsp, Qp, 2, 3], {eqsp, Qp}];
  xsave[fneqs, eqsp, True]; xsave[fneqs, Qp, True];

  (* tm = Timing[ gt = Collect[ GroebnerBasis[eqsp, {g, t}, Qp][[1]], g, Factor]; ][[1]];
  Print[(gt // InputForm), ", time ", tm]; *)

  (* 4c. solve the equation of state f(g, t) == 0 *)
  tm = Timing[
    nXmax = If [k >= 2, 2^(k - 1), 1];
    dt = 1;
    nt = nXmax * 3;
    prc = (1 - t)^(2^k) t^(2^(k - 1) - 1) g^(2^(k - 1));
    gt = numXt2[eqsp, g, t, Qp, prc, nt, dt, nXmax, i0, i1];
  ][[1]];
  Print[ (If [d <= 9, gt, g^Exponent[gt, g]] // InputForm), ", time ", tm];
  (* convert to Zt *)
  Zt = FactorTermsList[ Numerator[ Together[ gt /. {g -> (Z - 1) / 2^(2k) / t} ] ] ][[2]];
  Zt = Collect[Zt, Z, Factor];
  xsave["Zt" <> ToString[d] <> "p.txt", Zt, False, True];

  (* 4d. compute the virial coefficients *)
  vircoef[d, Zt, nmax, Z, 1, t, "p", prec];
];
