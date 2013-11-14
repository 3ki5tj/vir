(* Exact solution of the PY equation for the d-dimensional
   hard-sphere system in by Groebner basis
   Copyright (c) Cheng Zhang 2013

   The analytical solution is based on a paper by E. Leutheusser,
   ``Exact solution of the Percus-Yevick equation for a hard-core fluid
     in odd dimensions'' Physica 127A (1984) 667-676 *)
Clear[Q, X, t, xsave, geteqsQ, getprc, numXt, vircoef];


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



(* get the reference polynomial at Q[[-1]] = 0
   this special polynomial determines the radius of convergence
   it also allows us to use polynomial interpolation for the full solution *)
getprc[eqs_, Q_, t_, d_, prec_ : 20] := Module[{p, ls, fn},
  (* the last Q is proportional to Q~(q = 0) *)
  MemoryConstrained[
    p = GroebnerBasis[eqs /. {Q[[-1]] -> 0}, {t}, Drop[Q, -1] ][[1]],
  2000000000];
  p = Factor[p];
  fn = "tconv" <> ToString[d] <> ".txt";
  xsave[fn, p, False, True];
  Print[(p // InputForm)];
  ls = t /. NSolve[p, t, Reals, WorkingPrecision -> prec];
  ls = Sort[ls, Abs[#1] < Abs[#2] &]; (* sort by magnitude *)
  xsave[fn, ls[[1]], True, False];
  Print["tc = ", ls[[1]]];
  p
];



(* compute the polynomial of X and t *)
numXt[eqs_, X_, t_, Q_, pref_, nt_, dt_, texcl_, nXmax_] := Module[
  {tv, ls = {}, ls1, p1, val, v0, vr, i, j, s, sr, sref = 1, fn},

  fn = "ls" <> ToString[X] <> ToString[t] <> ToString[nt] <> ".txt";
  (* compute the Groebner basis for each t *)
  For [ tv = -Ceiling[nt/2] * dt, Length[ls] < nt, tv += dt,
    If [ Count[texcl, tv] > 0, Continue[];];
    MemoryConstrained[
      val = GroebnerBasis[eqs /. {t -> tv}, {X}, Q],
    2000000000];
    If [ Length[val] > 0, val = val[[1]], Continue[] ];
    ls = Append[ls, {tv, val}];
    If [ nt >= 64,
      Print["t ", (tv // InputForm), " ", fn];
      xsave[fn, {tv, val}, True];
    ];
  ];
  (* Print[ls // InputForm]; *)

  (* compute the reference scaling sref *)
  s = Table[1, {i, Length[ls]} ];
  For [ i = 1, i <= Length[ls], i++,
    v0 = ls[[i, 2]] /. {X -> 0};
    vr = sref pref /. {t -> ls[[i, 1]]};
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

  For[p = 0; i = 0, i <= nXmax, i++,
    (* list of coefficients of X^i at different t's *)
    ls1 = Table[{ls[[j, 1]], Coefficient[ls[[j, 2]], X, i]}, {j, Length[ls]}];
    p += Factor[ InterpolatingPolynomial[ls1, t] ] X^i;
  ];
  (* remove common factors *)
  p = Collect[ FactorTermsList[ p ][[2]], X, Factor];
  If [ByteCount[p] <= 10240, Print[p // InputForm]; ];
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
  For [ app = X0; l = 1, l <= nmax, l++,
    (* Print["hahaha, ", app, " ", X0, " X: ", X]; *)
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



(* 2. list basic equations *)
{eqs, Q} = geteqsQ[k, t];
Print[(eqs // InputForm), " ", Q];
fneqs = "eqs" <> ToString[d] <> ".txt";
xsave[fneqs, eqs, False];
xsave[fneqs, Q, True];


(* 3. solve the compressibility equation *)
If [ StringCount[op, "c"] > 0,
  (* 3a. radius of convergence; use (eqs, Q) instead of (eqsc, Qc),
     for the latter will slow down the solution *)
  tm = Timing[ prc = getprc[eqs, Q, t, d]; ][[1]];
  xsave[fneqs, eqs /. {Q[[-1]] -> 0}, True]; xsave[fneqs, Drop[Q, -1], True];
  Print["prc time ", tm];

  (* 3b. replace the last Q by X, and remove linear equations *)
  eqsc = eqs /. { Q[[-1]] -> X / ((-1)^(k + 1) k! 2^k) };
  Qc = Drop[Q, -1];
  xsave[fneqs, eqsc, True]; xsave[fneqs, Qc, True];

  {eqsc, Qc} = If [k >= 3, trimc[eqsc, Qc, 1, 2], {eqsc, Qc}];
  Print[(eqsc // InputForm), " ", Qc];
  xsave[fneqs, eqsc, True]; xsave[fneqs, Qc, True];

  (* 3c. solve the equation of state f(X, t) == 0 *)
  tm = Timing[
    nXmax = If [k >= 2, 2^(k - 1), 1];
    dt = 1;
    nt = nXmax * 2 + 1;
    Xt = numXt[eqsc, X, t, Qc, prc, nt, dt, {0}, nXmax];
  ][[1]];
  Print[ (If [d <= 9, Xt, X^Exponent[Xt, X]] // InputForm), ", time ", tm];
  xsave["Xt" <> ToString[d] <> "c.txt", Xt, False, True];

  (* 3d. compute the virial coefficients *)
  vircoef[d, Xt, nmax, X, 1, t, "c", prec];
];


(* 4. pressure route *)
If [ StringCount[op, "p"] > 0,
  (* 4a. compute the reference polynomial with g or Q[[1]] == 0 *)
  tm = Timing[
    gt0 = GroebnerBasis[eqs /. {Q[[1]] -> 0}, {t}, Drop[Q, 1] ][[1]];
  ][[1]];
  Print["zero polynomial ", (gt0 // InputForm), ", time ", tm];

  (* 4b.  replace the first Q by g(1), remove the linear equations,
    g(1) = (-)^(k+1) k! Q[[1]], Z = 1 + 2^(d-1) t g(1) *)
  (* may leave a few fractions *)
  eqsp = eqs /. { Q[[1]] -> g / ((-1)^(k+1) k!) };
  Qp = Drop[Q, 1];
  (* Print[{eqsp, Qp} // InputForm]; *)
  xsave[fneqs, eqsp, True]; xsave[fneqs, Qp, True];

  {eqsp, Qp} = If [k >= 3, trimc[eqsp, Qp, 2, 3], {eqsp, Qp}];
  xsave[fneqs, eqsp, True]; xsave[fneqs, Qp, True];

  (* tm = Timing[ gt = Collect[ GroebnerBasis[eqsp, {g, t}, Qp][[1]], g, Factor]; ][[1]];
  Print[(gt // InputForm), ", time ", tm]; *)

  (* 4c. solve the equation of state f(g, t) == 0 *)
  tm = Timing[
    nXmax = If [k >= 2, 2^(k - 1), 1];
    dt = 1;
    nt = nXmax * 3;
    gt = numXt[eqsp, g, t, Qp, gt0, nt, dt, {0, 1}, nXmax];
  ][[1]];
  Print[ (If [d <= 9, gt, g^Exponent[gt, g]] // InputForm), ", time ", tm];
  (* convert to Zt *)
  Zt = Collect[ Numerator[ Together[ gt /. {g -> (Z - 1) / 2^(2k) / t} ] ], Z, Factor];
  xsave["Zt" <> ToString[d] <> "p.txt", Zt, False, True];

  (* 4d. compute the virial coefficients *)
  vircoef[d, Zt, nmax, Z, 1, t, "p", prec];
];


(* 5. direct correlation function c(r), don't use, very slow *)
If [ StringCount[op, "C"] > 0,
  (* 5a. construct c(r) *)
  cr = -(2^k k! Q[[k + 1]])^2;
  dQ1[m_] := Q[[m - k + 1]] m!;
  For [ n = 0, n <= k, n++,
    (* Cn = C^(2k + 2n + 1)(0) / t *)
    Cn = (-1)^(k + n) dQ1[k + n]^2 + 2 * If [ 2*n <= k,
        Sum[(-1)^(nu + k) dQ1[k + 2 n - nu] dQ1[k + nu], {nu, 0, n -1}],
        Sum[(-1)^nu dQ1[2 k - nu] dQ1[nu + 2 n], {nu, 0, k - n - 1}] ];
    cr += (-1)^k n! (t 4^k (2 k + 1)!!) Cn / (2^k (n + k)! (2 n + 1)!) r^(2 n + 1);
  ];
  eqscr = Append[eqs, cr - c];
  tm = Timing[  crt = Collect[ GroebnerBasis[eqscr, {c, r, t}, Q][[1]], {c, r}, Factor ]; ][[1]];
  Print[(crt // InputForm), ", time ", tm];
];

