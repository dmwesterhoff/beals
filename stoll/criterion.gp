\\ Stoll's "Selmer group Chabauty" criterion for C'_l : 5y^2 = 4x^l + 1
\\ (M. Stoll, Chabauty without the Mordell-Weil group, 2017, arXiv:1506.04286,
\\ Proposition 8.5 in the form of Corollary 8.7), reimplemented in PARI/GP.
\\
\\ If the criterion holds for an odd prime l >= 7, then C'_l(Q) = {oo, (1,1), (1,-1)}
\\ and (Dahmen-Siksek) x^5 + y^5 = z^l has only trivial primitive solutions.
\\
\\ Notation (as in Stoll): L = Q(lambda), lambda^l = 2 (variable t here);
\\ C'_l : y^2 = f(x), f = (4x^l + 1)/5, leading coefficient c = 4/5;
\\ theta = -lambda^(-2) is a root of f; descent map mu([a,b]) = (-c)^deg(a) a(theta).
\\
\\ The arithmetic of L enters only through L({5},2), i.e. the class group (2-part)
\\ and the units of L. With a GRH-conditional bnf the result is conditional on
\\ GRH unless Cl(L) and the units are certified (see certify/).

\\ ---------- 2-adic square classes via Hilbert symbols ----------

\\ F_2-vector of Hilbert symbols (u, b_j)_P for a list B of elements.
hilbvec(nf, u, B, P) = vector(#B, j, (1 - nfhilbert(nf, u, B[j], P)) / 2);

\\ A basis b_1..b_n of L_P^* / L_P^*2 (n = [L_P:Q_2] + 2) drawn from candidates.
\\ Let H be the Hilbert-pairing matrix of ALL candidates against each other. If
\\ rank H = n, the candidates span the n-dimensional space (H = X G X^T with G
\\ nondegenerate), and any n candidates whose rows of H are independent form a
\\ basis. Then u is a square in L_P iff hilbvec(u, basis) = 0.
sqclassbasis(nf, P) =
{
  my(n = P.e * P.f + 2, l = poldegree(nf.pol), cand, H, rows);
  cand = concat([t, -1, 3, 5], vector(2*l + 2, i, 1 + t^i));
  H = matrix(#cand, #cand, i, j, (1 - nfhilbert(nf, cand[i], cand[j], P)) / 2) * Mod(1, 2);
  if (matrank(H) != n, error("candidates do not span L_P^*/L_P^*2: rank ", matrank(H), " < ", n));
  rows = matindexrank(H)[1];
  vecextract(cand, rows);
}

\\ ---------- 2-adic arithmetic in L_2 = Q_2(lambda), done with exact rationals ----------
\\ Elements are exact t_POLMOD in t modulo t^l - 2 with rational coefficients;
\\ trunc2 replaces each coefficient by a 2-adically close rational (error < 2^-N)
\\ with power-of-2 denominator, which keeps sizes bounded.

trunc2c(a, N) =
{
  my(n = numerator(a), d = denominator(a), k = valuation(d, 2), d1 = d >> k);
  (lift(Mod(n, 2^(N + k)) / d1)) / 2^k;
}
trunc2(z, N) = Mod(Polrev(apply(a -> trunc2c(a, N), Vecrev(lift(z))), t), t^l_cur - 2);

\\ lambda-adic valuation of an exact element (oo for 0).
lval(z, l) =
{
  my(v = oo, c = Vecrev(lift(z)));
  for (i = 1, #c, if (c[i] != 0, v = min(v, valuation(c[i], 2) * l + (i - 1))));
  v;
}

\\ Square root in L_2 of z != 0, to 2-adic precision about N (as an exact element
\\ of L), or 0 if z is not a square in L_2. Digit-by-digit until the error is
\\ divisible by 4*lambda, then Newton x <- (x + u/x)/2 with truncation.
l2sqrt(z, l, N) =
{
  my(m, k, u, lam = Mod(t, t^l - 2), x = Mod(1, t^l - 2), e, v);
  l_cur = l;
  m = lval(z, l);
  if (m == oo, error("l2sqrt: zero"));
  if (m % 2, return(0));
  k = m / 2;
  u = trunc2(z / lam^m, N + 4);
  for (iter = 1, 4 * l + 10,
    e = u - x^2; v = lval(e, l);
    if (v > 2 * l, break);
    if (v == 2 * l || v % 2, return(0));
    x = x + lam^(v / 2));
  for (iter = 1, ceil(log(N * l + 1) / log(2)) + 5,
    x = trunc2((x + u / x) / 2, N + 8));
  if (lval(u - x^2, l) < (N - 4) * l, error("l2sqrt: precision not reached"));
  x * lam^k;
}

\\ ---------- fast 2-adic square classes ----------
\\ L_2 = Q_2(lambda) is totally ramified (e = l, residue field F_2). A basis of
\\ L_2^*/L_2^*2 is: lambda; 1 + lambda^w for odd w < 2l; one unramified class
\\ 1 + 4*eta. sq2(z) returns the coordinates [v(z) mod 2, a_1, a_3, ..., a_{2l-1}, b]
\\ by peeling off the levels of z / lambda^v(z) (digit-by-digit as in l2sqrt).
\\ Only z mod 4*lambda matters, so a few bits of 2-adic precision suffice.
sq2(z, l, N = 10) =
{
  my(lam = Mod(t, t^l - 2), v, eps, x = Mod(1, t^l - 2), a = vector(l), b = 0, e, w);
  l_cur = l;
  v = lval(z, l);
  if (v == oo, error("sq2: zero"));
  eps = trunc2(z / lam^v, N);
  for (iter = 1, 4 * l + 10,
    e = trunc2(eps - x^2, N); w = lval(e, l);
    if (w > 2 * l, break);
    if (w == 2 * l, b = 1; break);
    if (w % 2,
      a[(w + 1) / 2] = 1 - a[(w + 1) / 2];
      eps = trunc2(eps / (1 + lam^w), N),
      x = trunc2(x + lam^(w / 2), N)));
  concat([v % 2], concat(a, [b]));
}

\\ sq2 for an nf element (column or polmod) or a famat (homomorphism).
vec2fast(nf, u, l) = famatapply(z -> sq2(Mod(lift(nfbasistoalg(nf, z)), t^l - 2), l), u) % 2;

\\ ---------- halving: sigma, sigma' (Stoll, Lemma 8.4 and Corollary 5.4(2)) ----------

\\ theta-power basis <-> lambda basis (exact rationals).
thetabasis(l) =
{
  my(th = Mod(-t^(l - 2) / 2, t^l - 2));
  matrix(l, l, i, j, polcoef(lift(th^(j - 1)), i - 1));   \\ column j = theta^(j-1)
}

\\ For P1 = (xi1, eta1) in C'(Q_2) in the residue disk of (1,1) (x = 1 + 4u), the
\\ point i_(1,1)(P1) = [P1 - (1,1)] lies in 2J(Q_2). Following Corollary 5.4(2):
\\ s^2 = (x - xi1)(x - xi2) mod f; w = monic of minimal degree such that the
\\ reduction v of w*s mod f has degree <= g+1 and eta2 v(xi1) + eta1 v(xi2) = 0;
\\ then a half Q has mu(Q) = (-c)^deg(w) w(theta).
\\ The 2-adic objects are computed to 2-adic precision prec, then rounded to exact
\\ rationals; the linear system for w is solved exactly over Q and the defining
\\ conditions are re-checked 2-adically. Returns [sigma (global element of L), w].
sigmaclass(l, xi1, prec) =
{
  my(g = (l - 1) / 2, c = 4/5, f, eta1, eta2 = 1, xi2 = 1, th, gam, s, sq, Tinv, Ms, A, K, w, wpol, chk);
  f = (4*x^l + 1) / 5;
  eta1 = truncate(sqrt(subst(f, x, xi1) + O(2^(prec + 8))));
  if (eta1 % 4 != 1, eta1 = -eta1);
  th = Mod(-t^(l - 2) / 2, t^l - 2);
  gam = (th - xi1) * (th - xi2);
  s = l2sqrt(gam, l, prec);
  if (s == 0, error("i_(1,1)(P) not divisible by 2 in J(Q_2)"));
  sq = s;                                                          \\ exact rational approximation
  Tinv = thetabasis(l)^(-1);
  Ms = matrix(l, l, i, j, 0);
  for (j = 1, l, Ms[, j] = Tinv * Vecrev(lift(sq * th^(j - 1)), l)~);
  \\ g rows (theta-degrees g+2..2g of v, and the evaluation condition), g+1 unknowns w_0..w_g
  A = matrix(g, g + 1, i, j, 0);
  for (k = g + 2, l - 1, for (j = 1, g + 1, A[k - g - 1, j] = Ms[k + 1, j]));
  for (j = 1, g + 1, A[g, j] = sum(k = 0, l - 1, Ms[k + 1, j] * (eta2 * xi1^k + eta1 * xi2^k)));
  K = matker(A);
  if (#K != 1, error("halving: kernel dimension ", #K, " (expected 1)"));
  w = K[, 1];
  my(d = g + 1); while (w[d] == 0, d--);
  w = w / w[d];
  \\ only the 2-adic square class of w(theta) matters: round 2-adically
  w = apply(a -> trunc2c(a, prec), w);
  wpol = Polrev(w[1..d]~, x);
  \\ a-posteriori 2-adic check of the defining congruences with the true (2-adic) s
  chk = subst(wpol, x, th) * s;
  my(v = Tinv * Vecrev(lift(chk), l)~, minv = oo, ev);
  for (k = g + 2, l - 1, if (v[k + 1] != 0, minv = min(minv, valuation(v[k + 1], 2))));
  ev = sum(k = 0, l - 1, v[k + 1] * (eta2 * xi1^k + eta1 * xi2^k));
  if (ev != 0, minv = min(minv, valuation(ev, 2)));
  if (minv < prec / 2, error("halving: congruence check failed at 2-adic precision ", minv));
  [(-c)^poldegree(wpol) * subst(wpol, x, th), wpol];
}

\\ Stoll's Conjecture 8.6: sigma ~ prod_{n >= 1, v_2(n) odd} (1 + lambda^(l+n)).
\\ Terms with n > l + 1 are = 1 mod 4*lambda, hence squares, and can be dropped.
sigmaconj(l) =
{
  my(lam = Mod(t, t^l - 2), z = 1);
  for (n = 1, l + 2, if (valuation(n, 2) % 2 == 1, z *= 1 + lam^(l + n)));
  z;
}

\\ ---------- elements in compact form ----------
\\ An element is either an nf element or a factorization [gens, exps] (t_MAT with
\\ two columns, as returned by bnfunits); local square-class maps are
\\ homomorphisms, so they are evaluated factor by factor.
isfamat(u) = type(u) == "t_MAT";
famatapply(F, u) =
{
  if (!isfamat(u), return(F(u)));
  my(r = 0);
  for (i = 1, #u~, r += u[i, 2] * F(u[i, 1]));
  r % 2;
}

\\ 5-adic square-class vector: for each prime v | 5 (unramified, uniformizer 5):
\\ [v_v(u) mod 2, [u / 5^v_v(u) is a non-square in the residue field]].
vec5one(nf, u, S5, M5) =
{
  my(r = vector(2 * #S5));
  for (i = 1, #S5,
    my(e = nfeltval(nf, u, S5[i]), w = nfeltmul(nf, u, nfalgtobasis(nf, 5^(-e))), z = nfmodpr(nf, w, M5[i]));
    if (z == 0, error("vec5: not a unit after scaling"));
    r[2*i - 1] = e % 2; r[2*i] = if (issquare(z), 0, 1));
  r;
}
vec5(nf, u, S5, M5) = famatapply(x -> vec5one(nf, x, S5, M5), u) % 2;
vec2(nf, u, B2, P2) = famatapply(x -> hilbvec(nf, x, B2, P2), u) % 2;

\\ Is u (element or famat) in L(S,2)? Exponents of primes outside S must be even.
inLS2(nf, u, S) =
{
  my(fs = if (isfamat(u), u, Mat([u, 1])), Ps = List(), ex = List(), Skeys);
  Skeys = Set(apply(P -> Str(P), S));
  for (i = 1, #fs~,
    my(fa = idealfactor(nf, fs[i, 1]));
    for (k = 1, #fa~,
      my(key = Str(fa[k, 1]), pos = 0);
      for (m = 1, #Ps, if (Ps[m] == key, pos = m; break));
      if (pos, ex[pos] += fs[i, 2] * fa[k, 2], listput(Ps, key); listput(ex, fs[i, 2] * fa[k, 2]))));
  for (m = 1, #Ps, if (ex[m] % 2 && !setsearch(Skeys, Ps[m]), return(0)));
  1;
}

\\ F_2-rank of quadratic characters of a list of elements/famats at degree-1
\\ primes p > 1000 (stops at full rank or after 400 primes).
famatcharrank(nf, gens) =
{
  my(cols = List(), p = 1000, rk = 0, k = 0);
  while (rk < #gens && k < 400,
    p = nextprime(p + 1);
    if (nf.index % p == 0, next);
    foreach(idealprimedec(nf, p), P,
      if (P.f != 1 || P.e != 1, next);
      \\ skip P if it divides any factor of a compact representation (the
      \\ product is a unit/S-unit, but individual factors need not be)
      if (vecsum(apply(u -> if (isfamat(u), vecsum(apply(z -> nfeltval(nf, z, P) != 0, u[, 1]~)), nfeltval(nf, u, P) != 0), gens)), next);
      my(mp = nfmodprinit(nf, P), col);
      col = vector(#gens, i, famatapply(z -> my(v = nfmodpr(nf, z, mp)); if (v == 0, error("char: vanishes")); if (issquare(v), 0, 1), gens[i]));
      listput(cols, col); k++);
    if (#cols, rk = matrank(matrix(#cols, #gens, i, j, cols[i][j]) * Mod(1, 2))));
  rk;
}

\\ Does the 2-adic square class of the global element z lie in r(S)? (R from stollcriterion)
inimage2(R, z) =
{
  my(nf = mapget(R, "nf"), V2S = mapget(R, "V2S"), zv);
  zv = vec2fast(nf, nfalgtobasis(nf, lift(z)), poldegree(nf.pol))~ * Mod(1, 2);
  (zv == 0) || (matrank(concat(V2S, zv)) == matrank(V2S));
}

\\ ---------- the criterion ----------
\\ Returns a record (Map) with the verdict and all intermediate data.
\\ bnf: a bnf for t^l - 2 with h odd (GRH-conditional unless certified).
stollcriterion(l, bnf, prec = 120) =
{
  my(nf = bnf.nf, R = Map(), lam = Mod(t, t^l - 2), th, c = 4/5, F, S5, M5, P2, B2,
     gens, V5, V2, fac, hs, imgs, W, Sb, z1, z2, sg, sg2, Z, ok = 1);
  l_cur = l;
  if (!isprime(l) || l < 7, error("l must be a prime >= 7"));
  if (Mod(2, l^2)^(l - 1) == 1, error("Corollary 8.7 needs l^2 not dividing 2^(l-1) - 1"));
  if (bnf.no % 2 == 0, error("class number even: L({5},2) needs the 2-part of Cl"));
  th = Mod(-t^(l - 2) / 2, t^l - 2);
  F = 4*x^l + 1;
  S5 = idealprimedec(nf, 5); M5 = apply(P -> nfmodprinit(nf, P), S5);
  P2 = idealprimedec(nf, 2); if (#P2 != 1 || P2[1].e != l, error("2 not totally ramified"));
  P2 = P2[1]; B2 = 0;
  \\ L({5},2) basis: -1, fundamental units (compact), S-units for primes above 5
  my(U = bnfunits(bnf), su = bnfunits(bnf, S5));
  gens = concat([-1], vector(#su[1] - 1, i, su[1][i]));
  mapput(R, "dim_LS2", #gens);
  \\ Self-check that gens is a basis of L({5},2), independent of how PARI found
  \\ them: (a) each lies in L({5},2): every prime outside S5 occurs to an even
  \\ power; (b) they are independent mod squares (F_2-rank of quadratic
  \\ characters at degree-1 primes); (c) #gens = #S5 + r1 + r2, the dimension of
  \\ L({5},2) when the class number is odd.
  for (j = 1, #gens, if (!inLS2(nf, gens[j], S5), error("generator ", j, " not in L({5},2)")));
  my(cr = famatcharrank(nf, gens));
  if (cr != #gens || #gens != #S5 + nf.r1 + nf.r2, error("L({5},2) basis check failed: rank ", cr, ", #gens ", #gens));
  mapput(R, "LS2_basis_verified", 1);
  V5 = matrix(2 * #S5, #gens, i, j, 0); V2 = matrix(l + 2, #gens, i, j, 0);
  for (j = 1, #gens, V5[, j] = vec5(nf, gens[j], S5, M5)~; V2[, j] = vec2fast(nf, gens[j], l)~);
  \\ local image of J'(Q_5): images of the 2-torsion points [h, 0], h | F over Q_5
  fac = factorpadic(F, 5, 30)[, 1];
  hs = apply(h -> h / pollead(h), Vec(fac));
  if (#hs != #S5, error("factor/prime count mismatch at 5"));
  my(hint = apply(h -> Polrev(apply(a -> truncate(a), Vecrev(h)), x), hs), owner = vector(#hs));
  for (i = 1, #hint, for (k = 1, #S5,
      if (nfmodpr(nf, nfalgtobasis(nf, subst(hint[i], x, th)), M5[k]) == 0, owner[i] = k)));
  if (vecsort(owner) != vector(#S5, k, k), error("could not match Q_5-factors with primes above 5"));
  imgs = matrix(2 * #S5, #hs - 1, i, j, 0);
  for (j = 1, #hs - 1,
    my(h = hint[j], d = poldegree(h), rest, col = vector(2 * #S5));
    \\ F/h = 4 * (product of the other monic Q_5-factors), truncated 5-adically
    rest = Polrev(apply(a -> truncate(a), Vecrev(prod(i = 1, #hs, if (i == j, 1, hs[i])) * 4)), x);
    for (k = 1, #S5,
      my(e);
      if (owner[j] == k, e = (-c)^d * (-(1/5) * subst(rest, x, th)), e = (-c)^d * subst(h, x, th));
      my(vv = vec5one(nf, nfalgtobasis(nf, lift(e)), [S5[k]], [M5[k]]));
      col[2*k - 1] = vv[1]; col[2*k] = vv[2]);
    imgs[, j] = col~);
  if (matrank(imgs * Mod(1, 2)) != #hs - 1, error("images of J(Q_5)[2] dependent: order-4 points possible"));
  mapput(R, "dim_local5_image", #hs - 1);
  \\ S = { xi in L({5},2) : vec5(xi) in span(imgs) }
  my(Q = matker(imgs~ * Mod(1, 2))~);            \\ rows: linear forms vanishing on span(imgs)
  Sb = if (#Q, matker(Q * V5 * Mod(1, 2)), matid(#gens) * Mod(1, 2));
  mapput(R, "dim_S", #Sb);
  my(V2S = V2 * Mod(1, 2) * Sb);
  mapput(R, "V2S", V2S); mapput(R, "B2", B2); mapput(R, "P2", P2); mapput(R, "nf", nf);
  my(inj = matrank(V2S) == #Sb);
  mapput(R, "injective_at_2", inj); if (!inj, ok = 0);
  \\ the set Z' (Proposition 8.5)
  z1 = 1 + lam^(2*l - 1);
  z2 = 1 + lam^(l + 2) / (1 + lam^2);
  sg = sigmaclass(l, 5, prec)[1]; sg2 = sigmaclass(l, -3, prec)[1];
  Z = [z1, z2, sg, sg2];
  my(hits = vector(#Z, i, inimage2(R, Z[i])));
  mapput(R, "Zprime_hits", hits); if (vecsum(hits), ok = 0);
  \\ consistency: mu((1,1)) = (-c)(theta - 1) must lie in S (it is a global Selmer element)
  my(m11 = nfalgtobasis(nf, lift((-c) * (th - 1))), f11 = idealfactor(nf, m11), inLS = 1);
  for (i = 1, #f11~, if (f11[i, 2] % 2 && f11[i, 1].p != 5, inLS = 0));
  my(v11 = vec5(nf, m11, S5, M5)~ * Mod(1, 2));
  mapput(R, "mu11_in_S", inLS && (#Q == 0 || Q * v11 == 0));
  if (!mapget(R, "mu11_in_S"), error("consistency check failed: mu((1,1)) not in S"));
  mapput(R, "sigma_eq_conj86", sq2(sg, l) == sq2(sigmaconj(l), l));
  mapput(R, "PASS", ok);
  R;
}
