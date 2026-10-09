\\ Units of the field $FIELD (e.g. L11): prove that U = <-1, fu_1, ..., fu_11> (from the GRH-conditional
\\ bnf) is a finite-index, 2-saturated subgroup of the full unit group E.
\\
\\ 1. Each fu_i is an algebraic integer of norm +-1 (exact check).
\\ 2. Quadratic characters: for degree-1 primes P = (p, t - r) not dividing the
\\    index or 2, map u -> Legendre(u mod P / p). The 12 x k matrix over F_2 of
\\    these characters for (-1, fu_1..fu_11) has rank 12.
\\ Rank 12 implies -1, fu_1, ..., fu_11 are independent in L*/L*^2, hence
\\ (r1 + r2 - 1 units; torsion of L is {+-1} since L has a real place)
\\ U has finite index in E, and no element of U \ U^2 is a square in L, so
\\ [E : U] is odd.
\\ Output: results/$FIELD/units.out (last line "STATUS OK" on success).

\\ F_2-rank of the quadratic-character matrix of gens (columns: degree-1 primes
\\ above p > pstart not dividing the index), stopping once rank = #gens or after
\\ maxcols columns. Returns [rank, columns used, last p].
charrank(nf, gens, pstart = 1000, maxcols = 400) =
{
  my(M = matrix(#gens, 0), k = 0, rk = 0, p = pstart);
  while (rk < #gens && k < maxcols,
    p = nextprime(p + 1);
    if (nf.index % p == 0, next);
    my(L = idealprimedec(nf, p));
    for (j = 1, #L,
      my(P = L[j]);
      if (P.f != 1 || P.e != 1, next);
      my(mp = nfmodprinit(nf, P), col);
      col = vector(#gens, i, my(v = nfmodpr(nf, gens[i], mp)); if (v == 0, error("element vanishes mod P")); if (issquare(v), 0, 1))~;
      M = concat(M, col); k++);
    rk = matrank(M * Mod(1, 2)));
  [rk, k, p];
}

read("certify/fields.gp");
outpath() = my(o = getenv("UNITS_OUT")); if (o, o, Str(fielddir(), "/units.out"));

units_main() =
{
  my(B = read(bnfpath()), nf, fu, gens, cr, f);
  nf = B.nf;
  fu = B.fu;
  if (#fu != nf.r1 + nf.r2 - 1, error("wrong number of fundamental units"));
  if (B.tu[1] != 2, error("unexpected torsion"));
  for (i = 1, #fu,
    my(u = nfalgtobasis(nf, fu[i]));
    if (denominator(u) != 1, error("unit ", i, " not integral"));
    if (abs(nfeltnorm(nf, u)) != 1, error("unit ", i, " norm != +-1")));
  gens = concat([nfalgtobasis(nf, -1)], apply(x -> nfalgtobasis(nf, x), fu));
  cr = charrank(nf, gens);
  if (cr[1] != #gens, error("F_2 rank ", cr[1], " < ", #gens));
  f = fileopen(outpath(), "w");
  filewrite(f, Str("unit_rank ", #fu, " torsion ", B.tu[1]));
  filewrite(f, Str("all_units_integral_norm_pm1 1"));
  filewrite(f, Str("character_columns ", cr[2], " up_to_p ", cr[3], " F2_rank ", cr[1]));
  filewrite(f, "units_2_saturated 1");
  filewrite(f, "STATUS OK");
  fileclose(f);
}

if (!getenv("UNITS_LIBRARY_ONLY"), units_main(); quit);
