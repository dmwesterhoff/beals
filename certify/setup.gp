\\ Build the (GRH-conditional) bnf for L11 once and save it; prove that every
\\ prime in its factor base Vbase is principal by exhibiting a generator and
\\ checking ideal equality exactly.
\\ Output: results/L11.bnf, results/setup.out (last line "STATUS OK" on success).
\\ Run: gp -q -D parisizemax=2000000000 certify/setup.gp < /dev/null
read("certify/zimmert.gp");

setup_main() =
{
  my(g = t^22 + 2*t^11 - 4, B, nf, V, zb, Z, ZB, MB, n, f);
  B = bnfinit(g, 1);
  nf = B.nf;
  n = poldegree(g);
  if (B.no != 1, error("conditional class number is not 1"));
  \\ the integral basis (hence disc and the ring of integers) is proven correct
  if (nfcertify(nf) != [], error("nfcertify: maximal order not certified"));
  \\ writebin appends to an existing file: write a fresh temp file, then rename
  \\ atomically so concurrent readers never see a partial or doubled file.
  system("rm -f results/L11.bnf.tmp");
  writebin("results/L11.bnf.tmp", B);
  system("mv -f results/L11.bnf.tmp results/L11.bnf");

  zb = zimmertbest(nf.r1, nf.r2);
  Z = zimmertZ(nf.r1, nf.r2, zb[2], zimmertalpha(zb[2]));
  ZB = ceil(sqrt(abs(nf.disc)) * exp(-Z)) + 1;
  MB = ceil(2^(2*nf.r2) / Pi^nf.r2 * n! / n^n * sqrt(abs(nf.disc)));

  V = B[5];
  for (i = 1, #V,
    my(P = V[i], r = bnfisprincipal(B, P, 3)); \\ nf_GEN | nf_FORCE
    if (#r[1] && r[1] != vector(#r[1])~, error("nontrivial class: ", P));
    if (idealhnf(nf, r[2]) != idealhnf(nf, P), error("generator check failed: ", P)));

  f = fileopen("results/setup.out", "w");
  filewrite(f, Str("pari_version ", version()));
  filewrite(f, Str("polynomial ", g));
  filewrite(f, Str("disc ", nf.disc, " = ", factor(nf.disc)));
  filewrite(f, Str("signature ", nf.sign));
  filewrite(f, Str("zimmert_gamma ", zb[2], " Z ", Z));
  filewrite(f, Str("zimmert_bound ", ZB));
  filewrite(f, Str("primepi_zimmert_bound ", primepi(ZB)));
  filewrite(f, "maximal_order_certified 1");
  filewrite(f, Str("minkowski_bound ", MB));
  filewrite(f, Str("vbase_size ", #V, " max_norm ", vecmax(apply(P -> idealnorm(nf, P), V))));
  filewrite(f, "vbase_all_principal_with_verified_generators 1");
  filewrite(f, "STATUS OK");
  fileclose(f);
}

setup_main();
quit;
