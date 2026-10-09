\\ Independent cross-check of Phase 1: for every prime p in [P1_A, P1_B] and
\\ EVERY prime ideal P | p with N(P) <= P1_BOUND, find an explicit generator a
\\ (via the conditional bnf) and verify exactly that (a) = P as ideals.
\\ This uses neither the factor base nor the relation lattice: success means
\\ every such P is principal, unconditionally.
\\ Env: P1_A, P1_B, P1_BOUND, P1_OUT, optional BNF (default results/L11.bnf). Writes a count and a fingerprint
\\ (sum of generator coordinates times random-ish weights mod a 61-bit prime).
bnfpath() = my(b = getenv("BNF")); if (b, b, "results/L11.bnf");

gens_main() =
{
  my(a = eval(getenv("P1_A")), b = eval(getenv("P1_B")), bound = eval(getenv("P1_BOUND")));
  my(out = getenv("P1_OUT"), B = read(bnfpath()), nf, t0 = getabstime());
  my(nid = 0, fp = 0, q = 2^61 - 1, w, f);
  nf = B.nf;
  w = vector(poldegree(nf.pol), i, (1009^i) % q);
  forprime (p = a, b,
    my(L = idealprimedec(nf, p));
    for (j = 1, #L,
      my(P = L[j]);
      if (idealnorm(nf, P) > bound, next);
      my(r = bnfisprincipal(B, P, 3), g);
      if (#r[1] && r[1] != vector(#r[1])~, error("nonprincipal: ", P));
      g = r[2];
      if (idealhnf(nf, g) != idealhnf(nf, P), error("generator check failed: ", P));
      nid++;
      fp = (fp + sum(i = 1, #g, (g[i] % q) * w[i])) % q));
  f = fileopen(out, "w");
  filewrite(f, Str("range ", a, " ", b, " bound ", bound));
  filewrite(f, Str("ideals ", nid, " fingerprint ", fp, " ms ", getabstime() - t0));
  filewrite(f, "STATUS OK");
  fileclose(f);
}

gens_main();
quit;
