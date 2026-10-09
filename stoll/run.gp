\\ Run Stoll's criterion for x^5 + y^5 = z^l, l = $STOLL_L, and write
\\ results/stoll/l<l>.out (last line "STATUS OK" when the run completed; the
\\ verdict is the PASS line). Uses the bnf at $STOLL_BNF (e.g. a certified
\\ results/<field>/field.bnf) when given, otherwise a fresh GRH-conditional bnfinit.
read("stoll/criterion.gp");

stoll_main() =
{
  my(l = eval(getenv("STOLL_L")), path = getenv("STOLL_BNF"), B, src, R, t0 = getabstime(), f);
  if (path, B = read(path); src = path, B = bnfinit(t^l - 2, 1); src = "bnfinit (GRH-conditional)");
  if (B.nf.pol != t^l - 2, error("bnf polynomial mismatch"));
  R = stollcriterion(l, B);
  f = fileopen(Str("results/stoll/l", l, ".out"), "w");
  filewrite(f, Str("l ", l));
  filewrite(f, Str("pari_version ", version()));
  filewrite(f, Str("bnf_source ", src));
  filewrite(f, Str("class_number ", B.no));
  foreach(["dim_LS2", "LS2_basis_verified", "dim_local5_image", "dim_S", "injective_at_2", "Zprime_hits", "mu11_in_S", "sigma_eq_conj86"], k,
    filewrite(f, Str(k, " ", mapget(R, k))));
  filewrite(f, Str("PASS ", mapget(R, "PASS")));
  filewrite(f, Str("ms ", getabstime() - t0));
  filewrite(f, "STATUS OK");
  fileclose(f);
}

stoll_main();
quit;
