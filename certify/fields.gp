\\ Registry of fields certified by this repository, and shared path helpers.
\\ Every script selects its field with the environment variable FIELD.
\\   L11: Dahmen-Siksek's field for x^5 + y^5 = z^11 (algebra of C_{11,0} over Q(sqrt 5))
\\   L23: Stoll's field Q(2^(1/23)) for x^5 + y^5 = z^23 (curve 5y^2 = 4x^23 + 1)
fieldpol(name) =
{
  if (name == "L11", return(t^22 + 2*t^11 - 4));
  if (name == "L23", return(t^23 - 2));
  error("unknown FIELD: ", name);
}
fieldname() = my(f = getenv("FIELD")); if (!f, error("set FIELD (e.g. FIELD=L11)")); f;
fielddir() = Str("results/", fieldname());
bnfpath() = my(b = getenv("BNF")); if (b, b, Str(fielddir(), "/field.bnf"));
