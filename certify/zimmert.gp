\\ Zimmert's unconditional bound (Invent. Math. 62 (1981) 367-380, Satz 2 and proof).
\\ Every ideal class of a number field K (degree n = r1 + 2 r2, discriminant d)
\\ contains an integral ideal a with  log(|d|^(1/2) / N(a)) >= zimmertZ(r1, r2, gam, alp)
\\ for any real gam > alp > 0.
\\ Implemented from the a,b-form in the proof of Satz 2 (a = r1 + r2, b = r2,
\\ A = |d|^(1/2) pi^(-n/2)), which avoids relying on the simplified display.
zimmertZ(r1, r2, gam, alp) =
{
  my(a = r1 + r2, b = r2, n = r1 + 2*r2);
  a * (-psi((1+gam)/2) - lngamma(1/2+gam) + lngamma(1+gam))
  + b * (-psi(1+gam/2) - lngamma(1+gam) + lngamma(3/2+gam))
  - 2/(gam-alp)
  - log((1 + 1/alp) * (1 + 1/gam)^(-2) * (1 + 1/(2*gam-alp))^(-1))
  + n/2 * log(Pi);
}
\\ Same quantity from the printed statement of Satz 2 (cross-check of transcription).
zimmertZ_satz2(r1, r2, gam, alp) =
{
  r1 * (-psi((1+gam)/2) - lngamma(1/2+gam) + lngamma(1+gam) + log(Pi)/2)
  + r2 * (-2*psi(1+gam) + 2*log(2) + log(1/2+gam) + log(Pi))
  - 2/(gam-alp)
  - log((1 + 1/alp) * (1 + 1/gam)^(-2) * (1 + 1/(2*gam-alp))^(-1));
}
\\ Bemerkung 1: optimal alpha for fixed gamma.
zimmertalpha(gam) = gam - gam*(gam+1)/sqrt(1 + 3*gam + 3*gam^2);
\\ Best Z over a grid of gamma (any gamma > 0 is valid; the grid only affects sharpness).
zimmertbest(r1, r2) =
{
  my(best = -oo, bg = 0);
  for (i = 1, 4000, my(g = i/1000., z = zimmertZ(r1, r2, g, zimmertalpha(g)));
    if (z > best, best = z; bg = g));
  [best, bg];
}
