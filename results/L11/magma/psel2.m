SetClassGroupBounds("GRH");
P<t> := PolynomialRing(Rationals());
L<a> := NumberField(t^22 + 2*t^11 - 4);
O := MaximalOrder(L);
S := &cat[[Q[1] : Q in Decomposition(O, p)] : p in [2, 5, 11]];
G, m := pSelmerGroup(2, Set(S));
reps := [L!((G.i) @@ m) : i in [1..Ngens(G)]];
print "dim =", #reps;
// membership in L(S,2): every prime ideal outside S occurs to an even power in (r)
inLS2 := true;
for r in reps do
  for pe in Factorization(r*O) do
    if pe[1] notin S and IsOdd(pe[2]) then inLS2 := false; end if;
  end for;
end for;
print "all reps in L(S,2):", inLS2;
F2 := GF(2); cols := []; p := 1000; rk := 0;
while rk lt #reps do
  p := NextPrime(p);
  for QQ in Decomposition(O, p) do
    Q := QQ[1];
    if Degree(Q) ne 1 or QQ[2] ne 1 then continue; end if;
    k, mk := ResidueClassField(Q);
    Append(~cols, [ IsSquare(mk(r)) select F2!0 else F2!1 : r in reps ]);
  end for;
  if #cols gt 0 then rk := Rank(Matrix(F2, cols)); end if;
end while;
print "character columns =", #cols, " up to p =", p, " F2-rank =", rk;
