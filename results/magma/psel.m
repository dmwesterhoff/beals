SetClassGroupBounds("GRH");
P<t> := PolynomialRing(Rationals());
L<a> := NumberField(t^22 + 2*t^11 - 4);
O := MaximalOrder(L);
S := &cat[[Q[1] : Q in Decomposition(O, p)] : p in [2, 5, 11]];
print "#S =", #S;
time G, m := pSelmerGroup(2, Set(S));
print "dim L(S,2) =", Ngens(G);
print "expected (Cl=1, unit rank 11, torsion +-1) =", #S + 12;
