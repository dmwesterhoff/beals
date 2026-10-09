SetClassGroupBounds("GRH");
SetVerbose("Selmer", 2);
K5<t> := NumberField(PolynomialRing(Rationals())![-1,1,1]);
_<x> := PolynomialRing(K5);
f := (4*t+12)/5*x^11 - (4*t+7)/5;
J := Jacobian(HyperellipticCurve(f));
S, mS := TwoSelmerGroup(J);
print "dim Sel2 =", Ngens(S);
print "codomain:", Codomain(mS);
