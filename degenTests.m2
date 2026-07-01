load "affGrSchubert.m2"

a =3;
e =5;
lam = (a,0);
mu = (a-1,1);
F = QQ;

varY = toList (y_(0,0,1)..y_(1,1,a-1));
varT = toList (t_1..t_(e-1));

R = F[varY,varT];

S = R[getSymbol "u"];
u = S_0;

diffFactor = if e <= 1 then 1_S else product apply(toList(1..e-1), i -> u - sub(t_i, S));

-- Matrix of Nnabla on F[u]^2/u^a with respect to the basis
-- (e1, e2, u e1, u e2, ..., u^(a-1)e1, u^(a-1)e2).
Nnabla = matrix table(2*a, 2*a, (row, col) -> (
    i := row % 2;
    s := row // 2;
    j := col % 2;
    r := col // 2;

    nPart := if s <= r or s >= a then 0_S else sub(y_(i,j,s-r), S);
    dPartPolynomial := r * diffFactor * u^r;
    dPart := if i == j then coefficient(u^s, dPartPolynomial) else 0_S;

    sub(nPart + dPart, R)
));


K = kernel Nnabla

rankAtMostOneLocus = muInput -> (
    -- Choose a coweight mu <= lam and pull back the Schubert chart data.
    (Rchart, Schart, uChart, Echart, Ichart, f, g) := openSchubertChartDataGLn(2, toList lam, toList muInput);
    Rsch := ring f;

    -- Work over the product of the Schubert chart and the Nnabla coefficient space.
    Rtotal := Rsch ** R;
    fR := sub(f, Rtotal);
    gR := sub(g, Rtotal);
    NnablaR := sub(Nnabla, Rtotal);
    KRgens := sub(gens K, Rtotal);

    -- The composite K -> F[u]^2/u^a -> coker(fR), written using gR.
    compositeMap := gR * KRgens;
    rankAtMostOneIdeal := minors(2, compositeMap);

    (Rsch, Rtotal, fR, gR, NnablaR, KRgens, compositeMap, rankAtMostOneIdeal)
);

(Rsch, Rtotal, fR, gR, NnablaR, KRgens, compositeMap, rankAtMostOneIdeal) = rankAtMostOneLocus(mu);
