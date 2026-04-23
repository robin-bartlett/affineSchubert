-- affGrSchubert.m2
-- Schubert varieties in the affine Grassmannian for GL_n.
--
-- FUNCTIONS
--
-- (1) schubertGLn(n, lam)
--     Returns K, the ideal of the Schubert variety Gr_lam in Plucker
--     coordinates. K is an ideal in the polynomial ring Grass = F[p_I, h].
--     The ring Grass is recoverable as ring(K).
--
-- (2) openSchubertGLn(n, lam, mu)
--     Returns I, the ideal defining the open cell of Gr_lam around u^mu.
--     The companion function openSchubertChartDataGLn additionally returns
--     the chart data and the tautological frame f described below.
--     Uses exterior-power divisibility on a generic loop-group matrix.
--

F = QQ;

-------------------------------------------------------------------------------
-- AUXILIARY FUNCTIONS
-------------------------------------------------------------------------------

-- normalizeLambda(lam)
-- Input:   lam - coweight list
-- Output:  lambda = lam shifted so that lambda_{n-1} = 0
normalizeLambda = lam -> (
    minLam := lam#(#lam-1);
    apply(lam, l -> l - minLam)
);

-------------------------------------------------------------------------------
-- MAIN FUNCTION 1: schubertGLn
-- Global coordinate ring of the Schubert variety Gr_{\leq \lambda}
-- in Plucker coordinates.
-------------------------------------------------------------------------------

schubertGLn = (n, lam) -> (
    -- Validate input
    if #lam != n then error("schubertGLn: lambda must have length n");
    if not all(0..n-2, i -> lam#i >= lam#(i+1)) then
        error("schubertGLn: lambda must be weakly decreasing");

    -- Normalize so that lambda_{n-1} = 0
    lambda := normalizeLambda(lam);

    -- Derived parameters
    N := lambda#0;              -- loop length
    nN := n * N;                -- dimension of V
    D := sum lambda;            -- |lambda| = size of Plucker subsets
    d := nN - D;                -- dim of subspace in Gr(d, nN)

    if d <= 0 or d >= nN then
        error("schubertGLn: degenerate case, d = " | toString d);

    -- Build variable list: x_(i,j,c) for i < j, c = 0..lambda_i - lambda_j - 1 
    varList := flatten apply(n, i ->
                   flatten apply(toList(i+1..n-1), j ->
                       if lambda#i > lambda#j then
                           apply(lambda#i - lambda#j, c -> x_(i,j,c))
                       else {}));

    R := if #varList == 0 then F[] else F[varList];

    -- Build generic matrix Q of size nN x d
    -- Columns indexed by (j, k) for j = 0..n-1, k = 0..N-lambda_j-1
    -- Row s = i*N + r corresponds to u^r e_i
    colList := flatten apply(n, j ->
                   apply(N - lambda#j, k -> (j, k)));

    Q := matrix apply(nN, s -> (
        i := s // N;
        r := s % N;
        apply(colList, col -> (
            j := col#0;
            k := col#1;
            if i == j then (
                if r == k + lambda#j then 1_R else 0_R
            ) else if i < j and lambda#i > lambda#j then (
                c := r - k - lambda#j;
                if c >= 0 and c < lambda#i - lambda#j then
                    x_(i, j, c)
                else 0_R
            ) else 0_R
        ))
    ));

    -- Plucker ring
    Slist := subsets(nN, D);
    var := apply(Slist, I -> p_(toSequence I));
    Grass := F[var, symbol h,
                Degrees => toList apply(var, i -> {1}) | {{1}}];
    hGrass := h_Grass;

    -- Warn for large Plucker spaces
    pluckerCount := binomial(nN, D);
    if pluckerCount > 500 then
        << "Warning: " << pluckerCount << " Plucker coordinates. This may be slow." << endl;

    -- Compute K: kernel of Plucker embedding map f: R -> Grass
    -- p_I maps to det of complement-of-I rows of Q
    v := matrix({apply(Slist, I ->
                    det submatrix(Q,
                        select(toList(0..nN-1), j -> not member(j, set I)),
                        toList(0..d-1)))
                | {1_R}});
    f := map(R, Grass, v);
    K := homogenize(ker f, hGrass);

    K
);

-------------------------------------------------------------------------------
-- MAIN FUNCTION 2: openSchubertGLn
-- Affine open in Gr_{\leq \lambda} around the torus-fixed point u^mu
-- (via the u^mu-translate of the big cell).
--
-- Builds the generic loop-group matrix
--     E_{ij} = delta_{ij} u^{mu_i} + x_(i,j,0) + x_(i,j,1) u + ... + x_(i,j,mu_i-1) u^{mu_i-1}
-- where x_(i,j,r) is the coordinate variable for the coefficient of u^r
-- in the (i,j)-entry of E. The coordinate ring is R = F[x_(i,j,r)].
--
-- Imposes the Schubert condition: /\^j E divisible by u^{lam_{n-j+1}+...+lam_n}
-- for j = 1,...,n. 
--
-- More, specifically:
--
-- openSchubertChartDataGLn(n, lam, mu) returns (R, S, uVar, E, I, f) where:
--   R    = F[x_(i,j,r)], the coordinate ring of the chart,
--   S    = R[u],
--   uVar = the polynomial variable u in S,
--   E    = the generic n x n loop-group matrix over S defined above,
--   I    = the ideal in R cutting out the Schubert chart by the exterior-power
--          divisibility conditions,
--   f    = the nN x d matrix over Q = R/I whose columns are the images of
--          u^k E e_j in F[u]^n / u^N F[u]^n for 0 <= k < N - mu_j, where
--          N = lam_0 and d = sum_j (N - mu_j). This gives a description of
--          the restriction to Spec(Q) \subet Gr_{\leq \lam} of the 
--          map S -> F[u]^n / u^N F[u]^N \otimes_F O_{Gr_\leq lam} of vector 
--          vector bundles, where S is the tautological bundle on Gr_{\leq \lam}.
--
-- openSchubertGLn(n lam, mu) returns the ideal I in 
-- openSchubertChartDataGLn(n, lam, mu).
-------------------------------------------------------------------------------

openSchubertChartDataGLn = (n, lam, mu) -> (
    -- Validate input
    if #lam != n then error("openSchubertGLn: lambda must have length n");
    if #mu != n then error("openSchubertGLn: mu must have length n");
    if not all(0..n-2, i -> lam#i >= lam#(i+1)) then
        error("openSchubertGLn: lambda must be weakly decreasing");
    if not all(0..n-2, i -> mu#i >= mu#(i+1)) then
        error("openSchubertGLn: mu must be weakly decreasing");
    if mu#(n-1) < 0 then
        error("openSchubertGLn: mu entries must be nonnegative");
    if sum mu != sum lam then
        error("openSchubertGLn: |mu| must equal |lambda|");
    scan(n-1, k -> (
        if sum take(mu, k+1) > sum take(lam, k+1) then
            error("openSchubertGLn: mu is not <= lambda in dominance order");
    ));

    -- Build variable list: x_(i,j,r) for i = 0..n-1, j = 0..n-1, r = 0..mu_i-1
    varList := flatten apply(n, i ->
                   flatten apply(n, j ->
                       apply(mu#i, r -> x_(i,j,r))));

    -- Build R = F[x_(i,j,r)] and S = R[u]
    R := if #varList == 0 then F[] else F[varList];
    S := R[getSymbol "u"];
    uVar := S_0;

    -- Build E: E_{ij} = delta_{ij} u^{mu_i} + sum_{r=0}^{mu_i-1} x_(i,j,r) u^r
    E := matrix table(n, n, (i,j) -> (
        d := mu#i;
        diag := if i == j then uVar^d else 0_S;
        if d <= 0 then diag
        else diag + sum(d, r -> sub(x_(i,j,r), S) * uVar^r)
    ));

    -- Impose Schubert conditions: for j = 1,...,n, the j-th exterior power
    -- of E must be divisible by u^{lambda_{n-j+1} + ... + lambda_n}
    rels := {};
    for j from 1 to n do (
        expo := sum for k from n-j to n-1 list lam#k;
        wedgeE := exteriorPower(j, E);
        scan(flatten entries wedgeE, f -> (
            for k from 0 to expo - 1 do (
                c := coefficient(uVar^k, f);
                if c != 0 then rels = append(rels, sub(c, R));
            );
        ));
    );

    I := if #rels == 0 then ideal(0_R) else ideal rels;

    Q := R / I;
    N := lam#0;
    colList := flatten apply(n, j ->
                   apply(N - mu#j, k -> (j, k)));
    -- Build the frame of L/u^N inside F[u]^n/u^N.
    -- Row (i*N + r) corresponds to u^r e_i; column (j,k) to u^k E e_j.
    -- Entry (i*N + r, (j,k)) = coefficient of u^r in u^k E_{ij}, taken in Q.
    f := matrix apply(n * N, row -> (
        i := row // N;
        r := row % N;
        apply(colList, col -> (
            j := col#0;
            k := col#1;
            if r < k then 0_Q else sub(coefficient(uVar^(r-k), E_(i,j)), Q)
        ))
    ));

    (R, S, uVar, E, I, f)
);

openSchubertGLn = (n, lam, mu) -> (
    (R, S, uVar, E, I, f) := openSchubertChartDataGLn(n, lam, mu);
    I
);
