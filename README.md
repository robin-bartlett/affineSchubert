# Affine Schubert and Affine Springer Calculations

This repository contains Macaulay2 scripts for calculations in the affine
Grassmannian for `GL_n`.

- `affGrSchubert.m2` computes Schubert varieties and affine Schubert charts.
- `affSpringer.m2` computes affine Springer loci inside those Schubert charts.

The base field in the scripts is `QQ`.

## `affGrSchubert.m2`

Load the file in Macaulay2 with:

```m2
load "affGrSchubert.m2"
```

The file provides two main entry points.

### `schubertGLn(n, lam)`

`schubertGLn(n, lam)` returns the homogeneous ideal of the Schubert variety
`Gr_{\le lambda}` in Plucker coordinates.

Input:

- `n`: the rank for `GL_n`.
- `lam`: a decreasing list of integers of length  `n`, interpreted as a dominant coweight of `GL_n`.

The function first normalizes `lam` by subtracting its last entry. It then sets
`N = lambda_0`, embeds the calculation in the finite quotient
`F[u]^n / u^N F[u]^n`, and writes the Schubert variety as a closed subvariety
of the ordinary Grassmannian `Gr(d, nN)`, where:

```m2
D = sum lambda
d = n*N - D
```

The returned ideal lives in a polynomial ring with Plucker coordinates
`p_I`, one for each `D`-element subset `I` of `{0, ..., nN-1}`, together with a
homogenizing variable `h`.

### `openSchubertGLn(n, lam, mu)`

`openSchubertGLn(n, lam, mu)` returns the ideal cutting out the affine open
chart of `Gr_{\le lambda}` around the torus-fixed point `u^mu`.

Input:

- `n`: the rank for `GL_n`.
- `lam`: a decreasing list of integers of length  `n`.
- `mu`: a decreasing list of integers of length  `n`.

The function requires `mu <= lambda` in dominance order.

The chart is built from the generic loop-group matrix

```text
E_ij = delta_ij u^(mu_i) + x_(i,j,0) + x_(i,j,1) u + ... + x_(i,j,mu_i-1) u^(mu_i-1).
```

The Schubert condition is imposed by exterior-power divisibility: for
`j = 1, ..., n`, every entry of `wedge^j E` must be divisible by

```text
u^(lambda_(n-j) + ... + lambda_(n-1)).
```

The companion function `openSchubertChartDataGLn(n, lam, mu)` returns more data:

```m2
(R, S, uVar, E, I, f)
```

where:

- `R` is the coordinate ring of the affine chart.
- `S = R[u]`.
- `uVar` is the variable `u` in `S`.
- `E` is the generic loop matrix.
- `I` is the Schubert chart ideal in `R`.
- `f` is the tautological frame over `R/I`, expressed in
  `F[u]^n / u^N F[u]^n`.

## Example: `n = 2`, `lam = (2,0)`, `mu = (1,1)`

Run the following in Macaulay2:

```m2
load "affGrSchubert.m2"

n = 2;
lam = {2, 0};
mu = {1, 1};

K = schubertGLn(n, lam);
(R, S, uVar, E, I, f) = openSchubertChartDataGLn(n, lam, mu);
```

For the global Schubert variety, `lambda` is already normalized, so:

```text
N = 2
nN = 4
D = 2
d = 2
```

Thus the global model is a closed subvariety of `Gr(2,4)`. The Plucker ring is:

```m2
QQ[p_(0,1), p_(0,2), p_(1,2), p_(0,3), p_(1,3), p_(2,3), h]
```

The ideal returned by `schubertGLn(2, {2,0})` is:

```m2
ideal (
    0,
    p_(1,3),
    p_(1,2) + p_(0,3),
    p_(0,1) - h,
    p_(0,3)^2 - p_(2,3)*h
)
```

For the open chart around `mu = (1,1)`, the coordinate ring is:

```m2
QQ[x_(0,0,0), x_(0,1,0), x_(1,0,0), x_(1,1,0)]
```

The generic loop matrix is:

```m2
| u + x_(0,0,0)   x_(0,1,0)     |
| x_(1,0,0)       u + x_(1,1,0) |
```

The Schubert condition says that `det(E)` is divisible by `u^2`. Since

```text
det(E) = u^2 + (x_(0,0,0) + x_(1,1,0))u
         + x_(0,0,0)x_(1,1,0) - x_(0,1,0)x_(1,0,0),
```

the coefficients of `u^0` and `u^1` must vanish. Macaulay2 returns:

```m2
ideal (
    -x_(0,1,0)*x_(1,0,0) + x_(0,0,0)*x_(1,1,0),
     x_(0,0,0) + x_(1,1,0)
)
```

The tautological frame `f` has columns indexed by `u^k E e_j` for
`0 <= k < N - mu_j`. Here each `N - mu_j` is `1`, so there are two columns,
corresponding to `E e_0` and `E e_1`. Rows are ordered as
`e_0, u e_0, e_1, u e_1`. Over the quotient by `I`, Macaulay2 prints:

```m2
| -x_(1,1,0)  x_(0,1,0) |
|  1          0         |
|  x_(1,0,0)  x_(1,1,0) |
|  0          1         |
```

The entry `-x_(1,1,0)` appears because the chart ideal contains
`x_(0,0,0) + x_(1,1,0)`.
