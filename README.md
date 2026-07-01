# Affine Schubert and Affine Springer Calculations

This repository contains Macaulay2 scripts for calculations in the affine
Grassmannian for `GL_n`.

- `affGrSchubert.m2` computes Schubert varieties in the affine Grassmannian.
- `affSpringer.m2` computes families which degenerate to affine Springer loci inside those Schubert charts.

The base field `F` can be changed in `affGrSchubert.m2`.

## `affGrSchubert.m2`

Load the file in Macaulay2 with:

```m2
load "affGrSchubert.m2"
```

The file has two main functions.

### `schubertGLn(n, lam)`

This returns the homogeneous ideal of the Schubert variety
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

This returns the ideal cutting out the affine open
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
(R, S, uVar, E, I, f, g)
```

where:

- `R` is the coordinate ring of the affine chart.
- `S = R[u]`.
- `uVar` is the variable `u` in `S`.
- `E` is the generic loop matrix.
- `I` is the Schubert chart ideal in `R`.
- `f` describes the map of vector bundles `E_univ -> F[u]^n/u^NF[u]^n \otimes R/I` where `E_univ` denotes the universal lattice on the open chart `R/I`.
- `g` is the quotient map `F[u]^n / u^N F[u]^n \otimes R/I -> coker(f)` written in the explicit quotient basis `u^r e_i` for `0 <= r < mu_i`.

## Example: `n = 2`, `lam = (2,0)`, `mu = (1,1)`

Run the following in Macaulay2:

```m2
load "affGrSchubert.m2"

n = 2;
lam = {2, 0};
mu = {1, 1};

K = schubertGLn(n, lam);
(R, S, uVar, E, I, f, g) = openSchubertChartDataGLn(n, lam, mu);
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

## `affSpringer.m2`

Load the affine Springer file with:

```m2
load "affSpringer.m2"
```

This file loads `affGrSchubert.m2` automatically. Its main entry point is:

```m2
openNablaSpringer(n, lam, mu, e)
```

It returns a pair:

```m2
(J, f)
```

where:

- `J` is the ideal cutting out the closed locus where a connection preserves
  the tautological lattice on the affine Schubert chart around `u^mu`.
- `f` is the same tautological frame as in `openSchubertChartDataGLn`, but
  base-changed to the coordinate ring used for the affine Springer equations.

The input validation is the same as for `openSchubertGLn(n, lam, mu)`, with one
additional parameter:

- `e`: a positive integer.

The connection has the form

```text
nabla = N(u) + (u - pi_1)...(u - pi_(e-1)) u d/du,
```

where

```text
N(u) = Y_1 u + Y_2 u^2 + ... + Y_(N-1) u^(N-1),
N = lam_0,
Y_k = (y_(i,j,k))_(i,j).
```

The coordinate ring for `J` is obtained from the open Schubert chart ring
`R/I` by adjoining the connection coefficients `y_(i,j,k)`. If `e > 1`, it
also adjoins the parameters `pi_1, ..., pi_(e-1)`.

Internally, the script writes the generic chart matrix as `X`, computes its
classical adjugate `X^*`, and imposes preservation of the lattice by requiring
every entry of

```text
X^* nabla(X)
```

to have `u`-adic valuation at least `sum lam`. Equivalently, the coefficients
of `u^0, ..., u^(sum lam - 1)` in all entries are added to the ideal `J`.

## Example: affine Springer equations for `n = 2`

Continuing the previous example with `lam = {2,0}` and `mu = {1,1}`, take
`e = 1`:

```m2
load "affSpringer.m2"

n = 2;
lam = {2, 0};
mu = {1, 1};
e = 1;

(J, fA) = openNablaSpringer(n, lam, mu, e);
```

Here `N = lam_0 = 2`, so

```text
N(u) = Y_1 u.
```

The affine Springer coordinate ring is the Schubert chart quotient, with four
new connection variables adjoined:

```m2
(QQ[x_(0,0,0), x_(0,1,0), x_(1,0,0), x_(1,1,0)]
 / ideal(
     -x_(0,1,0)*x_(1,0,0) + x_(0,0,0)*x_(1,1,0),
      x_(0,0,0) + x_(1,1,0)
   )
)[y_(0,0,1), y_(0,1,1), y_(1,0,1), y_(1,1,1)]
```

For readability, write:

```text
b = x_(0,1,0),  c = x_(1,0,0),  d = x_(1,1,0),
A = y_(0,0,1),  B = y_(0,1,1),  C = y_(1,0,1),  D = y_(1,1,1).
```

Since the Schubert chart relation gives `x_(0,0,0) = -d`, Macaulay2 prints the
tautological frame as:

```m2
| -d  b |
|  1  0 |
|  c  d |
|  0  1 |
```

With these abbreviations, the ideal `J` returned by
`openNablaSpringer(2, {2,0}, {1,1}, 1)` is generated by:

```text
-d^2 A + c*d B + b*d C + d^2 D + d^2
 b*d A + d^2 B - b C - b*d D - b
 c*d A - c B + d C - c*d D - c
 d A - c*d B - b*d C - d*D - d
```

The four generators are the coefficients of powers below `u^2` in the entries
of `X^* nabla(X)`. The returned frame `fA` is the base change of the Schubert
chart frame to the affine Springer coordinate ring.
