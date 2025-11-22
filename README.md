# veigs

Verified eigenvalue solver for the generalized eigenproblem $A x = \lambda B x$ (MATLAB).

This toolbox provides verified (interval) bounds for eigenvalues of the symmetric generalized eigenproblem using LDL-based separation and the Lehmann–Behnke inclusion. It is intended for small-to-medium dense or sparse problems.

Authors:
- Xuefeng Liu
- Yuuka Yanagisawa


Key files:
- [`veigs`](veigs.m) — main high-level routine for cluster-aware verified eigenvalue bounds. See [veigs.m](veigs.m).
- [`veig`](veig.m)  — interval solver for small generalized eigenproblems used in Lehmann–Behnke steps. See [veig.m](veig.m).
- [`GetInertia`](GetInertia.m) — inertia counter for 1×1 and 2×2 diagonal blocks from LDL. See [GetInertia.m](GetInertia.m).
- License: [LICENSE](LICENSE)

Requirements
- INTLAB (interval arithmetic toolbox for MATLAB/Octave): `intval`, `mid`, `inf`, `sup`, `hull`, `isspd`, etc.

Quick usage
- Call the main function as:
```matlab
[lambdas, ind_range] = veigs(A, B, NumOfEigs, SIGMA)
```

% Exmaple 1:
n = 8;
A = eye(n);
B = infsup(hilb(n)-1E-13, hilb(n)+1E-13);

% Find eigenvalue(s) near 0.1
[bounds, ind] = veigs(A, B, 0.1);

% Smallest real
[bounds, ind] = veigs(A, B, 1, 'sa');

% Largest real
[bounds, ind] = veigs(A, B, 1, 'la');


Example 2 — banded test
```matlab
% Example 2:
n = 10;
A = diag(ones(n,1))*6; A(1,1)=5; A(n,n)=5;
A = A + diag( ones(n-1,1),1)*(-4) + diag( ones(n-1,1),-1)*(-4);
A = A + diag( ones(n-2,1),2) + diag( ones(n-2,1),-2);
B = hilb(n) * 232792560;

% Largest eigenvalue verified bounds
[bound, ind] = veigs(A, B, 1,'la');
```

Notes
- The algorithm builds on LDL factorizations of A - λB (midpoint) and uses verified interval arithmetic to certify bounds. For clustered eigenvalues it applies Lehmann–Behnke complementary variational principles.
- Designed for small-to-medium problems; performance depends on problem size and whether sparse routines (`eigs`) are applicable.
- The code uses INTLAB interval objects and utilities (e.g. intval, mid, inf, sup, hull, isspd).

References
- H. Behnke, "The calculation of Guaranteed bounds for eigenvalues using complementary variational principles", Computing 47, 11–27 (1991).
- N. Yamamoto, "A simple method for error bounds of eigenvalues of symmetric matrices", Linear Algebra and its Applications, 234 (2001), 227–234.
