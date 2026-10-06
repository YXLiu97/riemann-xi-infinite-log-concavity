# Reproducibility supplement

This supplement contains the source code, exact computational data, and reproducibility information for the computer-assisted parts of the manuscript. Its purpose is to support inspection and reproduction of the finite interval computations used in the proof.

## Contents

The directory is organized as follows:

```text
supplement/
├── README.md
├── analytic-constants.json
├── scripts/
│   ├── check_analytic_constants.py
│   └── xi_finite_certificate.py
├── mcp/
│   ├── __init__.py
│   └── server.py
└── downloads/
└── riemann_xi_infinite_log_concavity/
├── finite_certificate_current_review_nodes.json
├── finite_certificate_current_review_raw_moments.json
├── finite_certificate_current_review_errors.json
├── finite_certificate_current_review_lambda.json
├── finite_certificate_current_review_ratios.json
└── finite_certificate_current_review_report.json
```

The script `scripts/xi_finite_certificate.py` performs the finite interval verification used in the manuscript. It certifies the Gauss--Legendre nodes and weights, evaluates the sums $Q_n$, constructs interval enclosures for $\lambda_0,\ldots,\lambda_{12016}$, verifies positivity through 14 iterations, and checks the final ratio inequalities.

The script `scripts/check_analytic_constants.py` separately verifies the explicit scalar inequalities appearing in the analytic estimates.

The local `mcp` module is used only for progress logging and is not a mathematical input to the verification.

## Software requirements

The computations were carried out with Python 3.12 and the following packages:

```text
mpmath==1.3.0
gmpy2==2.3.1
```

The recorded computation used Python 3.12.14.

The interval calculations use `mpmath` interval arithmetic together with explicitly directed elementary operations where indicated in the source code.

## Reproduce the finite verification

From the `supplement` directory, run

```text
python scripts/xi_finite_certificate.py --nmax 12016 --q 256 --dps 130 --workers 8 --steps 14 --certify-to 12000 --tag reproduced
```

Python assertions must be enabled. In particular, do not use `-O`, `-OO`, or set `PYTHONOPTIMIZE`.

For a complete reproduction, do not use the optional `--reuse-nodes` flag. Without this flag, the program regenerates the 256 Gauss--Legendre root brackets and verifies the endpoint signs before carrying out the remaining interval computations.

The value of `--tag` affects only the names of the generated output files. The archived files supplied with this repository use the tag `current_review`; using `--tag reproduced` avoids replacing those files and does not alter the mathematical computation.

A complete run verifies:

- 256 pairwise disjoint brackets containing the zeros of $P_{256}$;
- positive interval enclosures for the corresponding Gauss--Legendre weights;
- outward interval enclosures of the sums $Q_n$ for $0\le n\le12016$;
- the bounds $E_n$ and $T_n$ used in the coefficient enclosures;
- positive intervals containing $\lambda_0,\ldots,\lambda_{12016}$;
- positivity of the iterated array through step 14; and
- all 12,000 final ratio inequalities required in the manuscript.

Parallel worker scheduling can slightly change the binary endpoints of some enclosing sums because interval additions may occur in a different order. Bit-for-bit agreement with the archived parallel run is therefore not required. The relevant mathematical conditions are interval inclusion, positivity of the certified entries, and the strict final ratio bounds.

## Archived finite verification

The six archived JSON files are stored in

```text
downloads/riemann_xi_infinite_log_concavity/
```

and have the prefix

```text
finite_certificate_current_review_
```

Their contents are as follows:

| File | Content |
|---|---|
| `finite_certificate_current_review_nodes.json` | Certified brackets for all 256 zeros of $P_{256}$, endpoint signs, derivative intervals, and positive Gauss--Legendre weight intervals |
| `finite_certificate_current_review_raw_moments.json` | Outward interval enclosures of the sums $Q_n$, $0\le n\le12016$ |
| `finite_certificate_current_review_errors.json` | Certified intervals $[l_n,u_n]$ and outward bounds for $M_n$ |
| `finite_certificate_current_review_lambda.json` | Positive interval enclosures for $\lambda_0,\ldots,\lambda_{12016}$ |
| `finite_certificate_current_review_ratios.json` | Certified ratio intervals at the final iteration |
| `finite_certificate_current_review_report.json` | Parameters, software information, retained index ranges, iteration summaries, and verification results |

An interval endpoint stored as

```text
[sign, mantissa, exponent, bitcount]
```

represents the exact binary value

```text
(-1)^sign * mantissa * 2^exponent.
```

The mantissas are arbitrary-precision integers. Hence the stored values are the actual binary interval endpoints produced by the computation, rather than decimal approximations reconstructed from printed output.

Displayed decimal endpoints are obtained from these exact binary values by outward rounding.

For the archived computation, the final verification establishes

```text
r_{14,n} > 10
```

for every

```text
1 <= n <= 12000.
```

The smallest certified lower endpoint at iteration 14 occurs at `n = 12000`. Its outward decimal enclosure is

```text
[10.9042985947999313085054,
10.9042985947999313085055].
```

The strict inequality `> 10` is checked using the underlying exact binary interval endpoints, not the displayed decimal digits.

## Source identification

The main program used for the finite verification is

```text
scripts/xi_finite_certificate.py
```

with SHA-256

```text
01816166fa0b317f0b116ca3ac45ffa198ec9348aa3c73d48f4d8a571e0211fe
```

This hash identifies the exact version of the program corresponding to the finite computation described in the manuscript.

## Reproduce the analytic scalar checks

The explicit scalar inequalities used in the analytic estimates can be verified separately by running

```text
python scripts/check_analytic_constants.py --output analytic-constants-reproduced.json
```

from the `supplement` directory.

The archived output is

```text
analytic-constants.json
```

and contains 34 scalar checks carried out with 100-digit `mpmath` interval arithmetic.

These checks include the numerical constants appearing in the analytic estimates, including

```text
12 / (1 - 2^(-1/3))^3 < 1500.
```

The script verifies only the stated scalar inequalities. It does not verify the analytic arguments that reduce the proof to those inequalities, such as holomorphic continuation, contour estimates, monotonicity arguments, or induction. Those parts are proved mathematically in the manuscript.

## Relation to the proof

The script `xi_finite_certificate.py` verifies the finite range of indices required in the manuscript. The complementary unbounded range is treated by the analytic lemmas.

Accordingly, the computational scripts establish only the finite numerical statements explicitly used in the proof. The deduction of positivity for all iterations and all indices from the finite and analytic estimates is a separate mathematical argument in the manuscript.

## Interval arithmetic and verification scope

The finite computation relies on the inclusion properties of the interval operations provided by `mpmath` 1.3.0 and on the explicitly directed elementary operations used in `xi_finite_certificate.py`. The recorded computation used `gmpy2` 2.3.1 as the integer backend.

Approximate `mpmath.mp` calculations are used only to select candidate Gauss--Legendre root brackets and candidate intervals $[l_n,u_n]$. These intervals are subsequently certified by outward interval inequalities. Thus no unchecked approximate value is used as an endpoint of a certified enclosure.

The `mpmath` documentation describes the inclusion properties of its interval arithmetic but also marks the interval module as experimental. Consequently, this supplement does not constitute a formal verification of the numerical library itself.
