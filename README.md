# Reproducibility supplement

This supplement contains the source code, exact computational artifacts, and reproducibility information for the computer-assisted parts of the manuscript. These materials are provided to support verification and reproduction of the certified computations used in the paper.

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

The script `xi_finite_certificate.py` produces the rigorous finite certificate used in the proof of the finite escape result. The script `check_analytic_constants.py` independently checks the explicit scalar inequalities appearing in the analytic estimates.

The local `mcp` module only implements progress logging calls used by the certificate script. It is not a mathematical input to the proof.

## Software requirements

The computations were carried out with Python 3.12 and the following packages:

```text
mpmath==1.3.0
gmpy2==2.3.1
```

The recorded run used Python 3.12.14.

The interval computations rely on `mpmath` interval arithmetic together with explicitly directed basic operations where indicated in the source code.

## Reproduce the finite certificate

Run the following command from the `supplement` directory:

```text
python scripts/xi_finite_certificate.py --nmax 12016 --q 256 --dps 130 --workers 8 --steps 14 --certify-to 12000 --tag reproduced
```

Python assertions must be enabled. In particular, do not use `-O`, `-OO`, or `PYTHONOPTIMIZE`.

For a full reproduction, do not add the optional `--reuse-nodes` flag. Without this flag, the program regenerates the Gauss--Legendre node brackets and certifies the endpoint signs before carrying out the remaining interval computations.

The tag `reproduced` affects only the names of the newly generated output files. It does not change the mathematical calculation. The archived certificate supplied with the repository uses the tag `current_review`.

A full run:

- certifies all 256 Gauss--Legendre node brackets and their positive weights;
- computes enclosing raw Gaussian moment intervals;
- encloses the coefficients \(\lambda_0,\ldots,\lambda_{12016}\);
- verifies positivity through 14 nonlinear iterations; and
- checks all 12,000 final ratio inequalities required by the finite certificate.

Worker scheduling may change the exact endpoints of some enclosing sums slightly because the order of interval additions may differ. Bit-for-bit identity with the archived parallel run is therefore not required. The relevant conditions are interval inclusion, positivity of all certified entries, and the strict final ratio bound.

## Archived finite certificate

The six archived certificate files are stored in

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
| `finite_certificate_current_review_nodes.json` | Certified brackets for all 256 Gauss--Legendre nodes, endpoint signs, derivative intervals, and weight intervals |
| `finite_certificate_current_review_raw_moments.json` | Enclosing raw Gaussian moment sums |
| `finite_certificate_current_review_errors.json` | Certified quadrature-error bounds, tail bounds, and phase-root brackets |
| `finite_certificate_current_review_lambda.json` | Positive interval enclosures for \(\lambda_0,\ldots,\lambda_{12016}\) |
| `finite_certificate_current_review_ratios.json` | Certified final ratio intervals |
| `finite_certificate_current_review_report.json` | Parameters, software information, retained index ranges, and the final pass summary |

An interval endpoint stored as

```text
[sign, mantissa, exponent, bitcount]
```

represents the exact dyadic rational

```text
(-1)^sign * mantissa * 2^exponent.
```

The mantissas in the JSON files are arbitrary-precision integers. Thus the stored endpoints are the actual binary interval endpoints used by the calculation, rather than decimal approximations reconstructed from printed output.

Displayed decimal endpoints are rounded outward.

For the archived certificate, the final ratio test verifies

```text
r_{14,n} > 10
```

for every

```text
1 <= n <= 12000.
```

The smallest computed lower endpoint at iteration 14 occurs at `n = 12000`. Its outward decimal enclosure is

```text
[10.9042985947999313085054,
 10.9042985947999313085055].
```

The assertion `> 10` is made using the underlying exact binary interval endpoints, not the displayed decimal digits.

## Reproduce the analytic scalar checks

The explicit scalar inequalities used in the analytic estimates can be checked separately by running

```text
python scripts/check_analytic_constants.py --output analytic-constants-reproduced.json
```

from the `supplement` directory.

The archived output is

```text
analytic-constants.json
```

and contains 34 scalar checks performed with 100-digit `mpmath` interval arithmetic.

These include, among others, the numerical bounds used in the saddle-point estimates, the Gaussian error estimates, the shrinking-disk argument, and the finite-kernel estimates.

The scalar-check script verifies only the stated numerical inequalities. It does not verify the analytic reductions leading to those inequalities, the contour deformations, holomorphic continuation, monotonicity arguments, or the induction arguments in the manuscript. Those parts are proved mathematically in the paper.

## Relation to the proof

The finite certificate produced by `xi_finite_certificate.py` is used for the bounded range of indices in the finite escape argument. The analytic estimates in the manuscript treat the complementary unbounded range.

The computational scripts therefore certify only the finite numerical statements explicitly identified in the paper. The deduction of infinite log-concavity from these finite and analytic inputs is a separate mathematical argument given in the manuscript.

## Arithmetic trust boundary

The certified computations rely on the correctness of the supported interval operations in `mpmath` 1.3.0 and on the explicitly directed elementary operations used by the certificate program. The computations may use the `gmpy2` 2.3.1 integer backend when available.

The supplement does not constitute a formal verification of the numerical library itself.
