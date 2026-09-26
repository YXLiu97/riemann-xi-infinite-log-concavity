# Reproducibility supplement

This supplement contains the source code, exact computational artifacts, and reproducibility information for the computer-assisted parts of the manuscript. These materials are provided to support verification and reproduction of the certified computations used in the paper.

## Reproduce the full certificate

Use Python 3.12 with `mpmath==1.3.0` and `gmpy2==2.3.1`. The recorded run used Python 3.12.14 and the gmpy integer backend. Install these packages in a suitable environment and run from this supplement directory:

```text
python scripts/xi_finite_certificate.py --nmax 12016 --q 256 --dps 130 --workers 8 --steps 14 --certify-to 12000 --tag reproduced
```

The stored evidence uses the tag `current_review`. The `reproduced` tag in the command above writes new files without replacing it. The tag affects file naming, not the mathematical calculation. Keep Python assertions enabled: do not use `-O`, `-OO`, or `PYTHONOPTIMIZE`. Do not add `--reuse-nodes` for a full reproduction; that option trusts saved node data instead of rechecking the endpoint signs. The local `mcp` module only implements the original program's progress logging calls.

The unchanged program's SHA-256 is

```text
01816166fa0b317f0b116ca3ac45ffa198ec9348aa3c73d48f4d8a571e0211fe
```

The full run freshly generated nodes and 40 quadrature segments, enclosed all 12017 initial coefficients, and checked 14 nonlinear iterations. It passed in 121.112 seconds in the review environment. Worker scheduling can change the exact endpoints of enclosing sums slightly; bit-for-bit identity with a previous parallel run is not required for a valid enclosure. The inclusion assertions and the strict final lower bound are the relevant checks.

## Exact certificate files

All six files are under `downloads/riemann_xi_infinite_log_concavity/` with prefix `finite_certificate_current_review_`:

| Suffix | Content |
|---|---|
| `nodes.json` | All 256 disjoint root brackets, certified opposite endpoint signs, derivative and weight intervals |
| `raw_moments.json` | 12017 enclosing raw Gaussian sums |
| `errors.json` | Cauchy error bounds, real-tail bounds, certified phase-root brackets, constant audit |
| `lambda.json` | 12017 positive starting coefficient enclosures |
| `ratios.json` | 12000 final ratio intervals |
| `report.json` | Parameters, software assumptions, source hash, every retained index range, and pass results |

An endpoint `[sign, mantissa, exponent, bitcount]` represents the exact dyadic rational `(-1)^sign * mantissa * 2^exponent`. Parse mantissas as arbitrary-precision integers. Displayed decimal endpoints are rounded outward; they are not standard deviations or errors estimated from repeated floating-point runs.

The final report must contain `status=passed`, `certified_rational_lower_bound="10"`, positive entries through index 12002 at iteration 14, and 12000 checked ratios. The smallest computed final lower endpoint occurs at index 12000. Its outward decimal interval is `[10.9042985947999313085054, 10.9042985947999313085055]`.

## Analytic scalar checks

```text
python scripts/check_analytic_constants.py --output analytic-constants-reproduced.json
```

This uses 100-digit mpmath interval arithmetic. The stored `analytic-constants.json` contains 34 passing checks, their exact formulas, dyadic endpoints, and positive margins. In particular, the shrinking-disk constant is less than 1500 and the gamma remainder's reduced bound is less than 4. These checks do not establish the analytic reductions to the listed scalar endpoints, nor do they verify contours, holomorphic branches, or induction hypotheses.

## Evidence and provenance scope

| Record | Meaning |
|---|---|
| `review-runtime.json` | Fresh full-run Python and package versions and original source hash |
| `computation-record.json` | Current full-run execution and verification scope, with artifact hashes |
| `triangle-replay.json` | Earlier separate 512-bit integer-dyadic replay of the archived initial intervals; conditional on those enclosures |
| `bundle-audit.json` | Earlier manifest and source consistency checks of the unchanged original bundle; its paths refer to that original bundle |
| `retrieval.json` | Public references consulted, their URLs and hashes; retrieval paths identify the review workspace cache |
| `historical-provenance/` | Supplied blueprints and historical automated verification report, with their original claims preserved as historical claims |

The fresh full run goes beyond the earlier conditional triangle replay by regenerating the moment enclosures. It uses the supplied implementation, not an independently written moment integrator. The integer replay uses separate recurrence arithmetic but accepts archived starting enclosures. Neither is an independent review of the whole paper.

The arithmetic trust boundary is mpmath 1.3.0's supported interval functions and explicit directed basic operations, with the gmpy2 2.3.1 integer backend. mpmath documents its interval support as experimental. This supplement does not certify the numerical library itself. The main theorem also depends on the mathematical arguments in Sections 2--4 and 6, which do not run inside the finite certificate script.

