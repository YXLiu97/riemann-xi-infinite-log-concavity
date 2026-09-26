"""Directed scalar checks for the xi manuscript; not a whole-proof verifier.

Run with mpmath 1.3.0. Exact dyadic endpoints are retained in the JSON output.
The analytic reductions to the endpoints below must be checked in the paper.
"""
import argparse
import json
from pathlib import Path
import mpmath
from mpmath import iv
from mpmath.libmp import mpf_cmp, fzero


def raw(value):
    return [[int(x) for x in endpoint] for endpoint in value._mpi_]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("analytic-constants.json"))
    args = parser.parse_args()
    iv.dps = 100
    v = iv.mpf
    pi = iv.pi
    exp, sqrt, log = iv.exp, iv.sqrt, iv.log
    rows = []

    def less(name, formula, lhs, rhs):
        lhs, rhs = v(lhs), v(rhs)
        margin = rhs - lhs
        passed = mpf_cmp(margin._mpi_[0], fzero) > 0
        rows.append({"id": name, "inequality": formula, "lhs": raw(lhs),
                     "rhs": raw(rhs), "positive_margin": raw(margin),
                     "lhs_display": str(lhs), "passed": passed})

    a, t, n = v(6), v(630), v(12000)
    b = pi / 4
    d = sqrt(6 * pi / 7)
    less("S01", "(pi/2) exp(6) sqrt(36+(pi/4)^2) < 4000",
         pi / 2 * exp(a) * sqrt(a*a+b*b), 4000)
    less("S02", "sqrt(36+(pi/4)^2)/6 < 1.009", sqrt(a*a+b*b)/a, "1.009")
    less("S03", "633 < (pi/2) exp(6)", 633, pi/2*exp(a))
    less("S04", "1.4 < sqrt(2)/1.009", "1.4", sqrt(2)/v("1.009"))
    less("S05", "4a/(a-.1)^3+2exp(.1) < 2.4 at a=6",
         4*a/(a-v(".1"))**3 + 2*exp(v(".1")), "2.4")
    less("S06", "12a/(a-.1)^4+2exp(.1) < 2.4 at a=6",
         12*a/(a-v(".1"))**4 + 2*exp(v(".1")), "2.4")
    less("S07", "exp((9/4)^2/(.2*630)) < 1.2", exp((v(9)/4)**2/(v(".2")*t)), "1.2")
    less("S08", "1.64 < sqrt(6pi/7)", "1.64", d)
    c1 = v(".1")*3*sqrt(pi)/(4*v(".7")**(v(5)/2)*d)
    c2 = v(".6")*3*(v(81)/16)*sqrt(pi)/(2*v(".6")**(v(3)/2)*d)
    c3 = v(".6")*v(".48")*15*sqrt(pi)/(8*v(".6")**(v(7)/2)*d)
    c4 = v(".6")*v(".03")*105*sqrt(pi)/(16*v(".6")**(v(9)/2)*d)
    less("S09", "first three Gaussian constants sum < 14.4", c1+c2+c3, "14.4")
    less("S10", "fourth Gaussian constant < 1.3", c4, "1.3")
    less("S11", "2exp(.225)*630/(d*(.126*630-2.25)) < 12.6",
         2*exp(v(".225"))*t/(d*(v(".126")*t-v("2.25"))), "12.6")
    less("S12", "sqrt(630)*(12.6exp(-.0063*630)+8.8exp(-.007*630)) < 9",
         sqrt(t)*(v("12.6")*exp(-v(".0063")*t)+v("8.8")*exp(-v(".007")*t)), 9)
    less("S13", "sqrt(2)log(.52)+2/6 < -.58", sqrt(2)*log(v(".52"))+2/a, "-.58")
    less("S14", "2 < sqrt(2)log(6)-(2+pi/(2sqrt(2)))/6",
         2, sqrt(2)*log(a)-(2+pi/(2*sqrt(2)))/a)
    less("S15", "(4/3)^4 exp(-7pi/sqrt(2)) < 1e-6", (v(4)/3)**4*exp(-7*pi/sqrt(2)), "1e-6")
    less("S16", "81exp(-5pi/sqrt(2)) < .002", 81*exp(-5*pi/sqrt(2)), ".002")
    less("S17", "9exp(-5pi/sqrt(2)) < .001", 9*exp(-5*pi/sqrt(2)), ".001")
    less("S18", "64*20^2 exp(-6.6*20) < .01", 64*v(20)**2*exp(-v("6.6")*20), ".01")
    vertical = (2*pi**2+3*pi)*exp(-pi/sqrt(2)) + (34*pi**2+15*pi)*exp(-4*pi/sqrt(2))
    less("S19", "vertical theta-kernel bound < 4", vertical, 4)
    less("S20", "(3/(2pi))*1.585*(1+24/630) < .8", 3/(2*pi)*v("1.585")*(1+24/t), ".8")
    less("S21", "30/(1-30/630) < 32", 30/(1-30/t), 32)
    less("S22", "32/((pi/2)exp(6)) < .051", 32/(pi/2*exp(a)), ".051")
    less("S23", "2*.051*49 < 5", 2*v(".051")*49, 5)
    less("S24", "1.8*12000^2/(12000/2-1)^2 < 8", v("1.8")*n*n/(n/2-1)**2, 8)
    # Multiplying by n^2 makes the gamma bound's endpoint margin readable.
    less("S25", "2/(1-2/n)+2/(3n(1-1/n^2)) < 4 at n=12000",
         2/(1-2/n)+2/(3*n*(1-1/(n*n))), 4)
    q = exp(-log(2)/3)
    less("S26", "12/(1-2^(-1/3))^3 < 1500", 12/(1-q)**3, 1500)
    less("S27", "20 < 6000^(2/3)/10", 20, v(6000)**(v(2)/3)/10)
    less("S28", "1/4 < (1-exp(-3/4))/2", v(1)/4, (1-exp(-v(3)/4))/2)
    less("S29", "5*(20/19)^2 < 6", 5*(v(20)/19)**2, 6)
    less("S30", "5sqrt(2)/3 < 5/2", 5*sqrt(2)/3, v(5)/2)
    less("S31", "(3+sqrt(5))/2 < 8/3", (3+sqrt(5))/2, v(8)/3)
    less("S32", "8/3 < e", v(8)/3, exp(v(1)))
    less("S33", "4*17^4 exp(-3*(17^2-1)) < 1e-300", 4*v(17)**4*exp(-3*(v(17)**2-1)), "1e-300")
    less("S34", "2*12016/5+9/2 < pi exp(10)", 2*v(12016)/5+v(9)/2, pi*exp(v(10)))
    result = {"status": "passed" if all(r["passed"] for r in rows) else "failed",
              "precision_decimal_digits": iv.dps, "mpmath_version": mpmath.__version__,
              "scope": "34 scalar endpoint checks only; analytic domain and monotonicity reductions are not verified by this script",
              "endpoint_encoding": "[sign,mantissa,exponent,bitcount] = (-1)^sign * mantissa * 2^exponent",
              "checks": rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({"status": result["status"], "checks": len(rows),
                      "failed": [r["id"] for r in rows if not r["passed"]]}))
    raise SystemExit(0 if result["status"] == "passed" else 1)


if __name__ == "__main__":
    main()
