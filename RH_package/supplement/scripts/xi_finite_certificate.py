"""Directed interval certificate for a finite xi coefficient escape triangle.

Only files named finite_certificate*.json are written.  All numerical input
endpoints in those files are exact mpmath binary tuples, not rounded decimals.

The positive moment accumulation uses the two endpoint operations underlying
mpmath.iv, with explicitly directed rounding, to avoid interval-object overhead.
The moments are integrated before division by (2n)!, saving one multiplication
per node and coefficient.  This is mathematically the same positive quadrature.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TASK = ROOT / "downloads" / "riemann_xi_infinite_log_concavity"
sys.path.insert(0, str(TASK / "deps"))
sys.path.insert(0, str(ROOT))
import mpmath as mp
from mpmath.libmp import (
    fzero, mpf_add, mpf_mul, mpf_cmp, round_floor, round_ceiling,
)
from mcp.server import memory_append, branch_update

PID = "riemann_xi_infinite_log_concavity"


def iv_from_raw(raw):
    return mp.iv.make_mpf((tuple(raw[0]), tuple(raw[1])))


def point_from_raw(raw):
    return mp.iv.make_mpf((tuple(raw), tuple(raw)))


def lower(v):
    return mp.mpf(v._mpi_[0])


def upper(v):
    return mp.mpf(v._mpi_[1])


def positive(v):
    return mpf_cmp(v._mpi_[0], fzero) > 0


def negative(v):
    return mpf_cmp(v._mpi_[1], fzero) < 0


def decimal_endpoint(raw, digits, upward):
    """An outward decimal bound, using exact integer division at the end."""
    sign, man, exponent, _ = raw
    if not man:
        return "0"
    # This only selects the display scale; exact division below guarantees
    # direction even if the approximate logarithm is near an integer.
    order = int(mp.floor(mp.log10(abs(mp.mpf(raw)))))
    scale = order-digits+1
    numerator = -int(man) if sign else int(man)
    denominator = 1
    if exponent>=0:
        numerator <<= int(exponent)
    else:
        denominator <<= int(-exponent)
    if scale>=0:
        denominator *= 10**scale
    else:
        numerator *= 10**(-scale)
    mantissa = -((-numerator)//denominator) if upward else numerator//denominator
    return str(mantissa)+"e"+str(scale)


def brief(v, digits=24):
    return [decimal_endpoint(v._mpi_[0],digits,False),decimal_endpoint(v._mpi_[1],digits,True)]


def dump(name, data):
    path = TASK / name
    # With the optional gmpy backend, exact mantissas are mpz objects.
    path.write_text(json.dumps(data, separators=(",", ":"), default=int), encoding="utf-8")
    return str(path.relative_to(ROOT))


def event(data):
    memory_append(PID, "events", data)


def legendre_pair(q, x):
    p0 = mp.iv.mpf(1)
    p1 = x
    if q == 0:
        return p0, mp.iv.mpf(0)
    for k in range(1, q):
        p0, p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1)
    return p1, p0


def certified_nodes(q, dps, radius_digits=100):
    # Point endpoint evaluations need substantial guard precision because
    # naive interval recurrence propagates rounding widths pessimistically.
    mp.mp.dps = max(160, radius_digits + 50)
    mp.iv.dps = max(dps, radius_digits + q + 100)
    eps = mp.iv.mpf("1e-" + str(radius_digits))
    roots, _ = mp.gauss_quadrature(q, "legendre")
    rows = []
    previous = mp.iv.mpf(-1)
    for i, guess in enumerate(roots):
        # The decimal midpoint is treated as an exact rational, with outward
        # conversion.  Existence comes from endpoint signs, not this guess.
        center = mp.iv.mpf(mp.nstr(guess, radius_digits+30))
        left, right = center-eps, center+eps
        assert lower(left) > -1 and upper(right) < 1
        assert lower(left) > upper(previous)
        pl, _ = legendre_pair(q, left)
        pr, _ = legendre_pair(q, right)
        assert (negative(pl) and positive(pr)) or (positive(pl) and negative(pr)), (i, brief(pl), brief(pr))
        node = mp.iv.make_mpf((left._mpi_[0], right._mpi_[1]))
        pc, pm = legendre_pair(q, center)
        derivative_at_center = q*(center*pc-pm)/(center*center-1)
        # |P_q''| <= q^4 on [-1,1].  One elementary derivation is
        # P_n' = sum_{j<n, n-j odd}(2j+1)P_j and |P_j|<=1; apply twice.
        # Center is an exact decimal enclosed by center; 2*eps also absorbs
        # its tiny directed-conversion width.
        derror = 2*(q**4)*eps
        derivative = derivative_at_center + mp.iv.make_mpf(((-derror)._mpi_[0], derror._mpi_[1]))
        assert positive(derivative) or negative(derivative)
        weight = 2/((1-node*node)*derivative**2)
        assert positive(weight)
        rows.append({"node":node._mpi_, "weight":weight._mpi_, "left_sign":-1 if negative(pl) else 1, "right_sign":-1 if negative(pr) else 1})
        previous = right
    mp.iv.dps = dps
    sumweights = sum((iv_from_raw(r["weight"]) for r in rows), mp.iv.mpf(0))
    assert lower(sumweights) <= 2 <= upper(sumweights)
    return rows, {"all_brackets_disjoint":True, "all_endpoint_signs_opposite":True, "root_count":q, "degree":q, "weight_sum_interval":brief(sumweights), "root_radius":"1e-"+str(radius_digits)}


def kernel_truncated(t, theta_terms):
    u = mp.iv.pi*mp.iv.exp(2*t)
    total = mp.iv.mpf(0)
    for j in range(1, theta_terms+1):
        v = (j*j)*u
        total += (2*v*v-3*v)*mp.iv.exp(-v)
    return 4*mp.iv.exp(t/2)*total


def segment_moments(payload):
    segment, rows, nmax, dps, theta_terms = payload
    mp.iv.dps = dps
    precision = mp.iv.prec
    half = mp.iv.mpf(1)/16
    midpoint = mp.iv.mpf(2*segment+1)/16
    lows, highs = [fzero]*(nmax+1), [fzero]*(nmax+1)
    add, mul = mpf_add, mpf_mul
    down, up = round_floor, round_ceiling
    for row in rows:
        node, weight = iv_from_raw(row["node"]), iv_from_raw(row["weight"])
        t = midpoint + half*node
        assert positive(t)
        t2 = (t*t)._mpi_
        term = (half*weight*kernel_truncated(t,theta_terms))._mpi_
        assert mpf_cmp(term[0],fzero) > 0
        lo, hi = term
        tlo, thi = t2
        for n in range(nmax+1):
            lows[n] = add(lows[n], lo, precision, down)
            highs[n] = add(highs[n], hi, precision, up)
            lo = mul(lo, tlo, precision, down)
            hi = mul(hi, thi, precision, up)
    return segment, list(zip(lows,highs))


def audit_constants(nmax):
    iv = mp.iv
    radius = iv.mpf(1)/8
    alpha = iv.pi*iv.cos(iv.mpf(1)/4)
    assert positive(iv.pi*iv.exp(-iv.mpf(1)/4)-2)
    assert positive(alpha*iv.exp(-iv.mpf(1)/4)-2)
    geometric_ratio = 16*iv.exp(-6)
    # (2j^4+3j\u00b2/u)<= (7/2)j^4 when u>=2;
    # j^4 exp[-2(j\u00b2-1)] successive ratio <=16exp(-6).
    assert positive(16-(iv.mpf(7)/2)/(1-geometric_ratio))
    theta_relative = 4*(17**4)*iv.exp(-3*(17**2-1))
    assert positive(iv.mpf("1e-300")-theta_relative)
    assert positive(iv.pi*iv.exp(10)-(iv.mpf(2*nmax)/5+iv.mpf(9)/2))
    return {
        "complex_kernel_bound":"|K(x+iy)|<=64*pi\u00b2*exp(9*x/2-alpha*exp(2*x)), alpha=pi*cos(1/4), for x>=-1/8 and |y|<=1/8",
        "complex_kernel_geometric_factor":brief((iv.mpf(7)/2)/(1-geometric_ratio)),
        "theta_relative_omission_upper":brief(theta_relative),
        "tail_phase_derivative_margin":brief(iv.pi*iv.exp(10)-(iv.mpf(2*nmax)/5+iv.mpf(9)/2)),
        "gl_error_proof":"On each segment halfwidth h=1/16, Cauchy coefficients on radius R=1/8 bound Taylor tail after degree2q-1 by M*4^-q/(1-1/2). Both integral and positive Gaussian sum have mass2h, hence error<=8h*M*4^-q.40 segments give20*M*4^-q.",
        "theta_error_proof":"For real t>=0, u=pi exp(2t)>3, K_j/K_1<=2*j^4*exp[-3(j\u00b2-1)].The sum j>=17 is at most4*17^4*exp[-3(17\u00b2-1)]<1e-300.Thus truncated kernel>= (1-1e-300) full kernel.",
        "tail_proof":"For t>=5 the phase2n log t+9t/2-pi exp(2t) has derivative<=-pi exp(10), yielding tail<=64*pi\u00b2*5^(2n)*exp(45/2-pi exp10)/(pi exp10).This is for the raw moment, before division by(2n)!.",
        "phase_max_proof":"For x in[-R,R], (x\u00b2+R\u00b2)^n exp(9x/2-alpha exp2x)<= (2R\u00b2)^n exp(9R/2).For x>=R, derivative2n*x/(x\u00b2+R\u00b2)+9/2-2alpha exp2x is strictly decreasing.The maximum is at R if its derivative is nonpositive there, or at its uniquely bracketed root.Upper phase uses bracket hi in positive terms and lo in the negative exponential.",
        "node_derivative_proof":"|P_j|<=1 on[-1,1] from Laplace integral representation. P_n'=sum_{j<n,n-j odd}(2j+1)P_j implies |P_n'|<=n\u00b2 and |P_n''|<=n^4. Apply midpoint derivative enclosure with root bracket radius2eps.",
        "external_reference":"NIST DLMF3.5(v), equations3.5.18--3.5.21: Gaussian quadrature at orthogonal polynomial roots has positive weights and is exact through degree2q-1; Legendre normalization P_q(1)=1 has weights2/((1-x_i\u00b2)*P'_q(x_i)\u00b2).",
    }


def phase_error_bounds(nmax, q):
    iv = mp.iv
    radius = iv.mpf(1)/8
    alpha = iv.pi*iv.cos(iv.mpf(1)/4)
    alpha_mid = (lower(alpha)+upper(alpha))/2
    rmp = mp.mpf(1)/8
    delta = iv.mpf("1e-35")
    previous = mp.mpf(1)
    errors = []
    brackets = []
    logfactor = iv.log(64*iv.pi**2)
    logqfactor = iv.log(20)-q*iv.log(4)
    tailbase = logfactor+iv.mpf(45)/2-iv.pi*iv.exp(10)-iv.log(iv.pi*iv.exp(10))
    log5 = iv.log(5)
    for n in range(nmax+1):
        def deriv(x):
            return 2*n*x/(x*x+radius*radius)+iv.mpf(9)/2-2*alpha*iv.exp(2*x)
        at_radius = deriv(radius)
        if negative(at_radius):
            lo = hi = radius
        else:
            # Untrusted approximation followed by certified endpoint signs.
            def d(x):
                return 2*n*x/(x*x+rmp*rmp)+mp.mpf(9)/2-2*alpha_mid*mp.exp(2*x)
            approx = mp.findroot(d, (max(rmp,previous-mp.mpf('.1')), previous+mp.mpf('.1')), solver="secant")
            center = iv.mpf(mp.nstr(approx,70))
            lo, hi = center-delta, center+delta
            assert lower(lo)>=lower(radius)
            assert positive(deriv(lo)) and negative(deriv(hi)), (n,brief(deriv(lo)),brief(deriv(hi)))
            previous = approx
        lo_pt, hi_pt = point_from_raw(lo._mpi_[0]), point_from_raw(hi._mpi_[1])
        phase_b = n*iv.log(hi_pt**2+radius**2)+iv.mpf(9)/2*hi_pt-alpha*iv.exp(2*lo_pt)
        phase_a = n*iv.log(2*radius**2)+iv.mpf(9)/2*radius
        phase_up = max(phase_a._mpi_[1],phase_b._mpi_[1],key=lambda raw: mp.mpf(raw))
        bound = iv.exp(logfactor+point_from_raw(phase_up)+logqfactor)
        tail = iv.exp(tailbase+2*n*log5)
        errors.append((bound._mpi_[1],tail._mpi_[1]))
        brackets.append((lo._mpi_[0],hi._mpi_[1]))
    return errors, brackets


def enclose_coefficients(raw_moments, errors):
    iv = mp.iv
    invfact = iv.mpf(1)
    coefficients = []
    max_relative_width = iv.mpf(0)
    max_gl_relative = iv.mpf(0)
    for n,(raw,error) in enumerate(zip(raw_moments,errors)):
        approx = iv_from_raw(raw)
        gl = point_from_raw(error[0])
        tail = point_from_raw(error[1])
        a = (approx-gl)*invfact
        b = ((approx+gl)/(1-iv.mpf("1e-300"))+tail)*invfact
        value = iv.make_mpf((a._mpi_[0],b._mpi_[1]))
        assert positive(value), (n,brief(value))
        width = (point_from_raw(value._mpi_[1])-point_from_raw(value._mpi_[0]))/point_from_raw(value._mpi_[0])
        relgl = gl/approx
        if upper(width)>upper(max_relative_width):
            max_relative_width=width
        if upper(relgl)>upper(max_gl_relative):
            max_gl_relative=relgl
        coefficients.append(value._mpi_)
        invfact /= (2*n+1)*(2*n+2)
    return coefficients, {"max_initial_relative_width":brief(max_relative_width), "max_GL_relative_error":brief(max_gl_relative)}


def check_iterates(raw, steps, certify_to):
    a = [iv_from_raw(r) for r in raw]
    threshold = (3+mp.iv.sqrt(5))/2
    summaries = []
    for m in range(steps+1):
        assert all(positive(x) for x in a), ("positivity",m)
        ratios = [a[k]**2/(a[k-1]*a[k+1]) for k in range(1,min(certify_to,len(a)-2)+1)]
        least = min(range(len(ratios)),key=lambda i: lower(ratios[i]))
        summaries.append({"iteration":m,"positive_through":len(a)-1,"ratios_checked_through":len(ratios),"minimum_ratio_lower_at":least+1,"minimum_ratio_interval":brief(ratios[least]),"all_checked_ratios_above_threshold":all(positive(r-threshold) for r in ratios),"all_checked_ratios_above_10":all(positive(r-10) for r in ratios)})
        print(json.dumps(summaries[-1]),flush=True)
        if m==steps:
            assert len(ratios)>=certify_to
            passed = all(positive(r-threshold) for r in ratios)
            return summaries, [r._mpi_ for r in ratios], passed
        a = [a[0]**2]+[a[k]**2-a[k-1]*a[k+1] for k in range(1,len(a)-1)]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--nmax",type=int,default=12016)
    ap.add_argument("--q",type=int,default=256)
    ap.add_argument("--dps",type=int,default=130)
    ap.add_argument("--workers",type=int,default=8)
    ap.add_argument("--steps",type=int,default=14)
    ap.add_argument("--certify-to",type=int,default=12000)
    ap.add_argument("--tag",default="full")
    ap.add_argument("--reuse-nodes",action="store_true")
    args=ap.parse_args()
    assert 0<=args.nmax<=12016
    assert args.certify_to+args.steps+1<=args.nmax
    start=time.time()
    mp.mp.dps=160
    mp.iv.dps=args.dps
    prefix="finite_certificate_"+args.tag
    branch_update(PID,"finite_certificate",{"status":"running","parameters":vars(args)})
    audit=audit_constants(args.nmax)
    nodes_path=TASK/(prefix+"_nodes.json")
    if args.reuse_nodes and nodes_path.exists():
        saved=json.loads(nodes_path.read_text(encoding="utf-8"))
        assert saved["q"]==args.q
        rows,node_report=saved["rows"],saved["report"]
    else:
        rows,node_report=certified_nodes(args.q,args.dps)
        dump(prefix+"_nodes.json",{"q":args.q,"rows":rows,"report":node_report})
    event({"event_type":"finite_certificate_nodes","q":args.q,"report":node_report,"elapsed_seconds":time.time()-start})
    print(json.dumps({"phase":"nodes_certified","elapsed_seconds":time.time()-start,**node_report}),flush=True)
    precision=mp.iv.prec
    combined=[(fzero,fzero) for _ in range(args.nmax+1)]
    payloads=[(seg,rows,args.nmax,args.dps,16) for seg in range(40)]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(segment_moments,p) for p in payloads]
        for count,future in enumerate(as_completed(futures),1):
            segment,raw=future.result()
            combined=[(mpf_add(a[0],b[0],precision,round_floor),mpf_add(a[1],b[1],precision,round_ceiling)) for a,b in zip(combined,raw)]
            print(json.dumps({"phase":"quadrature","segments_complete":count,"last_segment":segment,"elapsed_seconds":time.time()-start}),flush=True)
    raw_path=dump(prefix+"_raw_moments.json",{"nmax":args.nmax,"q":args.q,"dps":args.dps,"moments":combined})
    print(json.dumps({"phase":"quadrature_complete","elapsed_seconds":time.time()-start}),flush=True)
    errors,brackets=phase_error_bounds(args.nmax,args.q)
    error_path=dump(prefix+"_errors.json",{"nmax":args.nmax,"q":args.q,"raw_moment_GL_and_tail_uppers":errors,"phase_max_brackets":brackets,"audit":audit})
    coefficients,error_report=enclose_coefficients(combined,errors)
    lambda_path=dump(prefix+"_lambda.json",{"nmax":args.nmax,"lambda_intervals":coefficients,"report":error_report})
    print(json.dumps({"phase":"coefficients_certified",**error_report,"elapsed_seconds":time.time()-start}),flush=True)
    summaries,ratios,passed=check_iterates(coefficients,args.steps,args.certify_to)
    if args.steps==14 and args.certify_to<=12000:
        assert summaries[-1]["all_checked_ratios_above_10"]
    ratio_path=dump(prefix+"_ratios.json",{"steps":args.steps,"certify_to":args.certify_to,"ratio_intervals":ratios})
    script_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report={"event_type":"finite_interval_escape_certificate","status":"passed" if passed else "failed","parameters":vars(args),"elapsed_seconds":time.time()-start,"node_report":node_report,"error_report":error_report,"iteration_summaries":summaries,"threshold":"(3+sqrt(5))/2","certified_rational_lower_bound":"10" if summaries[-1]["all_checked_ratios_above_10"] else None,"artifacts":{"raw_moments":raw_path,"errors":error_path,"lambda":lambda_path,"ratios":ratio_path},"script_sha256":script_hash,"endpoint_encoding":"Each endpoint is an exact binary tuple [sign,mantissa,exponent,bitcount], representing (-1)^sign*mantissa*2^exponent. Displayed decimal intervals are rounded outward using exact integer division.","trusted_arithmetic":"mpmath1.3.0 iv; positive accumulation uses its libmp mpf_add/mpf_mul with explicit round_floor/round_ceiling.","scope":"Finite escape triangle only; no whole-proof verification claim."}
    report_path=dump(prefix+"_report.json",report)
    event(report)
    branch_update(PID,"finite_certificate",{"status":report["status"],"report_path":report_path})
    print(json.dumps({"status":report["status"],"report_path":report_path,"elapsed_seconds":report["elapsed_seconds"]}),flush=True)


if __name__=="__main__":
    main()
