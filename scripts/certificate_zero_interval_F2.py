# /// script
# requires-python = ">=3.13"
# dependencies = ["mpmath==1.3.0"]
# ///
"""Interval certificate for Delta_2 = 0 (CIR variance, F = 2, rational eta, menu {1,2,3}).

For every real g in a closed box [g_lo, g_hi] the programme verifies the ordering
0 < gamma < a < u < zeta < b < beta, a < 2 eta/5, beta > 2 eta/3, and the switching condition
M_* - m_* <= 2 lambda of Proposition 6.1 of the paper, using exact rational inputs, outward-
rounded
interval arithmetic (mpmath, iv context), interval bisection that stops when a sign is
indeterminate (keeping a valid bracket), and an explicit tail bound for the series in E.  The join
point gamma is chosen per box from a fixed rational list.  All recorded endpoints are exact
binary rationals of the internal intervals; decimal versions are rounded outward.
Usage: certificate_zero_interval_F2.py [eta g_lo g_hi n_boxes tag]   (defaults: 3/2 23/20 7/5 160)
"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction as Fr
from pathlib import Path
import mpmath as mm

iv = mm.iv
iv.dps = 50
NSER, NBIS = 80, 80                                  # series terms, bisection steps
GAMMAS = [Fr(3, 100), Fr(1, 20), Fr(7, 100), Fr(1, 10), Fr(1, 8),
          Fr(3, 20), Fr(7, 40), Fr(1, 5), Fr(9, 40), Fr(1, 4)]
ETA = Fr(sys.argv[1]) if len(sys.argv) > 1 else Fr(3, 2)
G_LO = Fr(sys.argv[2]) if len(sys.argv) > 2 else Fr(23, 20)
G_HI = Fr(sys.argv[3]) if len(sys.argv) > 3 else Fr(7, 5)
NBOX = int(sys.argv[4]) if len(sys.argv) > 4 else 160
TAG = sys.argv[5] if len(sys.argv) > 5 else f'eta{ETA.numerator}_{ETA.denominator}'
MIN_WIDTH = Fr(1, 20480)                             # no further splitting below this width
OUT = Path(__file__).resolve().parents[1] / 'artifacts/paper_v05_certificate'
B1 = ETA - Fr(1, 2)                                  # mean reward of product 1 (F = 2, no fees)


def R(x: Fr):                     # exact rational -> point interval
    return iv.mpf(x.numerator) / x.denominator

def rat(raw):                     # raw mpf tuple (sign, man, exp, bc) -> exact Fraction
    sgn, man, e, _ = raw
    v = Fr(man) * (Fr(2) ** e)
    return -v if sgn else v

def lo(v): return rat(iv.mpf(v)._mpi_[0])            # exact lower endpoint
def hi(v): return rat(iv.mpf(v)._mpi_[1])            # exact upper endpoint

def dec(q: Fr, nd=25, up=False):  # decimal string rounded outward (down unless up=True)
    n = q * 10**nd
    k = -((-n) // 1) if up else n // 1                # ceil or floor
    s = f'{abs(int(k)):0{nd+1}d}'
    return ('-' if k < 0 else '') + s[:-nd] + '.' + s[-nd:]

def ends(v):                      # exact rationals and outward decimals of an interval
    a, b = lo(v), hi(v)
    return dict(lo=str(a), hi=str(b), dec=[dec(a), dec(b, up=True)])

def ex(x):  return iv.exp(x)
def P(x):   return 1 - ex(-2*x)*(1+2*x)
def q(x):   return 2*x*x*ex(-2*x)
eta = R(ETA); b1 = R(B1)
def G12(x): return b1 + (eta - R(Fr(3, 2)))*P(x) + R(Fr(3, 2))*q(x)   # G for the pair {1,2}
def G13(x): return b1 + 2*(eta - 2)*P(x) + 4*q(x)                       # G for the pair {1,3}
def K23(x): return (eta - R(Fr(5, 2)))*P(x) + R(Fr(5, 2))*q(x)          # G13 - G12

def series(x):
    """sum_{n>=1} (2x)^n/(n n!), tail enclosed by a_{N+1}/(1 - 2x/(N+2)); needs 2x < N+2."""
    z = 2*x
    tot = iv.mpf(0); term = iv.mpf(1)                # term = z^n/n!
    for n in range(1, NSER+1):
        term = term*z/n
        tot += term/n
    a_next = term*z/(NSER+1)/(NSER+1)
    ratio = z/(NSER+2)
    assert ratio.b < 1
    return tot + iv.mpf([0, (a_next/(1-ratio)).b])

def E_point(x): return iv.log(x) + series(x) - ex(2*x)/(2*x)          # E' = 1/q; narrow x only

def E(x):
    """E on an interval x by the mean-value form E(m) + (x-m)E'([x]), m the exact midpoint.
    (The direct formula loses the cancellation between the series and e^{2x}/(2x) on wide x.)"""
    xl, xh = lo(x), hi(x)
    if xl == xh:
        return E_point(x)
    m = R((xl + xh)/2)
    return E_point(m) + (x - m)*(ex(2*x)/(2*x*x))

W0 = lambda x: 1/(2*x) - iv.log(x)                                     # W0' = P/q - 1/q
B2, B3 = 2*eta - 2, 3*eta - R(Fr(9, 2))                                 # mean rewards of 2, 3
def Lg(x, g):                                                           # Lg' = (G13 - g)/q
    return (B3 - g)*E(x) + 2*(eta - 2)*W0(x) + 4*x
def L12(x, h):                                                          # L12' = (G12 - h)/q
    return (B2 - h)*E(x) + (eta - R(Fr(3, 2)))*W0(x) + R(Fr(3, 2))*x
def K(x):                                                               # K' = K23/q
    return (eta - R(Fr(5, 2)))*(E(x) + W0(x)) + R(Fr(5, 2))*x


def sign(v):
    if v.a > 0: return 1
    if v.b < 0: return -1
    return 0

def bisect(f, a: Fr, b: Fr):
    """Enclose the root of f in [a,b]; the endpoint signs must be definite and opposite."""
    sa, sb = sign(f(R(a))), sign(f(R(b)))
    if sa == 0 or sb == 0 or sa == sb:
        raise ValueError('no valid bracket')
    for _ in range(NBIS):
        mid = (a+b)/2
        sm = sign(f(R(mid)))
        if sm == 0:
            break                                    # indeterminate: keep the valid bracket
        if sm == sa: a = mid
        else: b = mid
    return iv.mpf([R(a).a, R(b).b])

XC13, XC12, XC23 = ETA/2, 2*ETA/3, 2*ETA/5           # maxima of G13, G12, K23

def roots_at(g, k):
    """Tight enclosures of a(g), b(g), u(g), beta(g) for a point value g (k = K23(gamma))."""
    a  = bisect(lambda x: G13(x)-g, Fr(1, 1000), XC13)
    b  = bisect(lambda x: G13(x)-g, XC13, Fr(12))
    u  = bisect(lambda x: G12(x)-(g-k), Fr(1, 1000), XC12)
    be = bisect(lambda x: G12(x)-(g-k), XC12, Fr(12))
    return a, b, u, be

def hull(p, r): return iv.mpf([min(p.a, r.a), max(p.b, r.b)])


def try_box(glo: Fr, ghi: Fr, gamma: Fr):
    """Evaluate the sufficient condition for all g in [glo, ghi] with join point gamma.

    a(g), u(g) increase and b(g), beta(g) decrease in g (monotone branches of G13, G12), so their
    ranges over the box are the hulls of the endpoint enclosures; zeta does not depend on g;
    lambda(g) is strictly decreasing (Proposition 4.2).  All other quantities are evaluated with
    g as an interval and the roots as their hull enclosures."""
    gm = R(gamma); k = K23(gm)
    g = iv.mpf([R(glo).a, R(ghi).b])
    a0, b0, u0, be0 = roots_at(R(glo), k)
    a1, b1_, u1, be1 = roots_at(R(ghi), k)
    a, b, u, be = hull(a0, a1), hull(b0, b1_), hull(u0, u1), hull(be0, be1)
    ze = bisect(lambda x: K23(x)-k, XC23, Fr(12))
    lam_at_ghi = (Lg(b1_, R(ghi))-Lg(a1, R(ghi)))/2
    lam_at_glo = (Lg(b0, R(glo))-Lg(a0, R(glo)))/2
    lam = iv.mpf([lam_at_ghi.a, lam_at_glo.b])
    Va  = L12(be, g-k) - L12(a, g-k)                 # V(a) = int_a^beta (G12-(g-k))/q
    Vu  = L12(be, g-k) - L12(u, g-k)                 # V(u)
    Vg  = Va + K(gm) - K(a) - k*(E(gm)-E(a))         # V(gamma)
    Tz  = Va - lam + K(ze) - K(a) - k*(E(ze)-E(a))   # T(zeta)
    zero = iv.mpf(0)
    # enclosures of m_* = min(V(gamma)-lambda, 0) and M_* = max(V(u), T(zeta), lambda, 0)
    m_star = iv.mpf([min((Vg-lam).a, zero.a), min((Vg-lam).b, zero.b)])
    M_star = iv.mpf([max(Vu.a, Tz.a, lam.a, zero.a), max(Vu.b, Tz.b, lam.b, zero.b)])
    margin = 2*lam - (M_star - m_star)       # encloses 2 lambda - (M_* - m_*) over the box
    x25, x23 = R(XC23), R(XC12)
    order = (gm.b < a.a) and (a.b < u.a) and (u.b < ze.a) and (ze.b < b.a) and (b.b < be.a) \
            and (a.b < x25.a) and (be.a > x23.b)
    # definite violation of the ordering for every g in the box (not merely uncertified)
    violated = (gm.a >= a.b) or (a.a >= u.b) or (u.a >= ze.b) or (ze.a >= b.b) \
               or (b.a >= be.b) or (a.a >= x25.b) or (be.b <= x23.a)
    passed = bool(order and margin.a > 0)
    hopeless = violated or (margin.b < 0)      # fails for every g in the box with this gamma
    return dict(g_lo=str(glo), g_hi=str(ghi), gamma=str(gamma), a=ends(a), u=ends(u),
                zeta=ends(ze), b=ends(b), beta=ends(be), lam=ends(lam), V_gamma=ends(Vg),
                V_u=ends(Vu), T_zeta=ends(Tz), m_star_lower=dec(lo(m_star)),
                M_star_upper=dec(hi(M_star), up=True), margin_lower=str(lo(margin)),
                margin_lower_dec=dec(lo(margin)), margin_upper_dec=dec(hi(margin), up=True),
                ordering_ok=bool(order), passed=passed), hopeless


def certify_box(glo: Fr, ghi: Fr):
    """Try the join points in order; return the first passing record, else the best failure."""
    best = None; all_hopeless = True
    for gamma in GAMMAS:
        try:
            rec, hopeless = try_box(glo, ghi, gamma)
        except (ValueError, AssertionError):   # no valid root bracket: unknown, not hopeless
            all_hopeless = False
            continue
        if rec['passed']:
            return rec, False
        all_hopeless = all_hopeless and hopeless
        if best is None or Fr(rec['margin_lower']) > Fr(best['margin_lower']):
            best = rec
    if best is None:
        best = dict(g_lo=str(glo), g_hi=str(ghi), gamma=None, passed=False, ordering_ok=False,
                    margin_lower='-1', margin_lower_dec='nan')
    return best, all_hopeless


def cost_at(g: Fr):
    """Enclosure of lambda(g) for the pair {1,3}."""
    a = bisect(lambda x: G13(x)-R(g), Fr(1, 1000), XC13)
    b = bisect(lambda x: G13(x)-R(g), XC13, Fr(12))
    return (Lg(b, R(g)) - Lg(a, R(g)))/2

def crit13():
    """crit_13 = (1/2)[L_{B1}(b0) - lim_{x->0} L_{B1}(x)], the limit being 2(2-eta)."""
    b0 = bisect(lambda x: G13(x)-b1, XC13, Fr(12))
    return (Lg(b0, b1) - 2*(2 - eta))/2, b0

def crit12():
    """crit_12 = (1/2)[L12_{B1}(b0) - lim_{x->0} L12_{B1}(x)], the limit being 3/2 - eta."""
    b0 = bisect(lambda x: G12(x)-b1, XC12, Fr(12))
    return (L12(b0, b1) - (R(Fr(3, 2)) - eta))/2, b0

def cost_at_B2():
    """lambda(g13 = B2), 3/2 < eta < 5/2: the cost beyond which holding 2 beats the {1,3} band."""
    if not (Fr(3, 2) < ETA < Fr(5, 2)):
        return None
    top = G13(R(XC13))                                # maximum of G13
    if B2.a >= top.b:
        return iv.mpf(0)                              # holding 2 beats every band
    if B2.b >= top.a:
        return None                                   # undecided
    return cost_at(2*ETA - 2)


def main():
    t0 = time.perf_counter(); OUT.mkdir(parents=True, exist_ok=True)
    width = (G_HI-G_LO)/NBOX
    queue = [(G_LO+i*width, G_LO+(i+1)*width) for i in range(NBOX)]
    done, parents, split_initial = [], 0, set()
    while queue:
        glo, ghi = queue.pop(0)
        rec, hopeless = certify_box(glo, ghi)
        if not rec['passed'] and not hopeless and (ghi-glo) > MIN_WIDTH:
            mid = (glo+ghi)/2
            queue[:0] = [(glo, mid), (mid, ghi)]
            parents += 1
            split_initial.add(int((glo-G_LO)//width))   # index of the initial box
            continue
        done.append(rec)
        print(str(glo), str(ghi), 'pass' if rec['passed'] else 'FAIL', 'gamma', rec['gamma'],
              'margin_lo', rec['margin_lower_dec'], flush=True)
    done.sort(key=lambda r: Fr(r['g_lo']))
    runs, cur = [], []                                # maximal contiguous runs of passing boxes
    for r in done:
        if r['passed'] and (not cur or Fr(cur[-1]['g_hi']) == Fr(r['g_lo'])):
            cur.append(r)
        else:
            if cur: runs.append(cur)
            cur = [r] if r['passed'] else []
    if cur: runs.append(cur)
    run = max(runs, key=lambda rr: Fr(rr[-1]['g_hi'])-Fr(rr[0]['g_lo'])) if runs else []
    g_start, g_end = (Fr(run[0]['g_lo']), Fr(run[-1]['g_hi'])) if run else (None, None)
    endpoints = {name: dict(g=str(gp), lam=ends(cost_at(gp)), lam_iv=cost_at(gp))
                 for name, gp in (('g_start', g_start), ('g_end', g_end)) if gp is not None}
    c13, b0 = (crit13() if ETA < 2 else (None, None))
    c12, b0_12 = (crit12() if ETA < Fr(3, 2) else (None, None))
    lamB2 = cost_at_B2()
    inward = ([dec(hi(endpoints['g_end']['lam_iv']), up=True),
               dec(lo(endpoints['g_start']['lam_iv']))] if run else None)   # certified for sure
    for e in endpoints.values():
        del e['lam_iv']
    rec = dict(software=dict(mpmath=mm.__version__, iv_dps=iv.dps), F=2, eta=str(ETA),
               series_terms=NSER, bisection_steps=NBIS, gammas=[str(x) for x in GAMMAS],
               g_range=[str(G_LO), str(G_HI)], initial_boxes=NBOX, min_box_width=str(MIN_WIDTH),
               terminal_boxes=len(done), subdivision_operations=parents,
               initial_boxes_split=len(split_initial), box_evaluations=len(done)+parents,
               terminal_boxes_passed=sum(r['passed'] for r in done),
               certified_g_run=[str(g_start), str(g_end)], boxes_in_run=len(run),
               min_margin_lower_over_run=(dec(min(Fr(r['margin_lower']) for r in run))
                                          if run else None),
               gammas_used_in_run=sorted({r['gamma'] for r in run}, key=Fr),
               cost_endpoint_enclosures=endpoints, certified_costs_inward=inward,
               lambda_B2_enclosure=(ends(lamB2) if lamB2 is not None else None),
               crit13_enclosure=(ends(c13) if c13 is not None else 'infinite (eta >= 2)'),
               b0_enclosure=(ends(b0) if b0 is not None else None),
               crit12_enclosure=(ends(c12) if c12 is not None else 'infinite (eta >= 3/2)'),
               b0_12_enclosure=(ends(b0_12) if b0_12 is not None else None),
               boxes_detail=done, seconds=round(time.perf_counter()-t0, 1))
    (OUT/f'certificate_F2_{TAG}.json').write_text(json.dumps(rec, indent=1))
    print(json.dumps({k: v for k, v in rec.items() if k != 'boxes_detail'}, indent=1))


if __name__ == '__main__':
    main()
