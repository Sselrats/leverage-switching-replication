"""Replay the opening example with outward interval arithmetic.

This is a separate F=2 positive-series implementation, using exact rational
inputs. Point arithmetic is used only to propose critical-root brackets.
"""
from pathlib import Path
from fractions import Fraction as Q
import json
from mpmath import mp, iv

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts/v40_rewrite_20260928"
mp.dps = 80
iv.dps = 80
N = 180


def I(x):
    x = Q(x)
    return iv.mpf(x.numerator) / x.denominator


def frac(endpoint):
    sign, man, exponent, bc = endpoint
    return Q((-1 if sign else 1) * man) * Q(2) ** exponent


def lo(x):
    return frac(x._mpi_[0])


def hi(x):
    return frac(x._mpi_[1])


def enc(x):
    return {"lo": str(lo(x)), "hi": str(hi(x)),
            "display": [mp.nstr(mp.mpf(lo(x).numerator)/lo(x).denominator, 50),
                        mp.nstr(mp.mpf(hi(x).numerator)/hi(x).denominator, 50)]}


def E(x):
    # sum (2x)^n/(n n!), with a geometric upper bound on its positive tail
    total = iv.mpf(0)
    factorial_term = iv.mpf(1)
    for n in range(1, N + 1):
        factorial_term *= 2*x/n
        total += factorial_term/n
    next_term = factorial_term * (2*x)/(N+1)**2
    ratio = 2*x/(N+2)
    assert hi(ratio) < 1
    tail = next_term/(1-ratio)
    enclosure = iv.mpf([0, str(hi(tail).numerator)]) / hi(tail).denominator
    return iv.log(x) + total + enclosure - iv.exp(2*x)/(2*x)


def A(x):
    return E(x) + 1/(2*x) - iv.log(x)


def theta(x):
    return (iv.exp(2*x)-1-2*x)/(2*x*x)


eta = Q(61, 50)
pairs = []
costs = []
for i,j,bracket in [(Q(1),Q(3,2),(Q(3),Q(4))),
                    (Q(1),Q(3),(Q(1),Q(2))),
                    (Q(3,2),Q(3),(Q(9,10),Q(11,10)))]:
    m=(i+j)/2
    # Locate the positive zero numerically, then check a rational enclosure.
    mm=mp.mpf(m.numerator)/m.denominator
    ee=mp.mpf(eta.numerator)/eta.denominator
    f=lambda x:mm-(mm-ee)*(mp.exp(2*x)-1-2*x)/(2*x*x)
    lower,upper=mp.mpf('1e-8'),mp.mpf(8)
    for _ in range(300):
        mid=(lower+upper)/2
        if f(mid)>0: lower=mid
        else: upper=mid
    center=Q(str((lower+upper)/2))
    radius=Q(1,10**45)
    zl,zh=center-radius,center+radius
    h=lambda x:I(m)-(I(m)-I(eta))*theta(x)
    assert lo(h(I(zl)))>0 and hi(h(I(zh)))<0
    # Construct the hull from rational intervals without binary-float conversion.
    z=iv.mpf([I(zl).a,I(zh).b])
    cost=I(j-i)*(I(m)*z-I(m-eta)*(A(z)+1))/2
    assert hi(cost)<Q(89,100)
    costs.append(cost)
    pairs.append({"pair":[str(i),str(j)],"root":enc(z),"cost":enc(cost)})

path=2*(costs[0]+costs[2])/3
assert lo(path)>Q(89,100)

entry=I(Q(1,10000))
z=I(Q(497477,500000))
w=I(Q(3455291,1000000))
def edge(i,j,end):
    m=(i+j)/2
    return I(j-i)*(I(m)*(end-entry)-I(m-eta)*(A(end)-A(entry)))
cycle=(edge(Q(1),Q(3,2),w)+edge(Q(3,2),Q(3),z)-3*I(Q(89,100)))/(E(w)-E(entry))
assert Q('0.00001789780378448029607') < lo(cycle)
assert hi(cycle) < Q('0.00001789780378448029608')
record={"status":"passed","F":"2","eta":str(eta),"lambda":"89/100",
        "precision":80,"positive_series_terms":N,"pairs":pairs,
        "path_threshold":enc(path),"cycle_excess":enc(cycle)}
OUT.mkdir(parents=True,exist_ok=True)
(OUT/"cycle_certificate.json").write_text(json.dumps(record,indent=2)+"\n")
print("All pairs inactive; full cycle positive:",record["cycle_excess"]["display"])
