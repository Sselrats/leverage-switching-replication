#!/usr/bin/env python3
"""Certified cost interval for the adjacent reduction (Theorem 4.1).

For a zero-fee menu, F>1 rational, eta>0 rational, the proof of
Theorem 4.1 gives the exact identity
    g_M = B_{l_1} + sum_r (g_r - B_{l_r})
at a cost lambda whenever
  (i)  every adjacent pair r is active and its optimal band {G_r > g_r}
       lies in a compact interval I_r = [a_r, b_r] around x_r; the I_r are
       disjoint, and min(R_{l_r}, R_{l_{r+1}}) - R_{l_k} >= delta_r > 0 on I_r
       for every k outside {r, r+1};
  (ii) (b_r/F) / lambda * ((G_r(x_r) - g_r) / min_{I_r} q)^2 < delta_r.
This program certifies (i)-(ii) for EVERY lambda in (0, lambda_bar].

Scaling.  With q0 = F^(F-1)/Gamma(F), put qt = q/q0 = x^F e^(-Fx),
Gt_r = (G_r - B_{l_r})/q0 = (j-i) qt ((eta-m) Theta + m), Theta = P/q, and
gt = (g - B_{l_r})/q0.  The area equation, the band, and the residual in (ii)
are invariant under this rescaling, so Gamma(F) is not needed for (i)-(ii).
It is used only in the optional cross-check at the physical cost.

Arithmetic.  Accepted inequalities use fixed-point interval arithmetic
(integers scaled by 2^-PREC, every operation rounded outward), exact
Fractions, and explicit geometric tail bounds for every series:
  Theta(x) = sum (Fx)^n/(F+1)_n,       A(x) = int_0^x Theta,
  E(a,b)   = int_a^b x^-F e^(Fx) dx = sum F^n/n! int_a^b x^(n-F) dx,
  exp by halving + Taylor with remainder, log by atanh series,
  log Gamma by Stirling with the first-neglected-term remainder (real z>0).
Floating point only proposes integration limits, levels and grid points.

Monotonicity used (Proposition 2.3): H_r(g) = int (G_r-g)_+/q is strictly
decreasing in g, and lambda -> g_r(lambda) is strictly decreasing
(d g_r/d lambda = -2/J_r < 0).  Hence bands grow with lambda, and for any
[alpha, beta],  int_alpha^beta (G_r-g)/q <= H_r(g)  gives lower bounds.

Run: python certify.py [--output PATH] [--skip-optional]
Standard library only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import time
from fractions import Fraction as Q
from pathlib import Path

PREC = 200
ONE = 1 << PREC


# ----------------------------------------------------------------------------
# Outward-rounded fixed-point intervals
# ----------------------------------------------------------------------------
def _fl(n: int, d: int) -> int:
    return n // d


def _ce(n: int, d: int) -> int:
    return -((-n) // d)


class IV:
    """Closed interval [lo, hi] * 2^-PREC with integer endpoints."""
    __slots__ = ('lo', 'hi')

    def __init__(self, lo: int, hi: int):
        if lo > hi:
            raise ArithmeticError('empty interval')
        self.lo, self.hi = lo, hi

    @staticmethod
    def of(x) -> 'IV':
        if isinstance(x, IV):
            return x
        x = Q(x)
        n = x.numerator << PREC
        return IV(_fl(n, x.denominator), _ce(n, x.denominator))

    def lower(self) -> Q:
        return Q(self.lo, ONE)

    def upper(self) -> Q:
        return Q(self.hi, ONE)

    def __add__(self, o):
        o = IV.of(o)
        return IV(self.lo + o.lo, self.hi + o.hi)
    __radd__ = __add__

    def __neg__(self):
        return IV(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-IV.of(o))

    def __rsub__(self, o):
        return IV.of(o) - self

    def __mul__(self, o):
        o = IV.of(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return IV(min(p) >> PREC, -((-max(p)) >> PREC))
    __rmul__ = __mul__

    def recip(self):
        if self.lo > 0:
            return IV(_fl(ONE * ONE, self.hi), _ce(ONE * ONE, self.lo))
        if self.hi < 0:
            return -((-self).recip())
        raise ZeroDivisionError('interval contains zero')

    def __truediv__(self, o):
        return self * IV.of(o).recip()

    def __rtruediv__(self, o):
        return IV.of(o) * self.recip()

    def ipow(self, k: int):
        assert k >= 0
        r = IV.of(1)
        for _ in range(k):
            r = r * self
        return r

    def mag(self) -> int:
        return max(abs(self.lo), abs(self.hi))

    def widen(self, t: int):
        """Add [-t, t] (t in scaled units, t >= 0)."""
        return IV(self.lo - t, self.hi + t)

    def hull(self, o):
        return IV(min(self.lo, o.lo), max(self.hi, o.hi))

    def dec(self, places: int = 16):
        return [out_decimal(self.lower(), places, False), out_decimal(self.upper(), places, True)]


def out_decimal(x: Q, places: int, upper: bool) -> str:
    scale = 10 ** places
    n = -((-x.numerator * scale) // x.denominator) if upper else (x.numerator * scale) // x.denominator
    sign = '-' if n < 0 else ''
    n = abs(n)
    return f'{sign}{n // scale}.{n % scale:0{places}d}'


def need(cond: bool, msg: str) -> None:
    if not cond:
        raise ArithmeticError(msg)


def _tiny_mag(x: IV, e: int) -> bool:
    """True if |x| <= 2^-e."""
    return x.mag() <= (ONE >> e)


def exp_point(y: Q) -> IV:
    """Enclosure of exp(y) for rational y."""
    s = 0
    while abs(y) / (1 << s) > Q(1, 1 << 10):
        s += 1
    z = IV.of(y / (1 << s))
    term = IV.of(1)
    total = IV.of(1)
    K = 40
    for k in range(1, K + 1):
        term = term * z / k
        total = total + term
    # |remainder| <= |z|^(K+1)/(K+1)! * 1/(1-|z|) <= 2 |z|^(K+1)/(K+1)!
    zm = Q(z.mag(), ONE)
    rem = 2 * zm ** (K + 1) / math.factorial(K + 1)
    total = total.widen(_ce(rem.numerator << PREC, rem.denominator))
    need(total.lo > 0, 'exp enclosure not positive')
    for _ in range(s):
        total = total * total
    return total


def exp_iv(x: IV) -> IV:
    return IV(exp_point(x.lower()).lo, exp_point(x.upper()).hi)


def _atanh_series(z: IV, zmax: Q) -> IV:
    """atanh(z) for |z| <= zmax < 1, with remainder bound."""
    z2 = z * z
    term = z
    total = z
    k = 0
    while True:
        k += 1
        term = term * z2
        total = total + term / (2 * k + 1)
        # remainder: sum_{i>k} |z|^(2i+1)/(2i+1) <= zmax^(2k+3)/((2k+3)(1-zmax^2))
        rem = zmax ** (2 * k + 3) / ((2 * k + 3) * (1 - zmax * zmax))
        if rem < Q(1, 1 << (PREC + 8)):
            return total.widen(_ce(rem.numerator << PREC, rem.denominator))


_LOG2 = None


def log2_iv() -> IV:
    global _LOG2
    if _LOG2 is None:
        _LOG2 = 2 * _atanh_series(IV.of(Q(1, 3)), Q(1, 3))
    return _LOG2


def log_point(y: Q) -> IV:
    need(y > 0, 'log of nonpositive')
    k = 0
    m = Q(y)
    while m >= Q(3, 2):
        m /= 2
        k += 1
    while m < Q(3, 4):
        m *= 2
        k -= 1
    z = Q(m - 1, m + 1)  # |z| <= 1/5
    return k * log2_iv() + 2 * _atanh_series(IV.of(z), Q(1, 5))


def log_iv(x: IV) -> IV:
    need(x.lo > 0, 'log of interval containing nonpositive values')
    return IV(log_point(x.lower()).lo, log_point(x.upper()).hi)


def rpow(x: Q, s: Q) -> IV:
    """x^s for rational x > 0 and rational s."""
    x, s = Q(x), Q(s)
    need(x > 0, 'rpow base')
    if s.denominator == 1:
        v = x ** int(s)
        return IV.of(v)
    return exp_iv(s * log_point(x))


def atan_series(z: Q) -> IV:
    """atan(z) for 0 < z <= 1/5: alternating series, error <= next term."""
    zi = IV.of(z)
    total = IV.of(0)
    term = zi
    n = 0
    while True:
        total = total + (term / (2 * n + 1) if n % 2 == 0 else -(term / (2 * n + 1)))
        n += 1
        term = term * zi * zi
        nxt = z ** (2 * n + 1) / (2 * n + 1)
        if nxt < Q(1, 1 << (PREC + 8)):
            return total.widen(_ce(nxt.numerator << PREC, nxt.denominator))


def pi_iv() -> IV:
    return 4 * (4 * atan_series(Q(1, 5)) - atan_series(Q(1, 239)))


BERN = [Q(1, 6), Q(-1, 30), Q(1, 42), Q(-1, 30), Q(5, 66), Q(-691, 2730), Q(7, 6),
        Q(-3617, 510), Q(43867, 798), Q(-174611, 330)]


def gamma_iv(z: Q) -> IV:
    """Gamma(z), z>0 rational: shift by N, Stirling with K terms.

    For real w>0 the Stirling remainder after K terms has the sign of, and is
    bounded in magnitude by, the first neglected term (DLMF 5.11.ii)."""
    z = Q(z)
    N = 30
    w = z + N
    K = 8
    lw = log_point(w)
    s = (w - Q(1, 2)) * lw - w + log_iv(2 * pi_iv()) / 2
    for k in range(1, K + 1):
        s = s + IV.of(BERN[k - 1] / (2 * k * (2 * k - 1) * w ** (2 * k - 1)))
    rem = abs(BERN[K]) / ((2 * K + 2) * (2 * K + 1) * w ** (2 * K + 1))
    s = s.widen(_ce(rem.numerator << PREC, rem.denominator))
    prod = Q(1)
    for k in range(N):
        prod *= z + k
    return exp_iv(s) / prod


# ----------------------------------------------------------------------------
# Model
# ----------------------------------------------------------------------------
class Case:
    def __init__(self, name: str, F: Q, eta: Q, levels, lambda_bar: Q, rho, label: str):
        self.name, self.F, self.eta = name, Q(F), Q(eta)
        self.levels = [Q(v) for v in levels]
        need(self.F > 1 and self.eta > 0, 'parameters')
        need(all(a < b for a, b in zip(self.levels, self.levels[1:])), 'menu order')
        self.lambda_bar = Q(lambda_bar)
        self.rho = [Q(r) for r in rho]
        self.label = label
        self.Ff = float(self.F)
        self.etaf = float(self.eta)
        self.pairs = []
        for r in range(len(self.levels) - 1):
            i, j = self.levels[r], self.levels[r + 1]
            m = (i + j) / 2
            need(self.eta <= m, 'this program assumes eta <= m (so B = B_i, Bt = 0)')
            self.pairs.append({'r': r + 1, 'i': i, 'j': j, 'm': m, 'x': self.eta / m})
        xs = [p['x'] for p in self.pairs]
        need(all(a > b for a, b in zip(xs, xs[1:])), 'crossings must decrease')
        self._theta_cache = {}
        self._qt_cache = {}

    # -- exact reward bookkeeping --------------------------------------------
    def R(self, k: Q, x: Q) -> Q:
        return k * self.eta - k * k * x / 2

    def B(self, k: Q) -> Q:
        return k * self.eta - k * k / 2

    def delta(self, p, a: Q, b: Q) -> Q:
        """min over I=[a,b] and labels k outside the pair of min(R_i,R_j)-R_k.

        Each map x -> min(R_i,R_j)(x) - R_k(x) is a minimum of affine
        functions, hence concave; its minimum over [a,b] is at an endpoint."""
        vals = []
        for k in self.levels:
            if k in (p['i'], p['j']):
                continue
            for x in (a, b):
                vals.append(min(self.R(p['i'], x), self.R(p['j'], x)) - self.R(k, x))
        return min(vals) if vals else Q(10 ** 9)

    # -- certified primitives --------------------------------------------------
    def qt(self, x: Q) -> IV:
        """qt(x) = x^F e^(-F x)."""
        x = Q(x)
        if x not in self._qt_cache:
            if self.F.denominator == 1:
                v = IV.of(x ** int(self.F)) * exp_point(-self.F * x)
            else:
                v = exp_iv(self.F * (log_point(x) - x))
            self._qt_cache[x] = v
        return self._qt_cache[x]

    def theta_A(self, x: Q):
        """Enclosures of Theta(x)=P/q and A(x)=int_0^x Theta."""
        x = Q(x)
        if x in self._theta_cache:
            return self._theta_cache[x]
        F = self.F
        xi = IV.of(x)
        Fx = IV.of(F * x)
        t = IV.of(1)
        th = IV.of(1)
        A = IV.of(x)
        n = 0
        while True:
            n += 1
            t = t * Fx / (F + n)
            th = th + t
            A = A + t * xi / (n + 1)
            r = F * x / (F + n + 1)
            if r < Q(1, 2):
                tu = t.upper()
                tail_th = tu * r / (1 - r)
                tail_A = tu * r * x / (n + 2) / (1 - r)
                if tail_th < Q(1, 1 << (PREC - 8)) and tail_A < Q(1, 1 << (PREC - 8)):
                    th = IV(th.lo, th.hi + _ce(tail_th.numerator << PREC, tail_th.denominator))
                    A = IV(A.lo, A.hi + _ce(tail_A.numerator << PREC, tail_A.denominator))
                    self._theta_cache[x] = (th, A)
                    return th, A
            need(n < 5000, 'Theta series did not converge')

    def Gt(self, p, x: Q) -> IV:
        th, _ = self.theta_A(x)
        return (p['j'] - p['i']) * self.qt(x) * ((self.eta - p['m']) * th + p['m'])

    def Eint(self, a: Q, b: Q) -> IV:
        """int_a^b x^-F e^(Fx) dx, 0<a<b, positive series with tail."""
        a, b = Q(a), Q(b)
        need(0 < a < b, 'Eint limits')
        F = self.F
        am = rpow(a, -F)
        bm = rpow(b, -F)
        total = IV.of(0)
        coef = Q(1)  # F^n/n!
        n = 0
        while True:
            e = n + 1 - F
            if e == 0:
                piece = log_point(b) - log_point(a)
            else:
                piece = (IV.of(b ** (n + 1)) * bm - IV.of(a ** (n + 1)) * am) / e
            total = total + coef * piece
            n += 1
            coef = coef * F / n
            r = F * b / (n + 1)
            if r < Q(1, 2):
                # tail sum_{k>=n} F^k/k! int x^(k-F) <= (b-a) a^-F (Fb)^n/n! /(1-r)
                tail = (b - a) * am.upper() * coef * b ** n / (1 - r)
                if tail < Q(1, 1 << (PREC - 8)):
                    return total + IV(0, _ce(tail.numerator << PREC, tail.denominator))
            need(n < 5000, 'E series did not converge')

    def integral(self, p, g: Q, a: Q, b: Q) -> IV:
        """int_a^b (Gt - g)/qt dx = (j-i)(eta-m)(A(b)-A(a)) + (j-i) m (b-a) - g E(a,b)."""
        d = p['j'] - p['i']
        _, Aa = self.theta_A(a)
        _, Ab = self.theta_A(b)
        return d * (self.eta - p['m']) * (Ab - Aa) + IV.of(d * p['m'] * (b - a)) - IV.of(g) * self.Eint(a, b)

    # -- floating-point helpers (proposals only) -------------------------------
    def theta_f(self, x: float) -> float:
        F = self.Ff
        s = t = 1.0
        n = 0
        while True:
            n += 1
            t *= F * x / (F + n)
            s += t
            if t < 1e-18 * s and F * x / (F + n + 1) < 0.5:
                return s

    def A_f(self, x: float) -> float:
        F = self.Ff
        t = 1.0
        s = x
        n = 0
        while True:
            n += 1
            t *= F * x / (F + n)
            s += t * x / (n + 1)
            if t < 1e-18 and F * x / (F + n + 1) < 0.5:
                return s

    def qt_f(self, x: float) -> float:
        return math.exp(self.Ff * (math.log(x) - x))

    def Gt_f(self, p, x: float) -> float:
        return float(p['j'] - p['i']) * self.qt_f(x) * (float(self.eta - p['m']) * self.theta_f(x) + float(p['m']))

    def E_f(self, a: float, b: float) -> float:
        F = self.Ff
        s = 0.0
        c = 1.0
        n = 0
        while True:
            e = n + 1 - F
            s += c * (math.log(b / a) if abs(e) < 1e-15 else (b ** e - a ** e) / e)
            n += 1
            c *= F / n
            if c * b ** (n + 1) < 1e-19 * max(s, 1.0) and F * b / (n + 1) < 0.5:
                return s

    def roots_f(self, p, g: float):
        xr = float(p['x'])
        f = lambda y: self.Gt_f(p, y) - g
        lo, hi = 1e-9, xr
        for _ in range(100):
            mid = (lo + hi) / 2
            if f(mid) > 0:
                hi = mid
            else:
                lo = mid
        u = (lo + hi) / 2
        lo, hi = xr, xr * 4 + 2
        need(f(hi) < 0, 'right bracket')
        for _ in range(100):
            mid = (lo + hi) / 2
            if f(mid) > 0:
                lo = mid
            else:
                hi = mid
        return u, (lo + hi) / 2

    def H_f(self, p, g: float) -> float:
        u, w = self.roots_f(p, g)
        d = float(p['j'] - p['i'])
        return d * float(self.eta - p['m']) * (self.A_f(w) - self.A_f(u)) + d * float(p['m']) * (w - u) - g * self.E_f(u, w)

    def g_f(self, p, lam: float) -> float:
        M = self.Gt_f(p, float(p['x']))
        lo, hi = 0.0, M
        for _ in range(90):
            mid = (lo + hi) / 2
            if self.H_f(p, mid) > 2 * lam:
                lo = mid
            else:
                hi = mid
        return lo

    # -- certified lower bound on g_r(lambda) ----------------------------------
    def g_lower(self, p, lam: Q):
        """Rational gh with H_r(gh) >= 2 lam certified, hence g_r(lam) >= gh."""
        gf = self.g_f(p, float(lam))
        M = self.Gt_f(p, float(p['x']))
        for back in (1e-11, 1e-9, 1e-7, 1e-5):
            gh = Q(gf - back * M).limit_denominator(10 ** 18)
            if gh <= 0:
                gh = Q(0)
                break
            u, w = self.roots_f(p, float(gh))
            a, b = Q(u).limit_denominator(10 ** 16), Q(w).limit_denominator(10 ** 16)
            Hl = self.integral(p, gh, a, b).lower()
            if Hl >= 2 * lam:
                return gh, (a, b), Hl
        raise ArithmeticError('could not certify a lower bound on g')

    # -- certified upper bound on g_r(lambda) (optional cross-check) ------------
    def g_upper(self, p, lam_hi: Q):
        """Rational gh with H_r(gh) < 2 lam_lo certified, hence g_r(lam) < gh."""
        gf = self.g_f(p, float(lam_hi))
        M = self.Gt_f(p, float(p['x']))
        for fwd in (1e-11, 1e-9, 1e-7):
            gh = Q(gf + fwd * M).limit_denominator(10 ** 18)
            u, w = self.roots_f(p, float(gh))
            eps = Q(1, 10 ** 12)
            a, b = Q(u).limit_denominator(10 ** 16), Q(w).limit_denominator(10 ** 16)
            need(0 < a - eps and a + eps < p['x'] < b - eps, 'root cells must avoid the crossing')
            # root enclosures: Gt increasing on (0,x_r), decreasing after
            if not (self.Gt(p, a - eps).upper() < gh < self.Gt(p, a + eps).lower()):
                continue
            if not (self.Gt(p, b + eps).upper() < gh < self.Gt(p, b - eps).lower()):
                continue
            corr = Q(0)
            for c0, c1 in ((a - eps, a + eps), (b - eps, b + eps)):
                gm = max(abs(self.Gt(p, c0).lower() - gh), abs(self.Gt(p, c0).upper() - gh),
                         abs(self.Gt(p, c1).lower() - gh), abs(self.Gt(p, c1).upper() - gh))
                qm = min(self.qt(c0).lower(), self.qt(c1).lower())
                corr += eps * gm / qm
            Hu = self.integral(p, gh, a, b).upper() + corr
            if Hu < 2 * lam_hi:
                return gh, Hu
        raise ArithmeticError('could not certify an upper bound on g')

    # -- sup bounds on a cell for the tail --------------------------------------
    def qt_sup(self, s: Q, t: Q) -> Q:
        """sup of qt on [s,t]; qt is unimodal with maximum at x=1."""
        if s <= 1 <= t:
            return self.qt(Q(1)).upper()
        return max(self.qt(s).upper(), self.qt(t).upper())

    def G2_sup(self, p, s: Q, t: Q) -> Q:
        """sup |Gt''| on [s,t]: Gt'' = (j-i) ft [((F-1)/x - F)(eta - m x) - m], ft = F qt/x."""
        F = self.F
        x = IV(IV.of(s).lo, IV.of(t).hi)
        br = ((F - 1) / x - F) * (self.eta - p['m'] * x) - p['m']
        ft = F * self.qt_sup(s, t) / s
        return (p['j'] - p['i']) * ft * Q(br.mag(), ONE)


# ----------------------------------------------------------------------------
# Certification of one case
# ----------------------------------------------------------------------------
def fr(x: Q) -> str:
    return f'{x.numerator}/{x.denominator}' if x.denominator != 1 else str(x.numerator)


def certify_case(C: Case, grid_den: int = 10 ** 6):
    t0 = time.time()
    lam_bar = C.lambda_bar
    out = {'name': C.name, 'label': C.label, 'F': fr(C.F), 'eta': fr(C.eta),
           'menu': [fr(v) for v in C.levels], 'lambda_bar': fr(lam_bar),
           'lambda_bar_decimal': out_decimal(lam_bar, 12, False), 'pairs': []}
    intervals = []
    for p in C.pairs:
        # (1) choose I_r from the floating band at lambda_bar, then certify
        gf = C.g_f(p, float(lam_bar))
        M = C.Gt_f(p, float(p['x']))
        gam = Q(gf - 1e-7 * M).limit_denominator(10 ** 15)       # gamma_r < g_r(lambda_bar)
        uf, wf = C.roots_f(p, float(gam))
        a = Q(math.floor(uf * grid_den), grid_den)
        b = Q(math.ceil(wf * grid_den), grid_den)
        need(0 < a < p['x'] < b, 'I_r must surround the crossing')
        Ga, Gb = C.Gt(p, a), C.Gt(p, b)
        need(Ga.upper() <= gam and Gb.upper() <= gam, 'Gt(a_r), Gt(b_r) must be <= gamma_r')
        # H_r(gamma_r) >= 2 lambda_bar  =>  g_r(lambda_bar) >= gamma_r > 0 = Bt
        ua, wa = C.roots_f(p, float(gam))
        ai, bi = Q(ua).limit_denominator(10 ** 16), Q(wa).limit_denominator(10 ** 16)
        H_gam = C.integral(p, gam, ai, bi).lower()
        need(gam > 0, 'activity requires gamma_r > Bt = 0')
        need(H_gam >= 2 * lam_bar, 'band containment at lambda_bar fails')
        delta = C.delta(p, a, b)
        need(delta > 0, 'reward gap delta_r must be positive')
        qmin = min(C.qt(a).lower(), C.qt(b).lower())       # qt unimodal, max at 1
        Mhi = C.Gt(p, p['x']).upper()
        supx = b / C.F
        # (2) residual on [lambda_a, lambda_bar]: greedy certified partition
        #     Phi(lam) <= supx * (Mhi - gh_k)^2 / (lam_{k-1} qmin^2) on [lam_{k-1}, lam_k]
        # (3) tail (0, lambda_a]: explicit bound (see memo)
        rho = C.rho[p['r'] - 1]
        xr = p['x']
        ncell = 400
        K2 = Q(0)
        Qs = Q(0)
        for c in range(ncell):
            s = xr - rho + 2 * rho * c / ncell
            t = xr - rho + 2 * rho * (c + 1) / ncell
            need(s > 0, 'tail neighbourhood must be positive')
            K2 = max(K2, C.G2_sup(p, s, t))
            Qs = max(Qs, C.qt_sup(s, t))
        Kc = K2 / 2
        Ccoef = Q(81, 16) * supx ** 3 * Kc ** 2 * Qs ** 4 / qmin ** 6
        lam_a_cap1 = delta ** 3 / Ccoef                  # Phi^3 <= Ccoef * lam
        lam_a_cap2 = Q(2, 3) * Kc * rho ** 3 / Qs          # 2 lam < (4/3) K rho^3 / Q
        lam_a_raw = min(lam_a_cap1, lam_a_cap2) * Q(999, 1000)
        lam_a = Q(math.floor(lam_a_raw * 10 ** 15), 10 ** 15)
        need(0 < lam_a, 'tail start must be positive')
        need(Ccoef * lam_a < delta ** 3, 'tail residual')
        need(2 * lam_a < Q(4, 3) * Kc * rho ** 3 / Qs, 'tail validity')
        tail_bound_cubed = Ccoef * lam_a
        # greedy partition downward from lambda_bar
        nodes = []
        lam_k = lam_bar
        while True:
            gh, lims, Hl = C.g_lower(p, lam_k)
            num = supx * (Mhi - gh) ** 2 / (delta * qmin ** 2)   # need lam_{k-1} > num
            need(num < lam_k, f'residual already >= delta at lambda={float(lam_k)}')
            lam_prev = num + min(num / 10000, (lam_k - num) / 2)
            lam_prev = Q(math.ceil(lam_prev * 10 ** 18), 10 ** 18)
            if lam_prev <= lam_a:
                lam_prev = lam_a
            phi = supx * (Mhi - gh) ** 2 / (lam_prev * qmin ** 2)
            need(phi < delta, 'subinterval residual')
            need(lam_prev < lam_k, 'partition stalled')
            nodes.append({'lambda_lo': fr(lam_prev), 'lambda_hi': fr(lam_k),
                          'lambda_lo_dec': out_decimal(lam_prev, 12, False),
                          'lambda_hi_dec': out_decimal(lam_k, 12, True),
                          'g_lower_at_lambda_hi': fr(gh), 'integration_limits': [fr(lims[0]), fr(lims[1])],
                          'H_lower': out_decimal(Hl, 16, False),
                          'residual_upper': fr(phi), 'residual_upper_dec': out_decimal(phi, 8, True),
                          'margin_delta_minus_residual_dec': out_decimal(delta - phi, 8, False)})
            if not nodes[:-1]:
                phi_bar = supx * (Mhi - gh) ** 2 / (lam_bar * qmin ** 2)
            if lam_prev == lam_a:
                break
            lam_k = lam_prev
        rec = {'pair_index': p['r'], 'levels': [fr(p['i']), fr(p['j'])], 'crossing_x': fr(xr),
               'I': [fr(a), fr(b)], 'I_dec': [out_decimal(a, 6, False), out_decimal(b, 6, True)],
               'delta': fr(delta), 'delta_dec': out_decimal(delta, 8, False),
               'gamma_band_level': fr(gam), 'Gt_at_a_upper': out_decimal(Ga.upper(), 16, True),
               'Gt_at_b_upper': out_decimal(Gb.upper(), 16, True),
               'H_gamma_lower': out_decimal(H_gam, 16, False), 'two_lambda_bar': out_decimal(2 * lam_bar, 16, True),
               'containment_integration_limits': [fr(ai), fr(bi)],
               'Gt_crossing_upper': out_decimal(Mhi, 16, True), 'qt_min_on_I_lower': out_decimal(qmin, 16, False),
               'sup_x_over_F': fr(supx),
               'residual_at_lambda_bar_upper_dec': out_decimal(phi_bar, 8, True),
               'residual_top_subinterval_upper_dec': nodes[0]['residual_upper_dec'],
               'partition_subintervals': len(nodes), 'partition': nodes,
               'tail': {'lambda_a': fr(lam_a), 'lambda_a_dec': out_decimal(lam_a, 12, False),
                        'rho': fr(rho), 'cells': ncell,
                        'K_upper_half_sup_abs_G2': out_decimal(Kc, 12, True),
                        'Q_upper_sup_qt': out_decimal(Qs, 12, True),
                        'C_coefficient': out_decimal(Ccoef, 12, True),
                        'residual_cubed_bound_at_lambda_a': out_decimal(tail_bound_cubed, 12, True),
                        'delta_cubed': out_decimal(delta ** 3, 12, False),
                        'validity_2lambda_a': out_decimal(2 * lam_a, 12, True),
                        'validity_rhs_4K_rho3_over_3Q': out_decimal(Q(4, 3) * Kc * rho ** 3 / Qs, 12, False)}}
        rec['tail']['residual_bound_at_lambda_a_upper'] = f'{float(tail_bound_cubed) ** (1 / 3):.6f} (cube root of certified cube bound, display only)'
        out['pairs'].append(rec)
        intervals.append((a, b))
    # disjointness: crossings decrease, so I_{r+1} must lie left of I_r
    for (a1, b1), (a2, b2) in zip(intervals, intervals[1:]):
        need(b2 < a1, 'I_r must be disjoint')
    out['disjoint_I'] = True
    out['binding_pair'] = min(out['pairs'], key=lambda r: Q(r['delta']) - Q(r['partition'][0]['residual_upper']))['pair_index']
    out['seconds'] = round(time.time() - t0, 2)
    return out


# ----------------------------------------------------------------------------
# Floating-point exploration (not certified; reported as estimates)
# ----------------------------------------------------------------------------
def float_limits(C: Case):
    """Largest lambda (float) for which the paper's check passes with I_r = band(lambda)."""
    def ok(lam):
        bands = []
        for p in C.pairs:
            g = C.g_f(p, lam)
            u, w = C.roots_f(p, g)
            M = C.Gt_f(p, float(p['x']))
            a, b = Q(u), Q(w)
            d = float(C.delta(p, a, b))
            phi = (w / C.Ff) / lam * ((M - g) / min(C.qt_f(u), C.qt_f(w))) ** 2
            if not phi < d:
                return False
            bands.append((u, w))
        return all(b2 < a1 for (a1, b1), (a2, b2) in zip(bands, bands[1:]))
    lo, hi = 1e-7, 1.0
    need(ok(lo), 'float check fails at 1e-7')
    for _ in range(50):
        mid = math.sqrt(lo * hi)
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo


def float_pointwise_limit(C: Case):
    """Largest lambda (float) for which the pointwise version of the smoothing
    inequality holds, with I_r = band(lambda):
       R_k - [(1-psi')R_i + psi' R_j] + (x/F) psi''(d) d'^2 <= 0 on each band,
    plus disjoint bands.  An estimate of how far THIS construction can reach."""
    def ok(lam):
        bands = []
        for p in C.pairs:
            g = C.g_f(p, lam)
            u, w = C.roots_f(p, g)
            bands.append((u, w))
            n = 3000
            h = (w - u) / n
            dval = lam
            prev = 0.0
            for s in range(n + 1):
                x = u + s * h
                dp = (g - C.Gt_f(p, x)) / C.qt_f(x)
                if s > 0:
                    dval += 0.5 * (prev + dp) * h
                prev = dp
                eps = lam / 2
                if dval >= eps:
                    psi1, psi2 = 1.0, 0.0
                elif dval <= -eps:
                    psi1, psi2 = 0.0, 0.0
                else:
                    psi1, psi2 = (dval + eps) / (2 * eps), 1 / lam
                xf = Q(x)
                Ri, Rj = float(C.R(p['i'], xf)), float(C.R(p['j'], xf))
                for k in C.levels:
                    if k in (p['i'], p['j']):
                        continue
                    res = float(C.R(k, xf)) - ((1 - psi1) * Ri + psi1 * Rj) + (x / C.Ff) * psi2 * dp * dp
                    if res > 0:
                        return False
        return all(b2 < a1 for (a1, b1), (a2, b2) in zip(bands, bands[1:]))
    lo, hi = 1e-6, 1.0
    for _ in range(30):
        mid = math.sqrt(lo * hi)
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo


def float_overlap_limit(C: Case):
    """Smallest lambda (float) at which two adjacent optimal bands meet."""
    def disjoint(lam):
        bands = [C.roots_f(p, C.g_f(p, lam)) for p in C.pairs]
        return all(b2 < a1 for (a1, b1), (a2, b2) in zip(bands, bands[1:]))
    lo, hi = 1e-6, 1.0
    for _ in range(40):
        mid = math.sqrt(lo * hi)
        if disjoint(mid):
            lo = mid
        else:
            hi = mid
    return lo


# ----------------------------------------------------------------------------
# Physical cost of the paper's numerical example
# ----------------------------------------------------------------------------
def physical_cost_check(C: Case, root: Path):
    kappa, theta = Q(31, 10), Q(1, 25)
    lam = (kappa / theta) * (-log_point(Q(4991, 5000)))   # lambda = kappa d / theta
    lam_lo, lam_hi = lam.lower(), lam.upper()
    res = {'lambda_phys': lam.dec(15), 'lambda_bar_over_lambda_phys_upper': out_decimal(C.lambda_bar / lam_lo, 6, True)}
    need(C.lambda_bar < lam_lo, 'unexpected: lambda_phys inside certified interval')
    res['inside_certified_interval'] = False
    # upper bounds on g_r(lambda_phys): H_r(gh) < 2 lam_lo  => g_r < gh
    ups = []
    for p in C.pairs:
        gh, Hu = C.g_upper(p, lam_lo)
        ups.append(gh)
    # bands contain {Gt_r > gh_r}; show x0 = 3/4 lies in both (overlap of adjacent bands)
    x0 = Q(3, 4)
    inside = [C.Gt(p, x0).lower() > gh for p, gh in zip(C.pairs, ups)]
    need(all(inside), 'overlap witness failed')
    res['overlap_witness_x'] = '3/4'
    res['both_optimal_bands_contain_witness'] = True
    res['g_tilde_upper'] = [out_decimal(g, 14, True) for g in ups]
    # adjacent sum B_1 + sum (g_r - B_r) = B_1 + q0 * sum gt_r  (B_{l_r} subtracted)
    q0 = rpow(C.F, C.F - 1) / gamma_iv(C.F)
    res['q0'] = q0.dec(15)
    res['gamma_F'] = gamma_iv(C.F).dec(15)
    adj_hi = C.B(C.levels[0]) + q0.upper() * sum(ups)
    res['adjacent_sum_upper'] = out_decimal(adj_hi, 12, True)
    rec_path = root / 'artifacts/cir_fullmenu_endpoints_20260927/full_menu_certificate.json'
    if rec_path.exists():
        d = json.loads(rec_path.read_text())
        nom = d['cases']['nominal']
        glo = Q(nom['gain_dim']['lo'])
        res['external_record'] = {'file_sha256': hashlib.sha256(rec_path.read_bytes()).hexdigest(),
                                  'branch': nom['branch'], 'g_M_lower': out_decimal(glo, 12, False),
                                  'record_lambda': d['parameters']['lambda']['dec']}
        res['identity_fails_at_lambda_phys_given_record'] = bool(glo > adj_hi)
        res['g_M_minus_adjacent_sum_lower'] = out_decimal(glo - adj_hi, 10, False)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    here = Path(__file__).resolve()
    root = here.parents[2]
    ap.add_argument('--output', type=Path, default=root / 'artifacts/adjacent_interval_20260929/record.json')
    ap.add_argument('--skip-optional', action='store_true')
    args = ap.parse_args()
    T0 = time.time()
    main_case = Case('menu123_F124over45_eta3over2', Q(124, 45), Q(3, 2), [1, 2, 3],
                     lambda_bar=Q(233, 100000), rho=[Q(1, 20), Q(1, 25)],
                     label='representative case (paper numerical example parameters)')
    cases = [main_case]
    if not args.skip_optional:
        cases.append(Case('menu1_3half_3_F2_eta61over50', Q(2), Q(61, 50), [1, Q(3, 2), 3],
                          lambda_bar=Q(9, 40000), rho=[Q(1, 20), Q(1, 50)],
                          label='optional case (Example 3.2 menu)'))
    records = []
    for C in cases:
        print(f'== {C.name}', flush=True)
        rec = certify_case(C)
        rec['float_estimates_not_certified'] = {
            'paper_check_limit_with_I_equal_band': f'{float_limits(C):.6g}',
            'pointwise_smoothing_limit_with_I_equal_band': f'{float_pointwise_limit(C):.6g}',
            'adjacent_bands_first_meet': f'{float_overlap_limit(C):.6g}'}
        if C is main_case:
            rec['physical_cost'] = physical_cost_check(C, root)
        records.append(rec)
        print(f"  lambda_bar = {rec['lambda_bar']} = {rec['lambda_bar_decimal']}")
        for pr in rec['pairs']:
            print(f"  pair {pr['levels']}: I={pr['I_dec']} delta={pr['delta_dec']} "
                  f"residual(lambda_bar)<={pr['residual_at_lambda_bar_upper_dec']} "
                  f"subintervals={pr['partition_subintervals']} lambda_a={pr['tail']['lambda_a_dec']}")
        print(f"  binding pair: {rec['binding_pair']}; float estimates: {rec['float_estimates_not_certified']}")
        if 'physical_cost' in rec:
            pc = rec['physical_cost']
            print(f"  lambda_phys = {pc['lambda_phys']}; inside certified interval: {pc['inside_certified_interval']}")
            print(f"  optimal bands overlap at x=3/4: {pc['both_optimal_bands_contain_witness']}")
            if 'identity_fails_at_lambda_phys_given_record' in pc:
                print(f"  adjacent sum <= {pc['adjacent_sum_upper']} < g_M >= {pc['external_record']['g_M_lower']}: "
                      f"{pc['identity_fails_at_lambda_phys_given_record']}")
    result = {'description': __doc__, 'python_version': platform.python_version(),
              'arithmetic': f'fixed-point outward-rounded intervals, 2^-{PREC} grid; exact Fractions for acceptance',
              'program_sha256': hashlib.sha256(here.read_bytes()).hexdigest(),
              'runtime_seconds': round(time.time() - T0, 1), 'cases': records}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print('Saved', args.output, f'({result["runtime_seconds"]} s)')


if __name__ == '__main__':
    main()
