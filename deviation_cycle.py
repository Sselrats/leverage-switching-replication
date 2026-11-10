"""Two numerical illustrations of the improving cycle in the deletion test.

Run with Python 3 and SciPy. The calculations use floating-point quadrature
and root finding. The nominal full-menu interval certificate is independent.
"""
from math import exp, log1p

from scipy.integrate import quad
from scipy.optimize import brentq, root
from scipy.special import gamma, gammainc, gammaincc



def integral(fun, start, end):
    return quad(fun, start, end, epsabs=2e-12, epsrel=2e-12)[0]


def calculate(shape, eta, cost=None, prescribed_outer_gain=None):
    def q(x):
        return x**shape * shape**(shape - 1) * exp(-shape * x) / gamma(shape)


    def G(i, j, x):
        midpoint = (i + j) / 2
        return i * eta - i**2 / 2 + (j - i) * (
            (eta - midpoint) * gammainc(shape, shape * x) + midpoint * q(x)
        )


    def K(x):
        return G(2, 3, x) - (2 * eta - 2)


    def outer_band(gain):
        peak = eta / 2
        level = lambda x: G(1, 3, x) - gain
        return brentq(level, 0, peak), brentq(level, peak, 8)


    def outer_area(gain):
        a, b = outer_band(gain)
        return integral(lambda x: (G(1, 3, x) - gain) / q(x), a, b)


    def green(start, end, label=None, upward=True):
        """Up/down Green primitive difference. label=None integrates time."""
        def density(x):
            mass = gammainc(shape, shape * x) if upward else gammaincc(shape, shape * x)
            if label is None:
                return mass / q(x)
            first_moment = mass - q(x) if upward else mass + q(x)
            return (label * eta * mass - label**2 / 2 * first_moment) / q(x)
        return integral(density, start, end)


    if prescribed_outer_gain is None:
        g = brentq(lambda h: outer_area(h) - 2 * cost, 1.00001, G(1, 3, eta / 2) - 1e-8)
    else:
        g = prescribed_outer_gain
        cost = outer_area(g) / 2
    a, b = outer_band(g)

    def fit(endpoints):
        low, high = endpoints
        k = K(low)
        area = integral(lambda x: (G(1, 2, x) - g + k) / q(x), a, high)
        area += integral(lambda x: (K(x) - k) / q(x), a, low)
        return [area - cost, G(1, 2, high) - g + k]

    if prescribed_outer_gain is None:
        solution = root(fit, [0.24, 1.565], tol=1e-12)
        if not solution.success or max(map(abs, solution.fun)) > 1e-11:
            raise RuntimeError("Middle contact fit failed")
        low, high = solution.x
    else:
        low = 0.0
        level = lambda x: G(1, 2, x) - g
        high = brentq(level, eta / 1.5, 12)
        u = brentq(level, 0, eta / 1.5)
        # Lower pair cannot improve on g, while the horizontal majorant
        # clears the left entrance obstacle strictly.
        lower_pair_margin = integral(lambda x: level(x) / q(x), u, high) - 2 * cost
        height = fit([low, high])[0]
        assert lower_pair_margin < 0 < height
        print(f"entrance height surplus: {height:.10f}")
        print(f"lower-pair area minus twice cost: {lower_pair_margin:.10f}")
    k = K(low)
    z = brentq(lambda x: K(x) - k, eta / 2.5, high)
    assert low < a < z < b < high
    scale = lambda x, y: integral(lambda t: 1 / q(t), x, y)
    p = scale(low, z) / scale(low, high) if low > 0 else 1.0

    # At z switch from 3 to 2, then leave 2 at low or high. The high exit
    # needs one additional pair switch at a before the cycle returns to z.
    def cycle_integral(middle, lower_exit, upper_exit):
        via_high = green(a, high, upper_exit, upward=False) + green(a, z, lower_exit)
        value = p * (green(z, high, middle) + via_high)
        if low > 0:
            value += (1 - p) * (green(z, low, middle) + green(low, z, lower_exit))
        return value

    duration = cycle_integral(None, None, None)
    net = cycle_integral(2, 3, 1) - cost * (2 + p)
    excess = net - g * duration
    outer_difference = cost - integral(lambda x: (G(1, 3, x) - g) / q(x), a, z)
    obstacle_excess = (
        integral(lambda x: (G(1, 2, x) - g + k) / q(x), z, high)
        - outer_difference - 2 * cost
    )
    assert abs(duration - p * scale(a, high)) < 1e-9
    assert abs(excess - obstacle_excess) < 1e-9 and excess > 0
    # The reported table rounds the obstacle excess to six decimal places.
    if prescribed_outer_gain is None:
        assert abs(excess - 0.015323) < 0.5e-6

    for name, value in [
        ("dimensionless cost", cost), ("outer pair gain", g),
        ("outer low", a), ("outer high", b), ("middle low", low), ("middle high", high),
        ("entry z", z), ("upper exit probability", p),
        ("mean cycle time", duration), ("excess per cycle", excess),
        ("gain improvement", excess / duration), ("deviation gain", net / duration),
    ]:
        print(f"{name}: {value:.10f}")


if __name__ == "__main__":
    print("Interior contact (b=0.060)")
    calculate(124 / 45, 1.5, cost=-(3.1 / 0.04) * log1p(-0.0018))
    print("\nEntrance contact (F=2, eta=2.1)")
    calculate(2, 2.1, prescribed_outer_gain=2.21)
