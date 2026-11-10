"""Rounded middle-level deletion table and cycle means at the premium boundary."""
from fractions import Fraction
from pathlib import Path
import sys

import mpmath as mp
from scipy.optimize import brentq

SCRIPTS = Path(__file__).resolve().parent / "scripts"
sys.path.insert(0, str(SCRIPTS / "cir_fullmenu_endpoints_20260927"))
sys.path.insert(0, str(SCRIPTS / "paper_v11_research/round12"))
from explore import Full, LAM
from endpoint_scan import NumericCIR


def deletion_row(eta):
    model = Full(float(eta))
    cir = NumericCIR(Fraction(124, 45), eta)
    seed = mp.mpf(model.sol[13]["g_dim"])
    cost = -mp.mpf("77.5") * mp.log(mp.mpf(".9982"))
    g = mp.findroot(lambda gain: cir.pair(13, gain)["H"] / 2 - cost,
                    (seed - mp.mpf("1e-7"), seed + mp.mpf("1e-7")))
    chord = cir.least(g)
    return (model.sol[12]["g_dim"] - float(g), float(chord["roots"]["gamma"]),
            float(chord["roots"]["beta"]), float(chord["k"]), float(-min(chord["margins"].values())))


def contact_residual(eta):
    model = Full(eta)
    outer, upper = model.sol[13], model.sol[23]
    k = upper["g_dim"] - model.B2
    a, gamma = outer["u_dim"], upper["u_dim"]
    _, beta = model.pairs[12].band(outer["g_dim"] - k)
    return model.intL(12, a, beta, outer["g_dim"] - k) + model.intK(a, gamma, k) - LAM


if __name__ == "__main__":
    expected = {"0.95": (-.009345, .059770, .998177, .001163, -.053247),
                "1.5": (-.049329, .242018, 1.565393, .048868, .015323)}
    for eta, target in expected.items():
        row = deletion_row(eta)
        print("eta =", eta, "  D, gamma, beta, k, margin =", row)
        assert tuple(round(x, 6) for x in row) == target

    eta = brentq(contact_residual, 1.39631804, 1.39631806, xtol=5e-15)
    model = Full(eta)
    branch = model.solve()
    gamma, a, zeta, beta = (branch["roots"][x] for x in ("gamma", "a", "zeta", "beta"))
    E, A = model.E, model.A
    p = (E(zeta) - E(gamma)) / (E(beta) - E(gamma))
    T = E(beta) - E(a)
    T3 = A(zeta) - A(a) + (1 - p) / p * (A(zeta) - A(gamma))
    T2 = A(beta) - A(a) - T3
    T1 = T - T2 - T3
    ell = (T1 + 2 * T2 + 3 * T3) / T
    outer = model.sol[13]["ell"]
    print("eta* =", eta, "  p =", p, "  T1, T2, T3 =", (T1, T2, T3))
    print("outer and four-threshold mean leverages =", (outer, ell))
    assert (round(outer, 8), round(ell, 8)) == (1.65926812, 1.69725183)
    assert abs(ell - branch["ell"]) < 1e-10
