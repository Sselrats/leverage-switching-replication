"""Print the two tangent-subtracted branches of Figure 2 as CSV."""
import csv
import sys

from scipy.optimize import brentq

from rounded_examples import Full, contact_residual


def main():
    eta_star = brentq(contact_residual, 1.39631804, 1.39631806, xtol=5e-15)
    at_contact = Full(eta_star).sol[13]
    writer = csv.writer(sys.stdout)
    writer.writerow(["eta", "outer_pair", "four_threshold"])
    for step in range(-20, 21):
        delta = step * 0.00005
        eta = eta_star + delta
        model = Full(eta)
        # Subtract the contact gain first to match the manuscript's plotted coordinates.
        tangent_increment = at_contact["ell"] * delta
        pair = 10000 * (model.sol[13]["g_dim"] - at_contact["g_dim"] - tangent_increment)
        four = 10000 * (model.solve()["gain_dim"] - at_contact["g_dim"] - tangent_increment)
        writer.writerow([f"{eta:.12f}", f"{pair:.10f}", f"{four:.10f}"])


if __name__ == "__main__":
    main()
