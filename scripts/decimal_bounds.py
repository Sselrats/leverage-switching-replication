"""Outward decimal representations of rational interval endpoints."""
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from fractions import Fraction


def decimal_bounds(lower: Fraction, upper: Fraction, digits: int) -> tuple[str, str]:
    """Enclose the exact interval using at most `digits` significant digits per endpoint."""
    floor = Context(prec=digits, rounding=ROUND_FLOOR)
    ceil = Context(prec=digits, rounding=ROUND_CEILING)
    left = floor.divide(Decimal(lower.numerator), Decimal(lower.denominator))
    right = ceil.divide(Decimal(upper.numerator), Decimal(upper.denominator))
    assert Fraction(left) <= lower <= upper <= Fraction(right)
    return str(left), str(right)
