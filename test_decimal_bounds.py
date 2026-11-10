import json
import unittest
from fractions import Fraction
from pathlib import Path

from scripts.decimal_bounds import decimal_bounds


REFERENCE = Path(__file__).resolve().parent / "reference"


def intervals(record, path=""):
    if isinstance(record, dict):
        if "lo" in record and "hi" in record:
            yield path, record
        for key, value in record.items():
            yield from intervals(value, f"{path}.{key}" if path else key)
    elif isinstance(record, list):
        for index, value in enumerate(record):
            yield from intervals(value, f"{path}[{index}]")


class DecimalBoundsTests(unittest.TestCase):
    def test_two_significant_digits(self):
        tiny = Fraction(1, 3 * 10**100)
        cases = [
            (Fraction(1, 3), Fraction(2, 3), "0.33", "0.67"),
            (Fraction(-100, 3), Fraction(-100, 3), "-34", "-33"),
            (Fraction(-2, 3), Fraction(1, 3), "-0.67", "0.34"),
            (Fraction(0), Fraction(0), "0", "0"),
            (Fraction(-1, 4), Fraction(1, 8), "-0.25", "0.13"),
            (Fraction(1, 8), Fraction(1, 4), "0.12", "0.25"),
            (-tiny, tiny, "-3.4e-101", "3.4e-101"),
            (tiny, tiny, "3.3e-101", "3.4e-101"),
        ]
        for low, high, expected_low, expected_high in cases:
            with self.subTest(low=low, high=high):
                displayed = decimal_bounds(low, high, 2)
                self.assertTrue(all(isinstance(endpoint, str) for endpoint in displayed))
                self.assertEqual(tuple(map(Fraction, displayed)),
                                 (Fraction(expected_low), Fraction(expected_high)))

    def test_archived_decimal_bounds_enclose_the_exact_bounds(self):
        for filename in ("contact.json", "cycle.json"):
            record = json.loads((REFERENCE / filename).read_text(encoding="utf-8"))
            boxes = list(intervals(record))
            self.assertTrue(boxes, filename)
            for path, box in boxes:
                if filename == "contact.json":
                    displayed = box["decimal_lo"], box["decimal_hi"]
                else:
                    displayed = box["display"]
                low, high = map(Fraction, displayed)
                exact_low, exact_high = Fraction(box["lo"]), Fraction(box["hi"])
                self.assertLessEqual(low, exact_low, f"{filename} {path} lower")
                self.assertLessEqual(exact_low, exact_high, f"{filename} {path} exact")
                self.assertLessEqual(exact_high, high, f"{filename} {path} upper")


if __name__ == "__main__":
    unittest.main()
