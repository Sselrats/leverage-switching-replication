"""Replay the numerical certificates for Optimal Switching among Leveraged Products."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"
RESEARCH = SCRIPTS / "paper_v11_research"
RESULTS = ROOT / "artifacts"
REFERENCE = ROOT / "reference"


def run(script, *args):
    subprocess.run([sys.executable, str(script), *map(str, args)], check=True, cwd=ROOT)


def read_record(path, reference):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as file:
        record = json.load(file)
    for key in ("elapsed_seconds", "runtime_seconds", "python_version"):
        record.pop(key, None)
    if reference == "adjacent.json":
        for case in record["cases"]:
            case.pop("seconds")
    return record


def compare(actual, reference):
    if read_record(actual, reference) != read_record(REFERENCE / reference, reference):
        raise SystemExit(f"Replay differs from reference: {actual}. Inspect both JSON records.")
    print(f"Matches {reference}", flush=True)


def check():
    for line in (ROOT / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise SystemExit(f"Checksum mismatch: {name}")
    print("All archived files match their SHA-256 checksums.")


def menus():
    out = RESULTS / "menus.json"
    run(SCRIPTS / "menu_enclosures_20260929/verify_s5_exact.py", "--output", out)
    compare(out, "menus.json")


def cycle():
    out = RESULTS / "v40_rewrite_20260928/cycle_certificate.json"
    out.parent.mkdir(exist_ok=True)
    run(SCRIPTS / "v40_rewrite_20260928/check_cycle.py")
    compare(out, "cycle.json")


def gains():
    out = RESULTS / "cir_fullmenu_endpoints_20260927/full_menu_certificate.json"
    out.parent.mkdir(exist_ok=True)
    run(SCRIPTS / "cir_fullmenu_endpoints_20260927/certify.py")
    compare(out, "gains.json")


def adjacent():
    gains()
    out = RESULTS / "adjacent_interval_20260929/record.json"
    run(SCRIPTS / "adjacent_interval_20260929/certify.py")
    compare(out, "adjacent.json")


def transition():
    targets = [
        ("round11/validate_contacts.py", ["nn"], "nn_contact.json", "contact.json"),
        ("round10/certify_small_cost.py", [], "small_cost_test.json", "small_cost.json.gz"),
        ("round10/certify_quadratic_smoothing.py", [], "quadratic_smoothing.json", "smoothing.json.gz"),
        ("round11/nn_transition_cover.py", [], "nn_transition_cover.json.gz", "transition.json.gz"),
    ]
    for source, args, actual, reference in targets:
        run(RESEARCH / source, *args)
        compare(RESEARCH / "records" / actual, reference)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=("check", "menus", "cycle", "gains", "adjacent", "transition", "points", "excursions"))
    task = parser.parse_args().task
    if not __debug__:
        raise SystemExit("Enable Python assertions to run the certificates.")
    if task not in ("check", "points", "excursions"):
        RESULTS.mkdir(exist_ok=True)
    if task == "points":
        run(ROOT / "rounded_examples.py")
    elif task == "excursions":
        run(ROOT / "deviation_cycle.py")
    else:
        {"check": check, "menus": menus, "cycle": cycle, "gains": gains,
         "adjacent": adjacent, "transition": transition}[task]()
