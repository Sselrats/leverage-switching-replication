# Optimal Switching among Leverage Levels: From Pairs to Menus

Numerical certificates for the v6.3 manuscript by Mugeon Kim. Section S6 of the mathematical supplement reports the calculations reproduced here.

Clone this repository at tag `v6.3` or extract the tagged source archive before running any program. An eight-file collection of
reference JSON records is in `reference/`. The Python programs in this archive generate the records from the stated
model inputs.

With Python 3.13, install the numerical dependencies in a virtual environment.
The replay driver and the checksum check use the Python standard library.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python replay.py check
```

Run any command from the extracted directory:

| Command | Result |
| --- | --- |
| `python replay.py menus` | Exact pair costs and increasing-path enclosures for three finite menus. |
| `python replay.py cycle` | Interval check of the example with an active cycle and inactive pairs. |
| `python replay.py gains` | Four-threshold gain, outer-pair comparison, deletion decision, premium-contact bracket and slope difference, at F = 124/45. |
| `python replay.py adjacent` | Replays the full-menu record first, then certifies both adjacent intervals and the gain difference at 18 bp. |
| `python replay.py transition` | Four interval computations for Proposition 4.3, from positive costs close to zero through 0.085 at F = 2. Allow several minutes. |
| `python replay.py excursions` | Green-integral check of the improving cycle for both an interior contact and an entrance contact. |
| `python replay.py points` | Rounded deletion table, cycle occupations and mean leverages at the premium boundary. |

The commands from `menus` through `transition` recompute the outputs, then compare **all** JSON
fields with the archived results except for the top-level Python version and
elapsed or runtime seconds and the two individual runtime fields `cases[*].seconds`
in the adjacent-interval record. The timing fields measure wall-clock time. The comparison reads gzip files as JSON. A difference
stops the run and leaves both files available for inspection. With different
libraries or numerical proposals, an enclosure may be valid yet differ from
the archived one. Run Python with assertions enabled.

The packaged checksums verify file identity. The mathematical acceptance tests
run inside the original programs. The JSON objects in `gains.json` include
`cases.lower`, `cases.nominal`, `nominal_exact_deletion_value` and
`premium_deletion_boundary`. The first-transition records, in execution order,
are `contact.json`, `small_cost.json.gz`, `smoothing.json.gz`, and
`transition.json.gz`. For the zero-fee menu {1, 2, 3} with F = 2 and eta = 3/2, an adjacent
policy attains g123 = gN on the initial interval (0, lambda_nn]. At lambda_nn
a four-threshold policy that uses the direct 1-to-3 switch in each cycle
also attains gN. On (lambda_nn, 0.085] the calibrated four-threshold
policy attains g123 > gN, so every adjacent policy is strictly suboptimal.
The band-touch bracket and transition bracket are recorded separately. The rounded examples use point arithmetic, with separate interval
checks for the two deletion signs and the transverse slope.

The replay stores generated files under `artifacts/` and
`scripts/paper_v11_research/records/`. The reference directory is read-only
from the perspective of the replay driver.
