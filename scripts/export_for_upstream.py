"""Assemble for_upstream/: this study laid out for the upstream quantum-sensing repository.

Upstream already holds the package (python/quantum_sensing) and the optimised-state
data, so neither is copied.  Everything else lands in one folder, backaction/, plus an
environment.yml for the repository root.  Run from this repository's root:

    python scripts/export_for_upstream.py

then copy for_upstream/* into the root of a clone of the upstream repository.
"""

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "for_upstream"
B = OUT / "backaction"

shutil.rmtree(OUT, ignore_errors=True)
B.mkdir(parents=True)

for d in ("backaction_qfi", "studies", "tests", "results", "figures"):
    shutil.copytree(ROOT / d, B / d, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))

(B / "pyproject.toml").write_text('''[build-system]
requires = ["setuptools>=64"]
build-backend = "setuptools.build_meta"

[project]
name = "backaction-qfi"
version = "0.1.0"
description = "QFI of optical probe states under radiation-pressure back-action (study layer over quantum_sensing)"
requires-python = ">=3.10"
dependencies = ["quantum-sensing", "matplotlib>=3.7", "pytest>=7.0"]

[tool.setuptools.packages.find]
include = ["backaction_qfi*"]

[tool.pytest.ini_options]
testpaths = ["tests"]
''')

(B / "requirements.txt").write_text('''# Run from the repository ROOT:   pip install -r backaction/requirements.txt
-e ./python          # the group's quantum_sensing package
-e ./backaction      # this study layer
''')

(OUT / "environment.yml").write_text('''# Reproducible environment for the back-action study.  From the repository root:
#   conda env create -f environment.yml
#   conda activate quantum-sensing-backaction
#   export QUANTUM_SENSING_DATA_DIR=<directory holding states.h5 and master_results.csv>
#   pytest backaction            # 98 tests
name: quantum-sensing-backaction
channels:
  - conda-forge
dependencies:
  - python=3.12
  - numpy>=2.0,<2.3
  - scipy>=1.13
  - pandas>=2.0
  - h5py
  - matplotlib>=3.7
  - pytest>=7.0
  - pip
  - pip:
      - -e ./python
      - -e ./backaction
''')

readme = (ROOT / "README.md").read_text()
(B / "README.md").write_text('''# backaction/ — QFI under radiation-pressure back-action

Study layer over the group's `quantum_sensing` package (installed from `../python`,
never modified here).  Everything below is run from the **repository root**.

## Set up

```bash
conda env create -f environment.yml            # environment.yml sits at the repo root
conda activate quantum-sensing-backaction
export QUANTUM_SENSING_DATA_DIR=<dir with states.h5 and master_results.csv>
pytest backaction                              # 98 tests, ~20 s
```

Without conda: any Python >= 3.10 environment, then `pip install -r backaction/requirements.txt`.

## Layout

| path | what it is |
|---|---|
| `backaction_qfi/` | `solver` (tight ODE tolerances injected without editing the package; `rp_channel`, `rp_qfi`), `probes` (coherent, squeezed, cat, squeezed cat, Fock at fixed photon number), `convergence` (Fock-cutoff ladder, truncation alarm), `gaussian` (covariance-matrix track, SQL / variational / FD-squeezing curves), `optimised` (loads the group's optimised probes, rotated onto the GW signal `epsilon_a`) |
| `studies/backaction_only.py` | first sweep: QFI vs frequency (10 Hz – 1 kHz) and vs photon number, back-action only |
| `results/backaction_only/` | its CSVs, `summary.md` and run log |
| `figures/backaction_only.png` | its figure |
| `tests/` | convergence, Gaussian benchmarks, optimised-state orientation, environment smoke tests |

## Reproduce the sweep

```bash
python backaction/studies/backaction_only.py --figure-only   # figure from the committed CSVs, seconds
python backaction/studies/backaction_only.py --n-max 400     # full recomputation, ~15 min
```

## Notes for the package author

* The optimised states reproduce their stored QFI for `epsilon_p`, not for the GW signal
  `epsilon_a`; the loader rotates them by a quarter turn (`c_n -> (-i)^n c_n`).
* Squeezed-family states reconstruct with a mean photon number ~1.7 % above the stored
  `N_actual` (cutoff independent) while their QFI reproduces; Fock superpositions match
  exactly.  See `tests/test_optimised_states.py`.
* `dynamics.py` uses QuTiP's default solver tolerances; the finite-difference QFI step is
  1e-5, close to the solver error (~4e-6).  `backaction_qfi.solver` tightens them from outside.
''')

print("bundle at", OUT)
for p in sorted(OUT.rglob("*")):
    if p.is_file() and "results" not in p.parts and "figures" not in p.parts:
        print("  ", p.relative_to(OUT))
