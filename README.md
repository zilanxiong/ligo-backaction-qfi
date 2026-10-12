# ligo-backaction-qfi

Quantum Fisher information (QFI) of optical probe states under
radiation-pressure back-action, built on the group's `quantum_sensing`
package.

This is a fresh start. The earlier, exploratory work lives in
[LIGO-backaction-simulation](https://github.com/zilanxiong/LIGO-backaction-simulation),
kept as a read-only reference.

## Layout

| path | what it is |
|---|---|
| `quantum_sensing/` | the group's single-mode sensing package, vendored unmodified from `shruti-jm/quantum-sensing` (`python/quantum_sensing/`). Not edited here |
| `pyproject.toml` | packaging for `quantum_sensing`, so `pip install -e .` installs it |
| `environment.yml` | conda environment that reproduces this project's setup |
| `requirements.txt` | the pip side of the same setup; installs the package above in editable mode |
| `tests/` | smoke tests that the environment is set up correctly |

## Set up

```bash
conda env create -f environment.yml
conda activate ligo-backaction-qfi
pytest            # should pass: the package is installed and importable
```

Without conda: create any Python 3.10+ virtual environment and run
`pip install -r requirements.txt`.

## Using the package

```python
import quantum_sensing as qs
```

`quantum_sensing/` is the upstream code as received. To
update it, replace that folder with the newer `python/quantum_sensing/` from
upstream; the editable install picks the change up without reinstalling.
