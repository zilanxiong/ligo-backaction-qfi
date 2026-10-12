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
| `backaction_qfi/` | the study layer, built *on top of* the package without editing it: `solver` (tight ODE tolerances injected from outside, `rp_channel`, `rp_qfi`), `probes` (coherent, squeezed, cat, Fock at a given cutoff), `convergence` (Fock-cutoff ladder, truncation alarm) |
| `pyproject.toml` | packaging for both, so `pip install -e .` installs them |
| `environment.yml` | conda environment that reproduces this project's setup |
| `requirements.txt` | the pip side of the same setup; installs the package above in editable mode |
| `tests/` | smoke tests for the environment, and the Fock-cutoff convergence suite |

## Set up

```bash
conda env create -f environment.yml
conda activate ligo-backaction-qfi
pytest            # ~10 s; passing means the package is installed and the cutoff checks hold
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

## Fock cutoff

Radiation pressure shears the state, so the package default `N_basis=20`,
adequate for the *input* probe, can be badly wrong for the *output*. Every
back-action number should come from `backaction_qfi.converged_qfi`, which
climbs a ladder of cutoffs until the QFI stops moving and reports whether it
did:

```python
import backaction_qfi as b
res = b.converged_qfi(lambda N: b.squeezed(N, 2.0), kappa_ba=1.0, eta_out=0.9)
res.value, res.cutoff, res.converged
```

`tests/test_convergence.py` pins this down, including the regression case
where 20 levels is more than 10 % off.
