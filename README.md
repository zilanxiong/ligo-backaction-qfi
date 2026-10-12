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
| `backaction_qfi/` | the study layer, built *on top of* the package without editing it: `solver` (tight ODE tolerances injected from outside, `rp_channel`, `rp_qfi`), `probes` (coherent, squeezed, cat, Fock at a given cutoff), `convergence` (Fock-cutoff ladder, truncation alarm), `gaussian` (covariance-matrix track: exact for Gaussian probes, SQL / variational / FD-squeezing curves, strain conversion), `optimised` (loads the group's optimised probes, rotated onto the GW signal) |
| `pyproject.toml` | packaging for both, so `pip install -e .` installs them |
| `environment.yml` | conda environment that reproduces this project's setup |
| `requirements.txt` | the pip side of the same setup; installs the package above in editable mode |
| `data/optimised_states/` | the group's previously optimised probes (no back-action), with the orientation note in its README |
| `tests/` | smoke tests for the environment, the Fock-cutoff convergence suite, the Gaussian benchmark suite, and the optimised-state loader |

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

## Gaussian track

`backaction_qfi.gaussian` carries the same channel at the level of first and
second moments. It is exact for Gaussian probes (vacuum, coherent, squeezed)
and Gaussian noise (loss, including loss *during* the shear), and it is the
reference the Fock-space track is checked against. It also holds the textbook
free-mass curves (KLMTV): the standard quantum limit `(1+κ²)/(2κ)`, variational
readout `1/(2κ)`, and frequency-dependent squeezing `e^{-2r}(1+κ²)/(2κ)`, all
pinned down in `tests/test_gaussian.py`. Phase noise and dephasing are
non-Gaussian and live only in the Fock track.

```python
import backaction_qfi.gaussian as g
g.gaussian_qfi(g.squeezed_at_fixed_nbar(2.0), kappa_ba=1.0, eta_out=0.9)   # F for epsilon_a
g.homodyne_noise_spectrum(kappa, r=1.0, squeeze_angle=g.fd_squeeze_angle(kappa))
```

## Optimised probes

`backaction_qfi.optimised.load_probe` reconstructs the group's optimised
states from `data/optimised_states/`. They were optimised for the `epsilon_p`
signal; the loader rotates them by a quarter turn so they are oriented for the
GW signal `epsilon_a` (see `data/optimised_states/README.md`).

```python
from backaction_qfi.optimised import load_probe
psi, entry = load_probe("fock_sup", "loss_ch", eta=0.9, pn=0.0, N_target=5.0)
```
