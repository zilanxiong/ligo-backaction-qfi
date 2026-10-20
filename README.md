# ligo-backaction-qfi

Quantum Fisher information (QFI) of optical probe states under
radiation-pressure back-action, built on the group's `quantum_sensing`
package.

## Layout

| path | what it is |
|---|---|
| `quantum_sensing/` | the group's single-mode sensing package, vendored unmodified from `shruti-jm/quantum-sensing` (`python/quantum_sensing/`). Not edited here |
| `backaction_qfi/` | the study layer, built *on top of* the package without editing it: `solver` (tight ODE tolerances injected from outside, `rp_channel`, `rp_qfi`, pure-state track), `branches` (input noise then the lossless unitary, exact at large cutoffs), `frame` (the channel in the frame moving with the shear: exact at any κ with ~200 levels), `probes` (coherent, squeezed, cat, Fock at a given cutoff), `convergence` (Fock-cutoff ladder, truncation alarm), `gaussian` (covariance-matrix track: exact for Gaussian probes, SQL / variational / FD-squeezing curves, strain conversion), `optimised` (loads the group's optimised probes, rotated onto the GW signal) |
| `pyproject.toml` | packaging for both, so `pip install -e .` installs them |
| `environment.yml` | conda environment that reproduces this project's setup |
| `requirements.txt` | the pip side of the same setup; installs the package above in editable mode |
| `data/optimised_states/` | the group's previously optimised probes (no back-action), with the orientation note in its README |
| `scripts/export_for_upstream.py` | assembles `for_upstream/`, this study laid out for the group's repository (see *Uploading to the group's repository*) |
| `studies/` | the sweeps. `backaction_only.py`: QFI vs frequency with back-action only. `input_loss.py`: input loss (transmission 0.95) + back-action, checked against the paper's numbers. `channel_loss.py`: loss concurrent with the shear, the first frequency-dependent case; `--pn-in 0.2` adds the paper's 200 mrad preparation phase noise |
| `results/`, `figures/` | outputs of the studies, one folder per study, with a `summary.md` each |
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

## Studies

### Back-action only (`studies/backaction_only.py`)

No loss, no phase noise; ⟨n⟩ = 5; 25 frequencies from 10 Hz to 1 kHz with
κ(Ω) from the package's LIGO defaults; eight probes (squeezed vacuum, coherent,
even cat, squeezed cat, Fock, and the three optimised families). Every point is
computed through the package's density-matrix channel with a convergence-checked
cutoff and compared with the exact lossless value 8 Var(x).

```bash
python studies/backaction_only.py --n-max 400     # ~15 min; --quick for 7 frequencies
python studies/backaction_only.py --figure-only   # re-render from saved CSVs
```

Result: with back-action alone the QFI is frequency independent for every probe
(the shear generator x² commutes with the signal generator x). See
`results/backaction_only/summary.md` and `figures/backaction_only_qfi_vs_frequency.png`.


### Input loss + back-action (`studies/input_loss.py`)

Same probes and frequencies, with 5 % input loss (the paper's headline transmission 0.95) and
the optimised states taken from the paper's matching `loss_in` slice. Input loss acts before
the sensing unitary, so the QFI is frequency independent and must equal its no-back-action
value: the sweep reproduces the paper's stored QFI and closed forms to 1e-7 for the Fock
superposition, squeezed vacuum and coherent light. Beyond the density-matrix cutoff budget
the exact branch method (`backaction_qfi.branches`) carries the low-frequency points. See
`results/input_loss/summary.md`.

### Channel loss + back-action (`studies/channel_loss.py`)

5 % loss *during* the sensing stage, in the same master-equation solve as the shear and
the signal. The first frequency-dependent result: every probe loses most of its QFI
below ~30 Hz (κ ≳ 1). Primary track: the exact shear-frame solver
(`backaction_qfi.frame`), which keeps the state compact at any κ; cross-checked against
the package's own density-matrix channel at the low-κ points, and at κ → 0 against the
paper's stored QFI and the β² closed forms. See `results/channel_loss/summary.md`.

### Phase noise + channel loss + back-action (`studies/channel_loss.py --pn-in 0.2`)

The paper's headline operating point: 200 mrad of input phase noise and 5 % channel loss,
with the optimised states from that exact slice. At high frequency the paper's ordering is
reproduced (the three optimised families at 30.4-30.5 against squeezed vacuum at 18.2, the
2.2 dB non-Gaussian advantage). Below ~30 Hz the ordering changes: the Fock state degrades
least and is the best probe at 10 Hz. See `results/channel_loss_pn0p2/summary.md`.
