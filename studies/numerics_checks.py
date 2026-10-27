"""Numerical checks of the QFI pipeline on real sweep points.

(a) finite-difference step: the QFI at steps 1e-4, 1e-5 (the package default) and 1e-6
    must sit on a plateau;
(b) exact derivative (sensitivity equation) vs finite difference at the same points;
(c) extended cutoff ladder (up to 560) at the hardest points, 10 Hz of the phase-noise
    sweep, with the exact derivative;
(d) sensitivity to the eigenvalue cutoff in the QFI sum (1e-6, 1e-8 default, 1e-10).
Writes results/numerics/summary.md.
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd

import backaction_qfi as b
import backaction_qfi.optimised as opt
from backaction_qfi.frame import qfi_from_rho_drho, unsheared_channel, rp_qfi_unsheared_exact
from quantum_sensing.sld import calculate_qfi

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "numerics"
OUT.mkdir(parents=True, exist_ok=True)
lines = ["# Numerical checks on sweep points", ""]


def probe(fam, N, slice_=None):
    if fam.startswith("opt_"):
        cfg, eta, pn = slice_
        return opt.load_probe(fam[4:], cfg, eta, pn, 5.0, N_basis=N)[0]
    return b.make_probe(fam, N, 5.0)


# ---- (a)+(b): step size and exact derivative, on points of the channel-loss + phase-noise
#               and three-loss sweeps
POINTS = [  # (label, family, slice, N, channel kwargs)
    ("channel loss + phase noise, 30 Hz", "squeezed", None, 220, dict(kappa_ba=0.949, eta_ch=0.95, pn_in=0.2)),
    ("channel loss + phase noise, 30 Hz", "opt_fock_sup", ("loss_ch", 0.95, 0.2), 220, dict(kappa_ba=0.949, eta_ch=0.95, pn_in=0.2)),
    ("channel loss + phase noise, 10 Hz", "fock", None, 220, dict(kappa_ba=9.53, eta_ch=0.95, pn_in=0.2)),
    ("three losses, 10 Hz", "squeezed", None, 220, dict(kappa_ba=9.53, eta_in=0.95, eta_ch=0.95, eta_out=0.95)),
    ("three losses, 30 Hz", "cat", None, 220, dict(kappa_ba=0.949, eta_in=0.95, eta_ch=0.95, eta_out=0.95)),
]
lines += ["## (a) finite-difference step and (b) exact derivative", "",
          "| point | probe | h = 1e-4 | h = 1e-5 (default) | h = 1e-6 | exact | (h=1e-5 − exact)/exact | time FD / exact |",
          "|---|---|---|---|---|---|---|---|"]
for label, fam, sl, N, kw in POINTS:
    psi = probe(fam, N, sl)
    rho1 = unsheared_channel(rho=psi, N_basis=N, eta_in=kw.get("eta_in", 1.0), pn_in=kw.get("pn_in", 0.0), kappa_ba=0.0)
    kw2 = {k: v for k, v in kw.items() if k not in ("eta_in", "pn_in")}
    vals = {}
    t0 = time.time()
    for h in (1e-4, 1e-5, 1e-6):
        vals[h] = float(calculate_qfi(unsheared_channel, param_value=0.0, param_type="epsilon_a", prec=h, rho=rho1, N_basis=N, **kw2))
    tfd = (time.time() - t0) / 3
    t0 = time.time(); ex = rp_qfi_unsheared_exact(psi, N_basis=N, **kw); tex = time.time() - t0
    lines.append(f"| {label} | {fam} | {vals[1e-4]:.8f} | {vals[1e-5]:.8f} | {vals[1e-6]:.8f} | {ex:.8f} | {(vals[1e-5]-ex)/ex:+.1e} | {tfd:.0f} s / {tex:.0f} s |")
    print(lines[-1], flush=True)

# ---- (c): extended ladder at the hardest points (phase-noise sweep, 10 Hz), exact derivative
lines += ["", "## (c) extended cutoff ladder at 10 Hz of the phase-noise sweep (exact derivative)", "",
          "| probe | N=160 | N=220 (sweep) | N=300 | N=400 | N=560 | (560 − 220)/220 |", "|---|---|---|---|---|---|---|"]
kw = dict(kappa_ba=9.53, eta_ch=0.95, pn_in=0.2)
for fam, sl in [("squeezed", None), ("fock", None), ("opt_fock_sup", ("loss_ch", 0.95, 0.2))]:
    row = {}
    for N in (160, 220, 300, 400, 560):
        t0 = time.time(); row[N] = rp_qfi_unsheared_exact(probe(fam, N, sl), N_basis=N, **kw)
        print(f"  ladder {fam} N={N}: {row[N]:.8f} ({time.time()-t0:.0f}s)", flush=True)
    lines.append(f"| {fam} | {row[160]:.6f} | {row[220]:.6f} | {row[300]:.6f} | {row[400]:.6f} | {row[560]:.6f} | {(row[560]-row[220])/row[220]:+.1e} |")

# ---- (d): eigenvalue cutoff in the QFI sum
lines += ["", "## (d) eigenvalue cutoff in the QFI sum", "",
          "| point | probe | eig_tol = 1e-6 | 1e-8 (default) | 1e-10 | (1e-10 − 1e-8)/1e-8 |", "|---|---|---|---|---|---|"]
for label, fam, sl, N, kw in POINTS[:3]:
    psi = probe(fam, N, sl)
    vals = {tol: rp_qfi_unsheared_exact(psi, N_basis=N, eig_tol=tol, **kw) for tol in (1e-6, 1e-8, 1e-10)}
    lines.append(f"| {label} | {fam} | {vals[1e-6]:.8f} | {vals[1e-8]:.8f} | {vals[1e-10]:.8f} | {(vals[1e-10]-vals[1e-8])/vals[1e-8]:+.1e} |")
    print(lines[-1], flush=True)

(OUT / "summary.md").write_text("\n".join(lines) + "\n")
print("summary:", OUT / "summary.md")
