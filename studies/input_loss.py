"""Input loss + back-action, nbar = 5, 10 Hz - 1 kHz, at the paper's headline transmission 0.95.

Baseline check against the paper (Maliakal et al.): input loss acts before the
sensing unitary and commutes with nothing that follows in a way that could
change the QFI, so the QFI must be frequency independent and equal to the
no-back-action value — the stored QFI of the optimised states (loss_in slice)
and the closed forms for squeezed vacuum, 4/(eta e^{-2r} + 1 - eta), and
coherent light, 4 eta.

Fock track: the package's density-matrix channel with a convergence-checked
cutoff where that fits the budget; beyond it (large kappa) the exact
branch method (backaction_qfi.branches), which runs the noise stage through
the package and the sensing unitary at a cutoff of ~1000-1600.

Usage:  python studies/input_loss.py [--n-max 400] [--etas 0.95 ...] [--quick] [--figure-only]
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

import backaction_qfi as b
import backaction_qfi.gaussian as g
import backaction_qfi.optimised as opt
from backaction_qfi.branches import _unitaries
from backaction_qfi.convergence import default_ladder
from quantum_sensing import rp_kappa

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "input_loss"
FIG = ROOT / "figures"

NBAR = 5.0
ETAS = [0.95]          # the paper's headline transmission (5 % loss); others via --etas
FREQS_HZ = np.logspace(1, 3, 25)
N_REF = 400
BRANCH_LADDER = [800, 1120, 1600, 2240]
PARAMETRIC = ["squeezed", "coherent", "cat", "squeezed_cat", "fock"]
OPTIMISED = ["opt_fock_sup", "opt_sqz_vac", "opt_sqz_coh"]
FAMILIES = PARAMETRIC + OPTIMISED
LEGEND = {"squeezed": "Squeezed vac. (p)", "coherent": "Coherent", "cat": "Even cat",
          "squeezed_cat": "Squeezed cat", "fock": "Fock",
          "opt_fock_sup": "Optimized Fock sup.", "opt_sqz_vac": "Optimized sqz-vac sup.",
          "opt_sqz_coh": "Optimized sqz-coh sup."}
GAUSSIAN = {"squeezed": lambda n: g.squeezed_at_fixed_nbar(n), "coherent": lambda n: g.coherent_gaussian(np.sqrt(n))}


def _available_etas(state_type):
    df = pd.read_csv(opt.DATA_DIR / "master_results.csv", dtype={"valid": str})
    df = df[(df.valid.fillna("").str.lower() != "false") & (df.state_type == state_type)
            & (df.loss_config == "loss_in") & (df.pn == 0.0) & (df.N_target == NBAR)]
    return sorted(df.eta.unique())


def opt_slice(family, eta):
    """(eta actually used, note): the matching loss_in slice, else the nearest available."""
    avail = _available_etas(family[4:])
    if eta in avail:
        return eta, ""
    near = min(avail, key=lambda e: abs(e - eta))
    return near, f"no loss_in entry at eta={eta}; nearest slice eta={near} used"


def probe(family, N, eta):
    if family.startswith("opt_"):
        e, _ = opt_slice(family, eta)
        psi, _ = opt.load_probe(family[4:], "loss_in", e, 0.0, NBAR, N_basis=N)
        return psi
    return b.make_probe(family, N, NBAR)


def reference(family, eta):
    """No-back-action reference the sweep must reproduce: paper/closed form where one exists."""
    r = np.arcsinh(np.sqrt(NBAR))
    if family == "squeezed":
        return 4 / (eta * np.exp(-2 * r) + 1 - eta), "closed form (paper Eq. 17)"
    if family == "coherent":
        # input loss leaves a coherent state pure (amplitude -> sqrt(eta) alpha); the signal
        # is applied afterwards, so F = 8 Var(x) = 4 for any eta_in (4 eta would be output loss)
        return 4.0, "closed form: 4, independent of input loss"
    if family.startswith("opt_"):
        e, _ = opt_slice(family, eta)
        return opt.best_entry(family[4:], "loss_in", e, 0.0, NBAR)["qfi"], f"stored QFI (paper data, eta={e})"
    # no closed form: exact branch method at kappa = 0 (no back-action)
    psi = probe(family, 200, eta)
    return b.rp_qfi_input_noise_branches(psi, kappa_ba=0.0, eta_in=eta, N_unitary=400), "no-back-action Fock value"


def sweep(etas, n_max, quick, checkpoint=None):
    freqs = FREQS_HZ[::4] if quick else FREQS_HZ
    kappas = rp_kappa(2 * np.pi * freqs)
    ladder = [N for N in default_ladder(40, n_max) if N <= n_max]
    rows = []
    for eta in etas:
        unit_cache = {}
        for fam in FAMILIES:
            ref, ref_src = reference(fam, eta)
            e_used, e_note = opt_slice(fam, eta) if fam.startswith("opt_") else (np.nan, "")
            nbar = b.mean_photon_number(probe(fam, N_REF, eta))
            print(f"[eta_in={eta}] [{fam}] <n>={nbar:.4f}  reference {ref:.4f} ({ref_src}) {e_note}", flush=True)
            for f, k in zip(freqs, kappas):
                t0 = time.time()
                gauss = g.gaussian_qfi(GAUSSIAN[fam](NBAR), kappa_ba=k, eta_in=eta) if fam in GAUSSIAN else np.nan
                predicted = b.suggested_cutoff(NBAR, k, safety=2.0)
                if predicted <= n_max:
                    start = max(40, int(0.6 * predicted))
                    rungs = [N for N in ladder if N >= start] or [n_max]
                    if len(rungs) < 2:
                        rungs = [start, n_max] if start < n_max else [int(0.75 * n_max), n_max]
                    res = b.converge(lambda N: b.rp_qfi(probe(fam, N, eta), kappa_ba=k, eta_in=eta, N_basis=N),
                                     cutoffs=rungs, rtol=1e-3)
                    track = "mesolve"
                else:
                    def value(N):
                        key = (k, N)
                        if key not in unit_cache:
                            unit_cache[key] = _unitaries(N, k, 1e-5)
                        return b.rp_qfi_input_noise_branches(probe(fam, 200, eta), kappa_ba=k, eta_in=eta,
                                                             N_unitary=N, unitaries=unit_cache[key])
                    res = b.converge(value, cutoffs=BRANCH_LADDER, rtol=1e-3)
                    track = "branches"
                rows.append(dict(eta_in=eta, family=fam, freq_hz=f, kappa_ba=k, nbar=nbar, qfi_fock=res.value,
                                 cutoff=res.cutoff, converged=res.converged, track=track, qfi_reference=ref,
                                 reference_source=ref_src, rel_err_vs_reference=(res.value - ref) / ref,
                                 qfi_gaussian=gauss, opt_eta_used=e_used, note=e_note, seconds=time.time() - t0))
                flag = "" if res.converged else "  NOT CONVERGED"
                print(f"   f={f:7.1f} Hz kappa={k:8.4f}  F={res.value:9.4f} @N={res.cutoff:<4d} [{track}] "
                      f"(ref {ref:.4f}, rel {rows[-1]['rel_err_vs_reference']:+.1e})  {rows[-1]['seconds']:.0f}s{flag}", flush=True)
            if checkpoint is not None:
                pd.DataFrame(rows).to_csv(checkpoint, index=False)
    return pd.DataFrame(rows)


def _kappa_axis(ax):
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    from quantum_sensing import rp_power_ratio
    gamma = 2 * np.pi * 450.0
    X = rp_power_ratio(gamma=gamma)

    def k_to_f(k):
        W2 = 0.5 * (-gamma**2 + np.sqrt(gamma**4 + 8 * X * gamma**4 / np.asarray(k, float)))
        return np.sqrt(W2) / (2 * np.pi)
    fmin, fmax = ax.get_xlim()
    ticks = [(k_to_f(10.0 ** e), e) for e in range(2, -6, -1) if fmin <= k_to_f(10.0 ** e) <= fmax]
    top = ax.twiny(); top.set_xscale("log"); top.set_xlim(fmin, fmax)
    top.xaxis.set_major_locator(FixedLocator([f for f, _ in ticks]))
    top.xaxis.set_major_formatter(FixedFormatter([f"$10^{{{e}}}$" for _, e in ticks]))
    top.xaxis.set_minor_locator(NullLocator()); top.set_xlabel(r"$\kappa_{ba}(\Omega)$"); top.grid(False)


def figure(df, eta, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcdefaults()
    plt.rcParams.update({"figure.figsize": (8, 5.3), "figure.dpi": 150, "savefig.dpi": 150,
                         "axes.grid": True, "grid.alpha": 0.3, "legend.framealpha": 0.9,
                         "lines.linewidth": 1.5, "lines.markersize": 4})
    fig, ax = plt.subplots()
    d_eta = df[df.eta_in == eta]
    for fam in [f for f in FAMILIES if f in set(d_eta.family)]:
        d = d_eta[d_eta.family == fam].sort_values("freq_hz")
        ok, bad = d[d.converged], d[~d.converged]
        label = LEGEND[fam]
        if fam.startswith("opt_") and d.note.iloc[0]:
            label += f" (η={d.opt_eta_used.iloc[0]:g} slice)"
        line, = ax.plot(ok.freq_hz, ok.qfi_fock, "-o", label=label)
        if len(bad):
            ax.plot(bad.freq_hz, bad.qfi_fock, "o", mfc="white", color=line.get_color())
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"Frequency  $\Omega/2\pi$  [Hz]"); ax.set_ylabel(r"QFI for $\varepsilon_a$")
    ax.set_title(f"input loss ($\\eta_{{in}}$ = {eta:g}) → backaction + displacement, ⟨n⟩ = 5")
    ax.grid(True, which="both", alpha=0.3)
    _kappa_axis(ax)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    print("figure:", path)


def write_summary(df, path):
    lines = ["# Input loss + back-action, ⟨n⟩ = 5, 10 Hz – 1 kHz", "",
             "Input loss acts before the sensing unitary, so the QFI must be frequency independent and",
             "equal to its no-back-action value. Reference: the stored QFI of the paper's optimised states",
             "(loss_in slice), closed forms for squeezed vacuum (paper Eq. 17) and coherent light (4), and the exact",
             "no-back-action Fock value for the rest.", ""]
    for eta in sorted(df.eta_in.unique()):
        d_eta = df[df.eta_in == eta]
        lines += [f"## η_in = {eta:g}", "",
                  "| probe | reference | source | Fock track over f (min – max) | worst rel. dev. | converged | tracks | largest cutoff |",
                  "|---|---|---|---|---|---|---|---|"]
        for fam in [f for f in FAMILIES if f in set(d_eta.family)]:
            d = d_eta[d_eta.family == fam]; ok = d[d.converged]
            tr = ", ".join(f"{t}×{n}" for t, n in d.track.value_counts().items())
            lines.append(f"| {LEGEND[fam]} | {d.qfi_reference.iloc[0]:.4f} | {d.reference_source.iloc[0]} | "
                         f"{ok.qfi_fock.min():.4f} – {ok.qfi_fock.max():.4f} | {ok.rel_err_vs_reference.abs().max():.1e} | "
                         f"{len(ok)}/{len(d)} | {tr} | {int(d.cutoff.max())} |")
        notes = d_eta[d_eta.note != ""].drop_duplicates("family")
        for _, r in notes.iterrows():
            lines.append(f"\n*{LEGEND[r.family]}: {r.note}.*")
        lines.append("")
    path.write_text("\n".join(lines) + "\n")
    print("summary:", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-max", type=int, default=400)
    ap.add_argument("--etas", nargs="*", type=float, default=ETAS)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--figure-only", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)
    csv = OUT / "qfi_vs_frequency.csv"
    if args.figure_only:
        df = pd.read_csv(csv); df["note"] = df["note"].fillna("")
    else:
        df = sweep(args.etas, args.n_max, args.quick, checkpoint=OUT / "checkpoint.csv")
        df.to_csv(csv, index=False)
        (OUT / "checkpoint.csv").unlink(missing_ok=True)
    write_summary(df, OUT / "summary.md")
    for eta in sorted(df.eta_in.unique()):
        figure(df, eta, FIG / (f"input_loss_eta{eta:g}".replace(".", "p") + ".png"))


if __name__ == "__main__":
    sys.exit(main())
