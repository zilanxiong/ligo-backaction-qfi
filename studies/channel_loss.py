"""Channel loss (concurrent with the shear) + back-action, nbar = 5, 10 Hz - 1 kHz, transmission 0.95.

Loss acts *during* the sensing stage, in the same master-equation solve as the shear
and the signal, so it does not commute with the back-action: this is the first sweep
where the QFI depends on frequency.  At kappa -> 0 it must reduce to the paper's
numbers for loss during displacement: the stored QFI of the loss_ch optimised states,
and for squeezed vacuum / coherent light the loss-then-displacement closed forms scaled
by beta^2, beta = 2(1 - sqrt(eta))/|ln eta| (paper supplement Eq. 16).

Fock track: the exact shear-frame solver (backaction_qfi.frame), which needs only
~200 levels at any kappa; cross-checked against the package's own density-matrix
channel at the frequencies where that is cheap (small kappa, column qfi_mesolve).

Usage:  python studies/channel_loss.py [--n-max 400] [--etas 0.95 ...] [--quick] [--figure-only]
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
from backaction_qfi.convergence import default_ladder
from quantum_sensing import rp_kappa

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "channel_loss"
FIG = ROOT / "figures"

NBAR = 5.0
ETAS = [0.95]
FREQS_HZ = np.logspace(1, 3, 25)
N_REF = 400
FRAME_LADDER = [120, 160, 220, 300, 400]
PARAMETRIC = ["squeezed", "coherent", "cat", "squeezed_cat", "fock"]
OPTIMISED = ["opt_fock_sup", "opt_sqz_vac", "opt_sqz_coh"]
FAMILIES = PARAMETRIC + OPTIMISED
LEGEND = {"squeezed": "Squeezed vac. (p)", "coherent": "Coherent", "cat": "Even cat",
          "squeezed_cat": "Squeezed cat", "fock": "Fock",
          "opt_fock_sup": "Optimized Fock sup.", "opt_sqz_vac": "Optimized sqz-vac sup.",
          "opt_sqz_coh": "Optimized sqz-coh sup."}
GAUSSIAN = {"squeezed": lambda n: g.squeezed_at_fixed_nbar(n), "coherent": lambda n: g.coherent_gaussian(np.sqrt(n))}
SLICE = "loss_ch"


def _available_etas(state_type):
    df = pd.read_csv(opt.DATA_DIR / "master_results.csv", dtype={"valid": str})
    df = df[(df.valid.fillna("").str.lower() != "false") & (df.state_type == state_type)
            & (df.loss_config == SLICE) & (df.pn == 0.0) & (df.N_target == NBAR)]
    return sorted(df.eta.unique())


def opt_slice(family, eta):
    avail = _available_etas(family[4:])
    if eta in avail:
        return eta, ""
    near = min(avail, key=lambda e: abs(e - eta))
    return near, f"no {SLICE} entry at eta={eta}; nearest slice eta={near} used"


def probe(family, N, eta):
    if family.startswith("opt_"):
        e, _ = opt_slice(family, eta)
        psi, _ = opt.load_probe(family[4:], SLICE, e, 0.0, NBAR, N_basis=N)
        return psi
    return b.make_probe(family, N, NBAR)


def beta(eta):
    """Loss during a displacement = loss, then the displacement scaled by beta (paper Eq. 16)."""
    return 2 * (1 - np.sqrt(eta)) / abs(np.log(eta))


def reference(family, eta):
    """The kappa -> 0 (no back-action) value the sweep must reduce to at high frequency."""
    r = np.arcsinh(np.sqrt(NBAR))
    if family == "squeezed":
        return beta(eta) ** 2 * 4 / (eta * np.exp(-2 * r) + 1 - eta), "closed form: beta² × paper Eq. 17"
    if family == "coherent":
        return 4 * beta(eta) ** 2, "closed form: 4 beta²"
    if family.startswith("opt_"):
        e, _ = opt_slice(family, eta)
        return opt.best_entry(family[4:], SLICE, e, 0.0, NBAR)["qfi"], f"stored QFI (paper data, eta={e})"
    return b.rp_qfi(probe(family, 200, eta), kappa_ba=0.0, eta_ch=eta), "no-back-action Fock value"


def sweep(etas, n_max, quick, checkpoint=None):
    freqs = FREQS_HZ[::4] if quick else FREQS_HZ
    kappas = rp_kappa(2 * np.pi * freqs)
    ladder = [N for N in default_ladder(40, n_max) if N <= n_max]
    rows = []
    for eta in etas:
        for fam in FAMILIES:
            ref, ref_src = reference(fam, eta)
            e_used, e_note = opt_slice(fam, eta) if fam.startswith("opt_") else (np.nan, "")
            nbar = b.mean_photon_number(probe(fam, N_REF, eta))
            print(f"[eta_ch={eta}] [{fam}] <n>={nbar:.4f}  kappa->0 reference {ref:.4f} ({ref_src}) {e_note}", flush=True)
            for f, k in zip(freqs, kappas):
                t0 = time.time()
                gauss = g.gaussian_qfi(GAUSSIAN[fam](NBAR), kappa_ba=k, eta_ch=eta) if fam in GAUSSIAN else np.nan
                # primary: the exact shear-frame solver (compact at any kappa); cross-check against
                # the package's own density-matrix channel where that is cheap (small kappa)
                res = b.converge(lambda N: b.rp_qfi_unsheared(probe(fam, N, eta), kappa_ba=k, eta_ch=eta, N_basis=N),
                                 cutoffs=FRAME_LADDER, rtol=1e-3, n_stable=2)
                track = "frame"
                predicted = b.suggested_cutoff(NBAR, k, safety=2.0)
                mes = np.nan
                if predicted <= n_max:
                    rungs = [N for N in ladder if N >= max(40, int(0.6 * predicted))] or [n_max]
                    if len(rungs) < 2:
                        rungs = [rungs[0], n_max] if rungs[0] < n_max else [int(0.75 * n_max), n_max]
                    chk = b.converge(lambda N: b.rp_qfi(probe(fam, N, eta), kappa_ba=k, eta_ch=eta, N_basis=N),
                                     cutoffs=rungs, rtol=1e-3)
                    mes = chk.value if chk.converged else np.nan
                rows.append(dict(eta_ch=eta, family=fam, freq_hz=f, kappa_ba=k, nbar=nbar, qfi_fock=res.value,
                                 cutoff=res.cutoff, converged=res.converged, track=track, qfi_reference=ref,
                                 reference_source=ref_src, ratio_to_reference=res.value / ref,
                                 qfi_gaussian=gauss, qfi_mesolve=mes, opt_eta_used=e_used, note=e_note, seconds=time.time() - t0))
                flag = "" if res.converged else "  NOT CONVERGED"
                chk_s = f"  mesolve {mes:.4f} ({(mes - res.value) / res.value:+.1e})" if np.isfinite(mes) else ""
                print(f"   f={f:7.1f} Hz kappa={k:8.4f}  F={res.value:9.4f} @N={res.cutoff:<4d} [{track}] "
                      f"(ref {ref:.4f}, ratio {rows[-1]['ratio_to_reference']:.4f}){chk_s}  {rows[-1]['seconds']:.0f}s{flag}", flush=True)
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
    d_eta = df[df.eta_ch == eta]
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
    ax.set_title(f"backaction + channel loss ($\\eta_{{ch}}$ = {eta:g}) → displacement, ⟨n⟩ = 5")
    ax.grid(True, which="both", alpha=0.3)
    _kappa_axis(ax)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    print("figure:", path)


def write_summary(df, path):
    lines = ["# Channel loss (concurrent with the shear) + back-action, ⟨n⟩ = 5, 10 Hz – 1 kHz", "",
             "Loss during the sensing stage does not commute with the shear, so the QFI depends on",
             "frequency. At high frequency (κ → 0) it reduces to the paper's loss-during-displacement",
             "numbers: stored QFI of the loss_ch optimised states; β² × loss-then-displacement closed",
             "forms for squeezed vacuum and coherent light; exact κ = 0 Fock values for the rest.", ""]
    for eta in sorted(df.eta_ch.unique()):
        d_eta = df[df.eta_ch == eta]
        lines += [f"## η_ch = {eta:g}", "",
                  "| probe | κ→0 reference | source | F at 1 kHz | F at 30 Hz (κ≈1) | F at 10 Hz (κ≈9.5) | 10 Hz / 1 kHz | converged | package-channel cross-check (points, worst rel.) | largest cutoff |",
                  "|---|---|---|---|---|---|---|---|---|---|"]
        for fam in [f for f in FAMILIES if f in set(d_eta.family)]:
            d = d_eta[d_eta.family == fam].sort_values("freq_hz"); ok = d[d.converged]
            at = lambda f0: d.iloc[(d.freq_hz - f0).abs().argmin()].qfi_fock
            m = d[np.isfinite(d.qfi_mesolve)]
            tr = f"{len(m)} pts, {((m.qfi_mesolve - m.qfi_fock) / m.qfi_fock).abs().max():.1e}" if len(m) else "–"
            lines.append(f"| {LEGEND[fam]} | {d.qfi_reference.iloc[0]:.4f} | {d.reference_source.iloc[0]} | "
                         f"{at(1000):.4f} | {at(30):.4f} | {at(10):.4f} | {at(10)/at(1000):.3f} | "
                         f"{len(ok)}/{len(d)} | {tr} | {int(d.cutoff.max())} |")
        notes = d_eta[d_eta.note != ""].drop_duplicates("family")
        for _, r in notes.iterrows():
            lines.append(f"\n*{LEGEND[r.family]}: {r.note}.*")
        lines.append("")
    path.write_text("\n".join(lines) + "\n")
    print("summary:", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-max", type=int, default=120, help="density-matrix cross-check budget")
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
    for eta in sorted(df.eta_ch.unique()):
        figure(df, eta, FIG / (f"channel_loss_eta{eta:g}".replace(".", "p") + ".png"))


if __name__ == "__main__":
    sys.exit(main())
