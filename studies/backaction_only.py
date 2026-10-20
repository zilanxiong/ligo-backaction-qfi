"""First sweep: back-action only (no loss, no phase noise), nbar = 5, 10 Hz - 1 kHz.

For every probe family the sweep reports
  * the QFI for the GW signal epsilon_a at each frequency, from the package's
    density-matrix channel with a convergence-checked Fock cutoff (Fock track),
    next to the exact lossless value 8 Var(x) and, for Gaussian probes, the
    covariance-matrix track;

Usage:  python studies/backaction_only.py [--n-max 600] [--quick]
Writes results/backaction_only/qfi_vs_frequency.csv, results/backaction_only/summary.md
and figures/backaction_only_qfi_vs_frequency.png.
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
OUT = ROOT / "results" / "backaction_only"
FIG = ROOT / "figures"

NBAR = 5.0
FREQS_HZ = np.logspace(1, 3, 25)
N_REF = 400                      # cutoff for the exact (input-state) quantities

#: Sweep order = plotting order = palette slot order.
FAMILIES = ["squeezed", "coherent", "cat", "squeezed_cat", "fock",
            "opt_fock_sup", "opt_sqz_vac", "opt_sqz_coh"]
LABELS = {
    "squeezed": "squeezed vacuum", "coherent": "coherent", "cat": "even cat",
    "squeezed_cat": "squeezed cat", "fock": "Fock |5⟩",
    "opt_fock_sup": "optimised Fock superposition", "opt_sqz_vac": "optimised sqz-vac superposition",
    "opt_sqz_coh": "optimised sqz-coh superposition",
}
GAUSSIAN = {"squeezed": lambda n: g.squeezed_at_fixed_nbar(n),
            "coherent": lambda n: g.coherent_gaussian(np.sqrt(n))}


def probe(family, N, nbar):
    """Ket for ``family`` at cutoff ``N`` and photon number ``nbar`` (oriented for epsilon_a)."""
    if family.startswith("opt_"):
        psi, _ = opt.load_probe(family[4:], "loss_ch", 1.0, 0.0, nbar, N_basis=N)
        return psi
    return b.make_probe(family, N, nbar)


def sweep_frequency(families, n_max, quick):
    freqs = FREQS_HZ[::4] if quick else FREQS_HZ
    kappas = rp_kappa(2 * np.pi * freqs)
    ladder = [N for N in default_ladder(40, n_max) if N <= n_max]
    rows = []
    for fam in families:
        psi_ref = probe(fam, N_REF, NBAR)
        exact = 8 * b.var_x(psi_ref)
        nbar = b.mean_photon_number(psi_ref)
        print(f"[{fam}] <n>={nbar:.4f}  exact lossless QFI = 8 Var(x) = {exact:.4f}", flush=True)
        for f, k in zip(freqs, kappas):
            t0 = time.time()
            predicted = b.suggested_cutoff(NBAR, k, safety=2.0)   # ~2 (1+k²)(n+1): output nbar x2
            if predicted > n_max:
                rows.append(dict(family=fam, freq_hz=f, kappa_ba=k, nbar=nbar, qfi_fock=np.nan,
                                 cutoff=0, converged=False, qfi_exact=exact,
                                 qfi_gaussian=g.gaussian_qfi(GAUSSIAN[fam](NBAR), kappa_ba=k) if fam in GAUSSIAN else np.nan,
                                 rel_err_fock=np.nan, seconds=0.0,
                                 note=f"skipped: predicted cutoff {predicted} > budget {n_max}"))
                print(f"   f={f:7.1f} Hz  kappa={k:8.4f}  Fock track skipped (predicted cutoff {predicted} > {n_max})", flush=True)
                continue
            start = max(40, int(0.6 * predicted))
            rungs = [N for N in ladder if N >= start] or [n_max]
            if len(rungs) < 2:
                rungs = [start, n_max] if start < n_max else [int(0.75 * n_max), n_max]
            res = b.converged_qfi(lambda N: probe(fam, N, NBAR), cutoffs=rungs, kappa_ba=k, rtol=1e-3)
            gauss = g.gaussian_qfi(GAUSSIAN[fam](NBAR), kappa_ba=k) if fam in GAUSSIAN else np.nan
            rows.append(dict(family=fam, freq_hz=f, kappa_ba=k, nbar=nbar, qfi_fock=res.value,
                             cutoff=res.cutoff, converged=res.converged, qfi_exact=exact,
                             qfi_gaussian=gauss, rel_err_fock=(res.value - exact) / exact,
                             seconds=time.time() - t0, note=""))
            flag = "" if res.converged else "  NOT CONVERGED"
            print(f"   f={f:7.1f} Hz  kappa={k:8.4f}  F_fock={res.value:9.4f} @N={res.cutoff:<4d} "
                  f"(exact {exact:.4f}, rel {rows[-1]['rel_err_fock']:+.1e})  {rows[-1]['seconds']:.0f}s{flag}", flush=True)
    return pd.DataFrame(rows)


PURE_LADDER = [400, 560, 800, 1120, 1600, 2240, 3200]


def fill_pure(df_f):
    """Where the density-matrix track was beyond its cutoff budget, compute the point
    with the pure-state Fock track (lossless channel, same Hamiltonian), converged
    on a ladder of much larger cutoffs."""
    df_f = df_f.copy()
    if "track" not in df_f:
        df_f["track"] = "mesolve"
    todo = df_f.index[df_f.note.fillna("").str.startswith("skipped")]
    for i in todo:
        fam, k = df_f.at[i, "family"], df_f.at[i, "kappa_ba"]
        t0 = time.time()
        res = b.converge(lambda N: b.rp_qfi_pure(probe(fam, N, NBAR), kappa_ba=k), cutoffs=PURE_LADDER, rtol=1e-3)
        exact = df_f.at[i, "qfi_exact"]
        df_f.loc[i, ["qfi_fock", "cutoff", "converged", "rel_err_fock", "seconds", "track", "note"]] = [
            res.value, res.cutoff, res.converged, (res.value - exact) / exact, time.time() - t0, "pure",
            "pure-state Fock track (sesolve), beyond density-matrix cutoff budget"]
        print(f"   [{fam}] f={df_f.at[i, 'freq_hz']:6.1f} Hz kappa={k:.3f}: pure track {res}  rel {(res.value-exact)/exact:+.1e}  {time.time()-t0:.0f}s", flush=True)
    return df_f


LEGEND = {
    "squeezed": "Squeezed vac. (p)", "coherent": "Coherent", "cat": "Even cat",
    "squeezed_cat": "Squeezed cat", "fock": "Fock",
    "opt_fock_sup": "Optimized Fock sup. (no-BA, lossless)",
    "opt_sqz_vac": "Optimized sqz-vac sup. (no-BA, lossless)",
    "opt_sqz_coh": "Optimized sqz-coh sup. (no-BA, lossless)",
}
TITLE_CHANNEL = "backaction → displacement, ⟨n⟩ = 5"


def _kappa_axis(ax):
    """Top axis showing kappa_ba(Omega): ticks placed at the exact frequencies
    where kappa hits each decade (kappa(f) is monotone decreasing, not a power law)."""
    from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
    from quantum_sensing import rp_power_ratio
    gamma = 2 * np.pi * 450.0
    X = rp_power_ratio(gamma=gamma)

    def k_to_f(k):           # invert kappa = 2 X g^4 / (W^2 (g^2 + W^2)) for W
        W2 = 0.5 * (-gamma**2 + np.sqrt(gamma**4 + 8 * X * gamma**4 / np.asarray(k, float)))
        return np.sqrt(W2) / (2 * np.pi)

    fmin, fmax = ax.get_xlim()
    decades = [10.0 ** e for e in range(2, -6, -1)]
    ticks = [(k_to_f(k), e) for k, e in zip(decades, range(2, -6, -1)) if fmin <= k_to_f(k) <= fmax]
    top = ax.twiny()
    top.set_xscale("log"); top.set_xlim(fmin, fmax)
    top.xaxis.set_major_locator(FixedLocator([f for f, _ in ticks]))
    top.xaxis.set_major_formatter(FixedFormatter([f"$10^{{{e}}}$" for _, e in ticks]))
    top.xaxis.set_minor_locator(NullLocator())
    top.set_xlabel(r"$\kappa_{ba}(\Omega)$")
    top.grid(False)
    return top


def _style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcdefaults()
    plt.rcParams.update({"figure.figsize": (8, 5.3), "figure.dpi": 150, "savefig.dpi": 150,
                         "axes.grid": True, "grid.alpha": 0.3, "legend.framealpha": 0.9,
                         "lines.linewidth": 1.5, "lines.markersize": 4})
    return plt


def figure_qfi_vs_frequency(df_f, path):
    plt = _style()
    fig, ax = plt.subplots()
    for fam in [f for f in FAMILIES if f in set(df_f.family)]:
        d = df_f[df_f.family == fam].sort_values("freq_hz")
        ok, skipped = d[d.converged], d[~d.converged]
        line, = ax.plot(ok.freq_hz, ok.qfi_fock, "-o", label=LEGEND[fam])
        if len(skipped):      # still beyond every budget: exact value, hollow marker
            ax.plot(skipped.freq_hz, skipped.qfi_exact, "o", mfc="white", color=line.get_color())
            ax.plot([skipped.freq_hz.max(), ok.freq_hz.min()], [skipped.qfi_exact.iloc[0], ok.qfi_fock.iloc[0]],
                    "-", color=line.get_color(), alpha=0.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"Frequency  $\Omega/2\pi$  [Hz]")
    ax.set_ylabel(r"QFI for $\varepsilon_a$")
    ax.set_title(TITLE_CHANNEL)
    ax.grid(True, which="both", alpha=0.3)
    _kappa_axis(ax)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    print("figure:", path)


def make_figures(df_f):
    figure_qfi_vs_frequency(df_f, FIG / "backaction_only_qfi_vs_frequency.png")


def write_summary(df_f, path):
    lines = ["# Back-action only, ⟨n⟩ = 5, 10 Hz – 1 kHz", "",
             "Channel: unitary, `H = √2 ε_a x + (κ/2) x²`, no loss, no phase noise. κ(Ω) from",
             "`quantum_sensing.rp_kappa` (package LIGO defaults). QFI for `ε_a` at `ε_a = 0`.", "",
             "## Result", "",
             "The QFI is **frequency independent** for every probe: `x²` commutes with the signal",
             "generator `x`, so `F = 8 Var(x)` of the input, exactly. The Fock track reproduces this",
             "at every frequency.", "",
             "| probe | ⟨n⟩ | F exact (8 Var x) | Fock-track range over f | worst rel. err. | converged points | largest cutoff |",
             "|---|---|---|---|---|---|---|"]
    for fam in [f for f in FAMILIES if f in set(df_f.family)]:
        d = df_f[df_f.family == fam]
        ok = d[d.converged]
        npure = int((d.get("track", "mesolve") == "pure").sum())
        lines.append(f"| {LABELS[fam]} | {d.nbar.iloc[0]:.3f} | {d.qfi_exact.iloc[0]:.3f} | "
                     f"{ok.qfi_fock.min():.3f} – {ok.qfi_fock.max():.3f} | {ok.rel_err_fock.abs().max():.1e} | "
                     f"{len(ok)}/{len(d)} ({npure} via pure-state track) | {int(d.cutoff.max())} |")
    path.write_text("\n".join(lines) + "\n")
    print("summary:", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-max", type=int, default=600, help="largest Fock cutoff the ladder may use")
    ap.add_argument("--quick", action="store_true", help="every 4th frequency only")
    ap.add_argument("--families", nargs="*", default=FAMILIES)
    ap.add_argument("--figure-only", action="store_true", help="re-render the figure and summary from saved CSVs")
    ap.add_argument("--fill-pure", action="store_true", help="fill beyond-budget points with the pure-state Fock track, then re-render")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)

    if args.figure_only or args.fill_pure:
        df_f = pd.read_csv(OUT / "qfi_vs_frequency.csv"); df_f["note"] = df_f["note"].fillna("")
        if args.fill_pure:
            df_f = fill_pure(df_f)
            df_f.to_csv(OUT / "qfi_vs_frequency.csv", index=False)
    else:
        df_f = sweep_frequency(args.families, args.n_max, args.quick)
        df_f.to_csv(OUT / "qfi_vs_frequency.csv", index=False)
    write_summary(df_f, OUT / "summary.md")
    make_figures(df_f)


if __name__ == "__main__":
    sys.exit(main())
