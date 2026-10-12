"""First sweep: back-action only (no loss, no phase noise), nbar = 5, 10 Hz - 1 kHz.

For every probe family the sweep reports
  * the QFI for the GW signal epsilon_a at each frequency, from the package's
    density-matrix channel with a convergence-checked Fock cutoff (Fock track),
    next to the exact lossless value 8 Var(x) and, for Gaussian probes, the
    covariance-matrix track;
  * the quantum Cramér-Rao bound on strain, sigma_h/h_SQL = sqrt(2/(kappa F)),
    which *is* frequency dependent because kappa(Omega) enters the conversion;
  * the scaling of the QFI with mean photon number.

Usage:  python studies/backaction_only.py [--n-max 600] [--quick]
Writes results/backaction_only/*.csv, results/backaction_only/summary.md and
figures/backaction_only.{png,svg}.
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
NBARS = [1, 2, 3, 5, 8, 12, 20]
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


def optimised_nbars(state_type):
    df = pd.read_csv(opt.DATA_DIR / "master_results.csv", dtype={"valid": str})
    df = df[(df.valid.fillna("").str.lower() != "false") & (df.state_type == state_type)
            & (df.loss_config == "loss_ch") & (df.eta == 1.0) & (df.pn == 0.0)]
    return sorted(df.N_target.unique())


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
                                 rel_err_fock=np.nan, sigma_h_over_sql=np.sqrt(2.0 / (k * exact)),
                                 sigma_h=np.sqrt(2.0 / (k * exact)) * g.h_sql(2 * np.pi * f), seconds=0.0,
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
                             sigma_h_over_sql=np.sqrt(2.0 / (k * exact)),
                             sigma_h=np.sqrt(2.0 / (k * exact)) * g.h_sql(2 * np.pi * f),
                             seconds=time.time() - t0, note=""))
            flag = "" if res.converged else "  NOT CONVERGED"
            print(f"   f={f:7.1f} Hz  kappa={k:8.4f}  F_fock={res.value:9.4f} @N={res.cutoff:<4d} "
                  f"(exact {exact:.4f}, rel {rows[-1]['rel_err_fock']:+.1e})  {rows[-1]['seconds']:.0f}s{flag}", flush=True)
    return pd.DataFrame(rows)


def sweep_nbar(families):
    rows = []
    for fam in families:
        if fam.startswith("opt_"):
            nbars = optimised_nbars(fam[4:])
        elif fam == "fock":
            nbars = [n for n in NBARS]
        else:
            nbars = NBARS
        for n in nbars:
            try:
                psi = probe(fam, N_REF, n)
            except (ValueError, LookupError) as e:
                print(f"   skip {fam} nbar={n}: {e}")
                continue
            tail = b.tail_population(psi)
            rows.append(dict(family=fam, nbar_target=n, nbar=b.mean_photon_number(psi),
                             var_x=b.var_x(psi), qfi_exact=8 * b.var_x(psi), tail_population=tail))
    return pd.DataFrame(rows)


def make_figure(df_f, df_n, path_png, path_svg):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    palette = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
    color = dict(zip(FAMILIES, palette))
    surface, ink, ink2, grid = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": ink2, "axes.labelcolor": ink,
                         "xtick.color": ink2, "ytick.color": ink2, "text.color": ink,
                         "axes.facecolor": surface, "figure.facecolor": surface,
                         "axes.grid": True, "grid.color": grid, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.6), gridspec_kw={"wspace": 0.75, "width_ratios": [1.15, 1, 1]})
    fams = [f for f in FAMILIES if f in set(df_f.family)]

    # (a) QFI vs frequency: exact line, Fock-track markers
    ax = axes[0]
    for fam in fams:
        d = df_f[df_f.family == fam].sort_values("freq_hz")
        ax.plot(d.freq_hz, d.qfi_exact, color=color[fam], lw=2, label=LABELS[fam])
        ok = d[d.converged]
        ax.plot(ok.freq_hz, ok.qfi_fock, "o", ms=5, color=color[fam], mec=surface, mew=1)
        bad = d[~d.converged]
        if len(bad):
            ax.plot(bad.freq_hz, bad.qfi_fock, "x", ms=6, color=color[fam])
    ax.set_xscale("log")
    ax.set_xlabel("sideband frequency  Ω/2π  [Hz]")
    ax.set_ylabel("QFI for the GW signal  F(ε_a)")
    ax.set_title("(a) QFI vs frequency, back-action only, ⟨n⟩ = 5", loc="left", fontsize=10)
    ax.text(0.02, 0.63, "lines: exact 8 Var(x)\n●  Fock track, converged cutoff\n×  Fock track, not converged",
            transform=ax.transAxes, fontsize=8, color=ink2, va="top")
    # direct labels at the right edge, stacked so they never overlap; a thin
    # leader joins a displaced label to its line
    vals = sorted(((df_f[df_f.family == f].qfi_exact.iloc[0], f) for f in fams), reverse=True)
    ymax = vals[0][0]
    ax.set_ylim(0, ymax * 1.08)
    spacing = ymax * 1.08 * 0.038            # ~ one 7.5-pt line of text in data units
    y_prev = np.inf
    for v, f in vals:
        y = min(v, y_prev - spacing)
        if abs(y - v) > 0.3 * spacing:
            ax.plot([1000, 1400], [v, y], color=grid, lw=0.8, clip_on=False, zorder=1)
        ax.annotate(f"{LABELS[f]}  {v:.1f}", xy=(1400, y), xytext=(4, 0), textcoords="offset points",
                    fontsize=7.5, color=ink, va="center", ha="left", annotation_clip=False)
        y_prev = y

    # (b) strain bound vs frequency
    ax = axes[1]
    for fam in fams:
        d = df_f[df_f.family == fam].sort_values("freq_hz")
        ax.plot(d.freq_hz, d.sigma_h_over_sql, color=color[fam], lw=2, label=LABELS[fam])
    kk = rp_kappa(2 * np.pi * FREQS_HZ)
    ax.plot(FREQS_HZ, np.ones_like(FREQS_HZ), "--", color=ink2, lw=1.2)
    ax.text(FREQS_HZ[-1], 1.05, "SQL", ha="right", va="bottom", fontsize=8, color=ink2)
    ax.plot(FREQS_HZ, np.sqrt((1 + kk**2) / (2 * kk)), ":", color=ink2, lw=1.2)
    ax.text(FREQS_HZ[0], np.sqrt((1 + kk[0]**2) / (2 * kk[0])) * 1.1, "vacuum, phase readout", fontsize=8, color=ink2)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("sideband frequency  Ω/2π  [Hz]")
    ax.set_ylabel("strain bound  σ_h / h_SQL  = √(2/(κF))")
    ax.set_title("(b) quantum Cramér-Rao bound on strain, relative to the SQL", loc="left", fontsize=10)

    # (c) QFI vs nbar
    ax = axes[2]
    for fam in [f for f in FAMILIES if f in set(df_n.family)]:
        d = df_n[df_n.family == fam].sort_values("nbar")
        ax.plot(d.nbar, d.qfi_exact, "-o", ms=4, lw=2, color=color[fam], mec=surface, mew=0.8, label=LABELS[fam])
    nn = np.logspace(np.log10(0.5), np.log10(40), 50)
    ax.plot(nn, 8 * nn + 4, "--", color=ink2, lw=1); ax.text(0.55, 8 * 0.55 + 4, "8n+4", fontsize=8, color=ink2, va="top", ha="left")
    ax.plot(nn, 4 * (1 + 2 * nn + 2 * np.sqrt(nn * (nn + 1))), ":", color=ink2, lw=1)
    ax.text(30, 4 * (1 + 60 + 2 * np.sqrt(30 * 31)), "4(1+2n+2√(n(n+1)))", fontsize=8, color=ink2, va="bottom", ha="right")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("mean photon number  ⟨n⟩")
    ax.set_ylabel("lossless QFI  F(ε_a) = 8 Var(x)")
    ax.set_title("(c) scaling with photon number, lossless", loc="left", fontsize=10)

    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 0.0))
    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.24, top=0.9, wspace=0.75)
    fig.savefig(path_png, dpi=160); fig.savefig(path_svg)
    print("figure:", path_png)


def write_summary(df_f, df_n, path):
    lines = ["# Back-action only, ⟨n⟩ = 5, 10 Hz – 1 kHz", "",
             "Channel: unitary, `H = √2 ε_a x + (κ/2) x²`, no loss, no phase noise. κ(Ω) from",
             "`quantum_sensing.rp_kappa` (package LIGO defaults). QFI for `ε_a` at `ε_a = 0`.", "",
             "## Result", "",
             "The QFI is **frequency independent** for every probe: `x²` commutes with the signal",
             "generator `x`, so `F = 8 Var(x)` of the input, exactly. The Fock track reproduces this",
             "at every frequency where its cutoff ladder converged. The strain bound",
             "`σ_h/h_SQL = √(2/(κF))` is frequency dependent through κ(Ω) only.", "",
             "| probe | ⟨n⟩ | F exact (8 Var x) | Fock-track range over f | worst rel. err. | converged points | largest cutoff | σ_h/h_SQL at 30 Hz |",
             "|---|---|---|---|---|---|---|---|"]
    for fam in [f for f in FAMILIES if f in set(df_f.family)]:
        d = df_f[df_f.family == fam]
        ok = d[d.converged]
        i30 = (d.freq_hz - 30).abs().idxmin()
        lines.append(f"| {LABELS[fam]} | {d.nbar.iloc[0]:.3f} | {d.qfi_exact.iloc[0]:.3f} | "
                     f"{ok.qfi_fock.min():.3f} – {ok.qfi_fock.max():.3f} | {ok.rel_err_fock.abs().max():.1e} | "
                     f"{len(ok)}/{len(d)} ({int((d.note != '').sum())} beyond budget) | {d.cutoff.max()} | {d.loc[i30, 'sigma_h_over_sql']:.3f} |")
    lines += ["", "## Photon-number scaling (exact, lossless)", "",
              "| probe | ⟨n⟩ values | F at ⟨n⟩ = 5 | F at largest ⟨n⟩ | local slope at top |", "|---|---|---|---|---|"]
    for fam in [f for f in FAMILIES if f in set(df_n.family)]:
        d = df_n[df_n.family == fam].sort_values("nbar")
        f5 = d.iloc[(d.nbar - 5).abs().argmin()].qfi_exact
        slope = np.polyfit(np.log(d.nbar.values[-2:]), np.log(d.qfi_exact.values[-2:]), 1)[0] if len(d) > 1 else np.nan
        lines.append(f"| {LABELS[fam]} | {', '.join(f'{v:g}' for v in d.nbar_target)} | {f5:.2f} | {d.qfi_exact.iloc[-1]:.2f} | {slope:.2f} |")
    path.write_text("\n".join(lines) + "\n")
    print("summary:", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-max", type=int, default=600, help="largest Fock cutoff the ladder may use")
    ap.add_argument("--quick", action="store_true", help="every 4th frequency only")
    ap.add_argument("--families", nargs="*", default=FAMILIES)
    ap.add_argument("--figure-only", action="store_true", help="re-render the figure and summary from saved CSVs")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)

    if args.figure_only:
        df_f = pd.read_csv(OUT / "qfi_vs_frequency.csv"); df_f["note"] = df_f["note"].fillna("")
        df_n = pd.read_csv(OUT / "qfi_vs_nbar.csv")
    else:
        df_f = sweep_frequency(args.families, args.n_max, args.quick)
        df_f.to_csv(OUT / "qfi_vs_frequency.csv", index=False)
        df_n = sweep_nbar(args.families)
        df_n.to_csv(OUT / "qfi_vs_nbar.csv", index=False)
    write_summary(df_f, df_n, OUT / "summary.md")
    make_figure(df_f, df_n, FIG / "backaction_only.png", FIG / "backaction_only.svg")


if __name__ == "__main__":
    sys.exit(main())
