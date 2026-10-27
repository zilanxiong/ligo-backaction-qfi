"""Output (detection) loss + back-action, nbar = 5, 10 Hz - 1 kHz, transmission 0.95.

Loss acts *after* the sensing stage, on the sheared state.  At kappa -> 0 it is exactly
eta times the input-loss value for any probe (loss maps a displaced state to one displaced
by sqrt(eta) as much): the stored QFI of the loss_out optimised states, eta x the
loss-then-displacement closed forms for squeezed vacuum and coherent light.  With the shear
on, the vacuum noise mixed in after the stretch leaves a residue of the shear in the phase
quadrature that no measurement can remove: F = 4 eta / (1 + eta (1 - eta) kappa^2) for the
vacuum.  This is the placement behind frequency-dependent squeezing.

Fock track: the exact shear-frame solver (backaction_qfi.frame), cross-checked against
the package's own density-matrix channel at the frequencies where that is cheap.

Usage:  python studies/output_loss.py [--pn-in 0.2] [--etas 0.95 ...] [--quick] [--figure-only]
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
OUT = ROOT / "results" / "output_loss"   # suffixed with the phase noise in main()
FIG = ROOT / "figures"

NBAR = 5.0
ETAS = [0.95]
PN_IN = 0.0            # set from --pn-in
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
SLICE = "loss_out"
FALLBACK_SLICE = "loss_ch"     # the paper has no loss_out optimisations for sqz_coh


def _available_etas(state_type, config):
    df = pd.read_csv(opt.DATA_DIR / "master_results.csv", dtype={"valid": str})
    df = df[(df.valid.fillna("").str.lower() != "false") & (df.state_type == state_type)
            & (df.loss_config == config) & (df.pn == PN_IN) & (df.N_target == NBAR)]
    return sorted(df.eta.unique())


def opt_slice(family, eta):
    """(config, eta used, note): the matching loss_out slice; else the fallback config
    (loss_ch) at the same eta; else the nearest eta."""
    for config in (SLICE, FALLBACK_SLICE):
        avail = _available_etas(family[4:], config)
        if not avail:
            continue
        if eta in avail:
            return config, eta, ("" if config == SLICE else f"no {SLICE} optimisation for this family; {config} state at eta={eta} used")
        near = min(avail, key=lambda e: abs(e - eta))
        return config, near, f"no {config} entry at eta={eta}, pn={PN_IN}; nearest slice eta={near} used"
    raise LookupError(family)


def probe(family, N, eta):
    if family.startswith("opt_"):
        cfg, e, _ = opt_slice(family, eta)
        psi, _ = opt.load_probe(family[4:], cfg, e, PN_IN, NBAR, N_basis=N)
        return psi
    return b.make_probe(family, N, NBAR)


def reference(family, eta):
    """The kappa -> 0 value: loss after the signal = eta x loss before the signal, exactly."""
    r = np.arcsinh(np.sqrt(NBAR))
    if family == "squeezed" and PN_IN == 0.0:
        return eta * 4 / (eta * np.exp(-2 * r) + 1 - eta), "closed form: eta × paper Eq. 17"
    if family == "coherent" and PN_IN == 0.0:
        return 4 * eta, "closed form: 4 eta"
    if family.startswith("opt_"):
        cfg, e, _ = opt_slice(family, eta)
        if cfg == SLICE:
            return opt.best_entry(family[4:], cfg, e, PN_IN, NBAR)["qfi"], f"stored QFI (paper data, {cfg}, eta={e}, pn={PN_IN})"
    return b.rp_qfi(probe(family, 200, eta), kappa_ba=0.0, eta_out=eta, pn_in=PN_IN), "no-back-action Fock value"


def sweep(etas, n_max, quick, checkpoint=None):
    freqs = FREQS_HZ[::4] if quick else FREQS_HZ
    kappas = rp_kappa(2 * np.pi * freqs)
    ladder = [N for N in default_ladder(40, n_max) if N <= n_max]
    rows = []
    for eta in etas:
        for fam in FAMILIES:
            ref, ref_src = reference(fam, eta)
            _, e_used, e_note = opt_slice(fam, eta) if fam.startswith("opt_") else (None, np.nan, "")
            nbar = b.mean_photon_number(probe(fam, N_REF, eta))
            print(f"[eta_out={eta} pn_in={PN_IN}] [{fam}] <n>={nbar:.4f}  kappa->0 reference {ref:.4f} ({ref_src}) {e_note}", flush=True)
            for f, k in zip(freqs, kappas):
                t0 = time.time()
                gauss = g.gaussian_qfi(GAUSSIAN[fam](NBAR), kappa_ba=k, eta_out=eta) if (fam in GAUSSIAN and PN_IN == 0.0) else np.nan
                # primary: the exact shear-frame solver (compact at any kappa); cross-check against
                # the package's own density-matrix channel where that is cheap (small kappa)
                res = b.converge(lambda N: b.rp_qfi_unsheared(probe(fam, N, eta), kappa_ba=k, eta_out=eta, pn_in=PN_IN, N_basis=N),
                                 cutoffs=FRAME_LADDER, rtol=1e-3, n_stable=2)
                track = "frame"
                predicted = b.suggested_cutoff(NBAR, k, safety=2.0)
                mes = np.nan
                if predicted <= n_max:
                    rungs = [N for N in ladder if N >= max(40, int(0.6 * predicted))] or [n_max]
                    if len(rungs) < 2:
                        rungs = [rungs[0], n_max] if rungs[0] < n_max else [int(0.75 * n_max), n_max]
                    chk = b.converge(lambda N: b.rp_qfi(probe(fam, N, eta), kappa_ba=k, eta_out=eta, pn_in=PN_IN, N_basis=N),
                                     cutoffs=rungs, rtol=1e-3, n_stable=2)
                    mes = chk.value if chk.converged else np.nan
                rows.append(dict(eta_out=eta, pn_in=PN_IN, family=fam, freq_hz=f, kappa_ba=k, nbar=nbar, qfi_fock=res.value,
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
    d_eta = df[df.eta_out == eta]
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
    pn = d_eta.pn_in.iloc[0] if "pn_in" in d_eta else 0.0
    head = f"phase noise ({1000 * pn:.0f} mrad) → " if pn > 0 else ""
    ax.set_title(head + f"backaction + displacement → output loss ($\\eta_{{out}}$ = {eta:g}), ⟨n⟩ = 5", fontsize=11)
    ax.grid(True, which="both", alpha=0.3)
    _kappa_axis(ax)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    print("figure:", path)


def write_summary(df, path):
    lines = ["# Output (detection) loss + back-action, ⟨n⟩ = 5, 10 Hz – 1 kHz", "",
             "Loss after the sensing stage acts on the sheared state. At κ → 0 it is exactly η × the",
             "input-loss value: stored QFI of the loss_out optimised states; η × the closed forms for",
             "squeezed vacuum and coherent light; exact κ = 0 Fock values for the rest.", ""]
    for eta in sorted(df.eta_out.unique()):
        d_eta = df[df.eta_out == eta]
        pn = d_eta.pn_in.iloc[0] if "pn_in" in d_eta else 0.0
        lines += [f"## η_out = {eta:g}" + (f", input phase noise {1000 * pn:.0f} mrad" if pn > 0 else ""), "",
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
    ap.add_argument("--pn-in", type=float, default=0.0, help="input phase noise rms [rad]")
    args = ap.parse_args()
    global PN_IN, OUT
    PN_IN = args.pn_in
    tag = f"_pn{PN_IN:g}".replace(".", "p") if PN_IN > 0 else ""
    OUT = ROOT / "results" / ("output_loss" + tag)
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)
    csv = OUT / "qfi_vs_frequency.csv"
    if args.figure_only:
        df = pd.read_csv(csv); df["note"] = df["note"].fillna("")
    else:
        df = sweep(args.etas, args.n_max, args.quick, checkpoint=OUT / "checkpoint.csv")
        df.to_csv(csv, index=False)
        (OUT / "checkpoint.csv").unlink(missing_ok=True)
    write_summary(df, OUT / "summary.md")
    for eta in sorted(df.eta_out.unique()):
        figure(df, eta, FIG / (f"output_loss_eta{eta:g}".replace(".", "p") + tag + ".png"))


if __name__ == "__main__":
    sys.exit(main())
