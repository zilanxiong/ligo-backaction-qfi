"""Combinations of input, channel and output loss with back-action, nbar = 5, 10 Hz - 1 kHz.

The paper has no input+channel loss results; what it provides is the state-independent
identity (supplement Eq. 16) that loss *during* the displacement equals loss *before* it
with the displacement scaled by beta = 2(1 - sqrt(eta))/|ln eta|.  So at kappa -> 0 this
channel must reduce, for every probe, to

    F = eta_out x beta(eta_ch)^2 x F_input-loss(eta_in * eta_ch * eta_out)

(every loss commutes with the displacement up to a rescaling of the displacement, so all
of them move in front of the signal and combine into one loss of product transmission,
while the displacement is scaled by beta(eta_ch) * sqrt(eta_out)).  For squeezed
vacuum and coherent light that is a closed form (paper Eq. 17); for the others the sweep
records both the beta-relation value and the direct kappa = 0 value of the package channel.
The optimised states are the paper's channel-loss states at eta_ch (its only loss slice
with phase noise); they were not optimised for the extra input loss.

Usage:  python studies/combined_loss.py [--eta-in 0.95] [--eta-out 0.95] [--pn-in 0.2] [--etas 0.95 ...] [--quick] [--figure-only]
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
OUT = ROOT / "results" / "combined_loss"   # suffixed in main()
FIG = ROOT / "figures"

NBAR = 5.0
ETAS = [0.95]
PN_IN = 0.0            # set from --pn-in
ETA_IN = 0.95          # set from --eta-in (1.0 = off)
ETA_OUT = 1.0          # set from --eta-out (1.0 = off)
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


def _available_etas(state_type, config):
    df = pd.read_csv(opt.DATA_DIR / "master_results.csv", dtype={"valid": str})
    df = df[(df.valid.fillna("").str.lower() != "false") & (df.state_type == state_type)
            & (df.loss_config == config) & (df.pn == PN_IN) & (df.N_target == NBAR)]
    return sorted(df.eta.unique())


def _slice_for(eta):
    """Which of the paper's loss slices to draw the optimised states from."""
    if eta < 1.0:
        return "loss_ch", eta
    if ETA_IN < 1.0:
        return "loss_in", ETA_IN
    if ETA_OUT < 1.0:
        return "loss_out", ETA_OUT
    return "loss_ch", 1.0


def opt_slice(family, eta):
    config, e_want = _slice_for(eta)
    avail = _available_etas(family[4:], config)
    if not avail and config == "loss_out":          # sqz_coh has no loss_out optimisations
        config, avail = "loss_ch", _available_etas(family[4:], "loss_ch")
    if e_want in avail:
        return e_want, ("" if config == "loss_ch" and eta < 1.0 else f"{config} state at eta={e_want} used")
    near = min(avail, key=lambda e: abs(e - e_want))
    return near, f"no {config} entry at eta={e_want}, pn={PN_IN}; nearest {config} slice eta={near} used"


def probe(family, N, eta):
    if family.startswith("opt_"):
        e, note = opt_slice(family, eta)
        config = note.split(" state")[0] if note and " state at" in note else (note.split("nearest ")[1].split(" slice")[0] if "nearest" in note else "loss_ch")
        psi, _ = opt.load_probe(family[4:], config, e, PN_IN, NBAR, N_basis=N)
        return psi
    return b.make_probe(family, N, NBAR)


def beta(eta):
    """Loss during a displacement = loss, then the displacement scaled by beta (paper Eq. 16)."""
    return 1.0 if eta >= 1.0 else 2 * (1 - np.sqrt(eta)) / abs(np.log(eta))


def input_loss_qfi(family, eta, eta_pre):
    """QFI of the sweep's probe (slice ``eta``) with a single loss of transmission
    ``eta_pre`` *before* the signal, no back-action."""
    r = np.arcsinh(np.sqrt(NBAR))
    if family == "squeezed" and PN_IN == 0.0:
        return 4 / (eta_pre * np.exp(-2 * r) + 1 - eta_pre)
    if family == "coherent" and PN_IN == 0.0:
        return 4.0
    return b.rp_qfi(probe(family, 200, eta), kappa_ba=0.0, eta_in=eta_pre, pn_in=PN_IN)


def reference(family, eta):
    """kappa -> 0 value from the paper's beta relation: eta_out × beta(eta_ch)² × F_input-loss(eta_in eta_ch)."""
    src = ("closed form: eta_out × beta² × paper Eq. 17 at eta_in·eta_ch·eta_out" if (family in ("squeezed", "coherent") and PN_IN == 0.0)
           else "eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out")
    # every loss commutes with the displacement up to rescaling: all of them move before the
    # signal (product transmission), the displacement picks up beta(eta_ch) * sqrt(eta_out)
    return ETA_OUT * beta(eta) ** 2 * input_loss_qfi(family, eta, ETA_IN * eta * ETA_OUT), src


def direct_kappa0(family, eta):
    """The package channel itself at kappa = 0 with both losses (what the beta relation predicts)."""
    return b.rp_qfi(probe(family, 200, eta), kappa_ba=0.0, eta_in=ETA_IN, eta_ch=eta, eta_out=ETA_OUT, pn_in=PN_IN)


def sweep(etas, n_max, quick, checkpoint=None):
    freqs = FREQS_HZ[::4] if quick else FREQS_HZ
    kappas = rp_kappa(2 * np.pi * freqs)
    ladder = [N for N in default_ladder(40, n_max) if N <= n_max]
    rows = []
    for eta in etas:
        for fam in FAMILIES:
            ref, ref_src = reference(fam, eta)
            k0 = direct_kappa0(fam, eta)
            e_used, e_note = opt_slice(fam, eta) if fam.startswith("opt_") else (np.nan, "")
            nbar = b.mean_photon_number(probe(fam, N_REF, eta))
            print(f"[eta_in={ETA_IN} eta_ch={eta} eta_out={ETA_OUT} pn_in={PN_IN}] [{fam}] <n>={nbar:.4f}  kappa->0: beta relation {ref:.4f} vs direct {k0:.4f} "
                  f"(rel {(k0 - ref) / ref:+.1e}) ({ref_src}) {e_note}", flush=True)
            for f, k in zip(freqs, kappas):
                t0 = time.time()
                gauss = g.gaussian_qfi(GAUSSIAN[fam](NBAR), kappa_ba=k, eta_in=ETA_IN, eta_ch=eta, eta_out=ETA_OUT) if (fam in GAUSSIAN and PN_IN == 0.0) else np.nan
                # primary: the exact shear-frame solver (compact at any kappa); cross-check against
                # the package's own density-matrix channel where that is cheap (small kappa)
                res = b.converge(lambda N: b.rp_qfi_unsheared(probe(fam, N, eta), kappa_ba=k, eta_in=ETA_IN, eta_ch=eta, eta_out=ETA_OUT, pn_in=PN_IN, N_basis=N),
                                 cutoffs=FRAME_LADDER, rtol=1e-3, n_stable=2)
                track = "frame"
                predicted = b.suggested_cutoff(NBAR, k, safety=2.0)
                mes = np.nan
                if predicted <= n_max:
                    rungs = [N for N in ladder if N >= max(40, int(0.6 * predicted))] or [n_max]
                    if len(rungs) < 2:
                        rungs = [rungs[0], n_max] if rungs[0] < n_max else [int(0.75 * n_max), n_max]
                    chk = b.converge(lambda N: b.rp_qfi(probe(fam, N, eta), kappa_ba=k, eta_in=ETA_IN, eta_ch=eta, eta_out=ETA_OUT, pn_in=PN_IN, N_basis=N),
                                     cutoffs=rungs, rtol=1e-3)
                    mes = chk.value if chk.converged else np.nan
                rows.append(dict(eta_in=ETA_IN, eta_ch=eta, eta_out=ETA_OUT, pn_in=PN_IN, family=fam, freq_hz=f, kappa_ba=k, nbar=nbar, qfi_fock=res.value,
                                 cutoff=res.cutoff, converged=res.converged, track=track, qfi_reference=ref, qfi_kappa0_direct=k0,
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
    pn = d_eta.pn_in.iloc[0] if "pn_in" in d_eta else 0.0
    head = f"phase noise ({1000 * pn:.0f} mrad) → " if pn > 0 else ""
    ein, eout = d_eta.eta_in.iloc[0], d_eta.eta_out.iloc[0]
    title = head + (f"input loss ($\\eta_{{in}}$ = {ein:g}) → " if ein < 1 else "")
    title += "backaction + displacement" + (f" + channel loss ($\\eta_{{ch}}$ = {eta:g})" if eta < 1 else "")
    title += (f" → output loss ($\\eta_{{out}}$ = {eout:g})" if eout < 1 else "") + ", ⟨n⟩ = 5"
    ax.set_title(title, fontsize=10 if len(title) < 95 else 8.6)
    ax.grid(True, which="both", alpha=0.3)
    _kappa_axis(ax)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    print("figure:", path)


def write_summary(df, path):
    lines = ["# Loss combinations with back-action, ⟨n⟩ = 5, 10 Hz – 1 kHz", "",
             "The paper has no multi-placement loss results; the check is its state-independent β relation",
             "(supplement Eq. 16) plus the exact η_out scaling of loss after the signal: at κ → 0,",
             "F = η_out × β(η_ch)² × F_input-loss(η_in·η_ch·η_out). Closed form (Eq. 17) for",
             "squeezed vacuum and coherent light; for the rest the β relation applied to the package's own",
             "input-loss value, next to the package channel evaluated directly at κ = 0 with both losses.", ""]
    for eta in sorted(df.eta_ch.unique()):
        d_eta = df[df.eta_ch == eta]
        pn = d_eta.pn_in.iloc[0] if "pn_in" in d_eta else 0.0
        ein, eout = d_eta.eta_in.iloc[0], d_eta.eta_out.iloc[0]
        lines += [f"## η_in = {ein:g}, η_ch = {eta:g}, η_out = {eout:g}" + (f", input phase noise {1000 * pn:.0f} mrad" if pn > 0 else ""), "",
                  "| probe | κ→0 by β relation | source | direct κ=0 (package) | rel. diff | F at 1 kHz | F at 30 Hz (κ≈1) | F at 10 Hz (κ≈9.5) | 10 Hz / 1 kHz | converged | cross-check (points, worst rel.) | largest cutoff |",
                  "|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for fam in [f for f in FAMILIES if f in set(d_eta.family)]:
            d = d_eta[d_eta.family == fam].sort_values("freq_hz"); ok = d[d.converged]
            at = lambda f0: d.iloc[(d.freq_hz - f0).abs().argmin()].qfi_fock
            m = d[np.isfinite(d.qfi_mesolve)]
            tr = f"{len(m)} pts, {((m.qfi_mesolve - m.qfi_fock) / m.qfi_fock).abs().max():.1e}" if len(m) else "–"
            k0 = d.qfi_kappa0_direct.iloc[0]
            lines.append(f"| {LEGEND[fam]} | {d.qfi_reference.iloc[0]:.4f} | {d.reference_source.iloc[0]} | {k0:.4f} | "
                         f"{(k0 - d.qfi_reference.iloc[0]) / d.qfi_reference.iloc[0]:+.1e} | "
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
    ap.add_argument("--eta-in", type=float, default=0.95, help="input-loss transmission (1 = off)")
    ap.add_argument("--eta-out", type=float, default=1.0, help="output-loss transmission (1 = off)")
    args = ap.parse_args()
    global PN_IN, OUT, ETA_IN, ETA_OUT
    PN_IN, ETA_IN, ETA_OUT = args.pn_in, args.eta_in, args.eta_out
    ch = args.etas[0] if len(args.etas) == 1 else None
    tag = ((f"_in{ETA_IN:g}" if ETA_IN < 1 else "") + (f"_ch{ch:g}" if (ch is not None and ch < 1) else "")
           + (f"_out{ETA_OUT:g}" if ETA_OUT < 1 else "") + (f"_pn{PN_IN:g}" if PN_IN > 0 else "")).replace(".", "p")
    OUT = ROOT / "results" / ("combined_loss" + tag)
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
        figure(df, eta, FIG / ("combined_loss" + tag + ".png"))


if __name__ == "__main__":
    sys.exit(main())
