# Loss combinations with back-action, ⟨n⟩ = 5, 10 Hz – 1 kHz

The paper has no multi-placement loss results; the check is its state-independent β relation
(supplement Eq. 16) plus the exact η_out scaling of loss after the signal: at κ → 0,
F = η_out × β(η_ch)² × F_input-loss(η_in·η_ch·η_out). Closed form (Eq. 17) for
squeezed vacuum and coherent light; for the rest the β relation applied to the package's own
input-loss value, next to the package channel evaluated directly at κ = 0 with both losses.

## η_in = 0.95, η_ch = 1, η_out = 0.95

| probe | κ→0 by β relation | source | direct κ=0 (package) | rel. diff | F at 1 kHz | F at 30 Hz (κ≈1) | F at 10 Hz (κ≈9.5) | 10 Hz / 1 kHz | converged | cross-check (points, worst rel.) | largest cutoff |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Squeezed vac. (p) | 27.4155 | closed form: eta_out × beta² × paper Eq. 17 at eta_in·eta_ch·eta_out | 27.4155 | -5.5e-08 | 27.4155 | 20.7081 | 0.8140 | 0.030 | 25/25 | 19 pts, 2.0e-03 | 220 |
| Coherent | 3.8000 | closed form: eta_out × beta² × paper Eq. 17 at eta_in·eta_ch·eta_out | 3.8000 | -1.4e-09 | 3.8000 | 3.6442 | 0.7151 | 0.188 | 25/25 | 21 pts, 6.2e-06 | 220 |
| Even cat | 13.5600 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 13.5600 | +5.8e-09 | 13.5600 | 8.0682 | 0.7152 | 0.053 | 25/25 | 21 pts, 1.6e-05 | 220 |
| Squeezed cat | 27.3332 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 27.3332 | +4.6e-11 | 27.3332 | 20.6493 | 0.8139 | 0.030 | 25/25 | 19 pts, 5.7e-04 | 300 |
| Fock | 18.3195 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 18.3195 | +1.0e-09 | 18.3195 | 14.7582 | 0.7512 | 0.041 | 25/25 | 21 pts, 1.9e-07 | 220 |
| Optimized Fock sup. | 22.9980 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 22.9980 | -5.3e-09 | 22.9980 | 16.9816 | 0.6967 | 0.030 | 25/25 | 21 pts, 1.6e-05 | 220 |
| Optimized sqz-vac sup. | 27.2984 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 27.2984 | +1.8e-10 | 27.2984 | 20.5748 | 0.8129 | 0.030 | 25/25 | 19 pts, 2.2e-04 | 220 |
| Optimized sqz-coh sup. | 27.6049 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 27.6049 | -4.0e-11 | 27.6049 | 20.8022 | 0.8137 | 0.029 | 25/25 | 19 pts, 1.5e-03 | 220 |

*Optimized Fock sup.: loss_in state at eta=0.95 used.*

*Optimized sqz-vac sup.: loss_in state at eta=0.95 used.*

*Optimized sqz-coh sup.: loss_in state at eta=0.95 used.*

