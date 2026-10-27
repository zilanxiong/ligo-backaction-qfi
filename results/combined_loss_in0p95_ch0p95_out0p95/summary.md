# Loss combinations with back-action, ⟨n⟩ = 5, 10 Hz – 1 kHz

The paper has no multi-placement loss results; the check is its state-independent β relation
(supplement Eq. 16) plus the exact η_out scaling of loss after the signal: at κ → 0,
F = η_out × β(η_ch)² × F_input-loss(η_in·η_ch·η_out). Closed form (Eq. 17) for
squeezed vacuum and coherent light; for the rest the β relation applied to the package's own
input-loss value, next to the package channel evaluated directly at κ = 0 with both losses.

## η_in = 0.95, η_ch = 0.95, η_out = 0.95

| probe | κ→0 by β relation | source | direct κ=0 (package) | rel. diff | F at 1 kHz | F at 30 Hz (κ≈1) | F at 10 Hz (κ≈9.5) | 10 Hz / 1 kHz | converged | cross-check (points, worst rel.) | largest cutoff |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Squeezed vac. (p) | 20.3877 | closed form: eta_out × beta² × paper Eq. 17 at eta_in·eta_ch·eta_out | 20.3877 | -6.3e-08 | 20.3877 | 15.3780 | 0.6020 | 0.030 | 25/25 | 19 pts, 1.5e-03 | 220 |
| Coherent | 3.7040 | closed form: eta_out × beta² × paper Eq. 17 at eta_in·eta_ch·eta_out | 3.7040 | -8.6e-10 | 3.7040 | 3.5125 | 0.5698 | 0.154 | 25/25 | 21 pts, 3.5e-06 | 220 |
| Even cat | 7.3707 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 7.3707 | -1.1e-08 | 7.3707 | 4.9804 | 0.5698 | 0.077 | 25/25 | 21 pts, 3.4e-06 | 220 |
| Squeezed cat | 20.3289 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 20.3289 | +1.3e-10 | 20.3289 | 15.3392 | 0.6020 | 0.030 | 25/25 | 19 pts, 3.9e-04 | 220 |
| Fock | 13.3137 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 13.3137 | -1.3e-09 | 13.3137 | 10.9915 | 0.5859 | 0.044 | 25/25 | 21 pts, 1.3e-06 | 220 |
| Optimized Fock sup. | 16.6654 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 16.6654 | +7.4e-09 | 16.6654 | 12.1736 | 0.5426 | 0.033 | 25/25 | 21 pts, 2.6e-06 | 220 |
| Optimized sqz-vac sup. | 20.3593 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 20.3593 | +9.4e-11 | 20.3593 | 15.3749 | 0.6024 | 0.030 | 25/25 | 19 pts, 2.7e-04 | 220 |
| Optimized sqz-coh sup. | 20.3441 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 20.3441 | +2.8e-10 | 20.3440 | 15.1392 | 0.5998 | 0.029 | 25/25 | 19 pts, 5.4e-04 | 220 |

