# Loss combinations with back-action, ⟨n⟩ = 5, 10 Hz – 1 kHz

The paper has no multi-placement loss results; the check is its state-independent β relation
(supplement Eq. 16) plus the exact η_out scaling of loss after the signal: at κ → 0,
F = η_out × β(η_ch)² × F_input-loss(η_in·η_ch·η_out). Closed form (Eq. 17) for
squeezed vacuum and coherent light; for the rest the β relation applied to the package's own
input-loss value, next to the package channel evaluated directly at κ = 0 with both losses.

## η_in = 1, η_ch = 0.95, η_out = 0.95

| probe | κ→0 by β relation | source | direct κ=0 (package) | rel. diff | F at 1 kHz | F at 30 Hz (κ≈1) | F at 10 Hz (κ≈9.5) | 10 Hz / 1 kHz | converged | cross-check (points, worst rel.) | largest cutoff |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Squeezed vac. (p) | 26.7228 | closed form: eta_out × beta² × paper Eq. 17 at eta_in·eta_ch·eta_out | 26.7228 | -5.6e-08 | 26.7228 | 18.7254 | 0.6061 | 0.023 | 25/25 | 19 pts, 1.9e-03 | 220 |
| Coherent | 3.7040 | closed form: eta_out × beta² × paper Eq. 17 at eta_in·eta_ch·eta_out | 3.7040 | -3.9e-10 | 3.7040 | 3.5125 | 0.5698 | 0.154 | 25/25 | 21 pts, 4.1e-06 | 220 |
| Even cat | 13.2173 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 13.2173 | +2.7e-10 | 13.2173 | 7.2047 | 0.5698 | 0.043 | 25/25 | 21 pts, 5.1e-06 | 220 |
| Squeezed cat | 26.6426 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 26.6426 | +6.9e-11 | 26.6426 | 18.6738 | 0.6061 | 0.023 | 25/25 | 19 pts, 5.4e-03 | 220 |
| Fock | 17.8566 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 17.8566 | +2.1e-10 | 17.8566 | 14.7437 | 0.8433 | 0.047 | 25/25 | 21 pts, 8.8e-07 | 220 |
| Optimized Fock sup. | 23.3728 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 23.3728 | -9.4e-09 | 23.3728 | 16.3701 | 0.5943 | 0.025 | 25/25 | 21 pts, 2.1e-05 | 220 |
| Optimized sqz-vac sup. | 26.7003 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 26.7003 | +7.2e-11 | 26.7003 | 18.7304 | 0.6065 | 0.023 | 25/25 | 19 pts, 5.3e-04 | 220 |
| Optimized sqz-coh sup. | 26.7708 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 26.7708 | -1.4e-09 | 26.7707 | 18.4671 | 0.6040 | 0.023 | 25/25 | 19 pts, 5.8e-04 | 220 |

