# Loss combinations with back-action, ⟨n⟩ = 5, 10 Hz – 1 kHz

The paper has no multi-placement loss results; the check is its state-independent β relation
(supplement Eq. 16) plus the exact η_out scaling of loss after the signal: at κ → 0,
F = η_out × β(η_ch)² × F_input-loss(η_in·η_ch·η_out). Closed form (Eq. 17) for
squeezed vacuum and coherent light; for the rest the β relation applied to the package's own
input-loss value, next to the package channel evaluated directly at κ = 0 with both losses.

## η_in = 0.95, η_ch = 0.99, η_out = 0.9, input phase noise 200 mrad

| probe | κ→0 by β relation | source | direct κ=0 (package) | rel. diff | F at 1 kHz | F at 30 Hz (κ≈1) | F at 10 Hz (κ≈9.5) | 10 Hz / 1 kHz | converged | cross-check (points, worst rel.) | largest cutoff |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Squeezed vac. (p) | 7.9604 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 7.9604 | +6.5e-10 | 7.9604 | 6.0467 | 0.3534 | 0.044 | 25/25 | 19 pts, 8.8e-05 | 220 |
| Coherent | 2.1752 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 2.1752 | -1.7e-08 | 2.1752 | 2.0709 | 0.3573 | 0.164 | 25/25 | 20 pts, 2.0e-08 | 220 |
| Even cat | 4.9049 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 4.9049 | +2.2e-08 | 4.9049 | 2.9595 | 0.3573 | 0.073 | 25/25 | 21 pts, 6.9e-06 | 220 |
| Squeezed cat | 7.8979 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 7.8979 | +8.4e-10 | 7.8979 | 5.9761 | 0.3531 | 0.045 | 25/25 | 19 pts, 2.5e-05 | 220 |
| Fock | 12.0761 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 12.0761 | -2.0e-09 | 12.0761 | 9.6434 | 0.4330 | 0.036 | 25/25 | 21 pts, 4.5e-09 | 220 |
| Optimized Fock sup. | 12.6005 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 12.6005 | -2.8e-09 | 12.6005 | 10.2367 | 0.4353 | 0.035 | 25/25 | 21 pts, 6.2e-09 | 220 |
| Optimized sqz-vac sup. | 12.3096 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 12.3096 | -4.1e-10 | 12.3096 | 9.4921 | 0.3036 | 0.025 | 25/25 | 21 pts, 6.8e-07 | 220 |
| Optimized sqz-coh sup. | 12.5768 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 12.5768 | -2.4e-08 | 12.5768 | 9.3263 | 0.3798 | 0.030 | 25/25 | 21 pts, 1.6e-07 | 220 |

*Optimized Fock sup.: loss_ch state at eta=0.85 used (slice matched to total transmission).*

*Optimized sqz-vac sup.: loss_ch state at eta=0.85 used (slice matched to total transmission).*

*Optimized sqz-coh sup.: loss_ch state at eta=0.85 used (slice matched to total transmission).*

