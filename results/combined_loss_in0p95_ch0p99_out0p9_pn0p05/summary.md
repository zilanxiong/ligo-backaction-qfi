# Loss combinations with back-action, ⟨n⟩ = 5, 10 Hz – 1 kHz

The paper has no multi-placement loss results; the check is its state-independent β relation
(supplement Eq. 16) plus the exact η_out scaling of loss after the signal: at κ → 0,
F = η_out × β(η_ch)² × F_input-loss(η_in·η_ch·η_out). Closed form (Eq. 17) for
squeezed vacuum and coherent light; for the rest the β relation applied to the package's own
input-loss value, next to the package channel evaluated directly at κ = 0 with both losses.

## η_in = 0.95, η_ch = 0.99, η_out = 0.9, input phase noise 50 mrad

| probe | κ→0 by β relation | source | direct κ=0 (package) | rel. diff | F at 1 kHz | F at 30 Hz (κ≈1) | F at 10 Hz (κ≈9.5) | 10 Hz / 1 kHz | converged | cross-check (points, worst rel.) | largest cutoff |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Squeezed vac. (p) | 15.6603 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 15.6603 | +2.0e-10 | 15.6603 | 11.0933 | 0.3755 | 0.024 | 25/25 | 18 pts, 4.1e-05 | 220 |
| Coherent | 3.4369 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 3.4369 | -1.0e-08 | 3.4369 | 3.1837 | 0.3808 | 0.111 | 25/25 | 20 pts, 1.3e-08 | 220 |
| Even cat | 6.2440 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 6.2440 | +1.7e-08 | 6.2440 | 3.9761 | 0.3808 | 0.061 | 25/25 | 21 pts, 3.6e-06 | 220 |
| Squeezed cat | 15.7451 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 15.7451 | +6.2e-10 | 15.7451 | 11.1176 | 0.3755 | 0.024 | 25/25 | 18 pts, 5.8e-06 | 220 |
| Fock | 12.0761 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 12.0761 | -2.0e-09 | 12.0761 | 9.6434 | 0.4330 | 0.036 | 25/25 | 21 pts, 4.5e-09 | 220 |
| Optimized Fock sup. | 15.1220 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 15.1220 | +1.7e-10 | 15.1220 | 11.5815 | 0.4193 | 0.028 | 25/25 | 20 pts, 2.1e-09 | 220 |
| Optimized sqz-vac sup. | 15.7640 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 15.7640 | +2.2e-10 | 15.7640 | 11.1089 | 0.3755 | 0.024 | 25/25 | 18 pts, 2.0e-07 | 220 |
| Optimized sqz-coh sup. | 15.7665 | eta_out × beta relation (paper Eq. 16) on the package's input-loss value at eta_in·eta_ch·eta_out | 15.7665 | +1.8e-09 | 15.7665 | 11.1103 | 0.3756 | 0.024 | 25/25 | 18 pts, 1.6e-06 | 220 |

*Optimized Fock sup.: loss_ch state at eta=0.85 used (slice matched to total transmission).*

*Optimized sqz-vac sup.: loss_ch state at eta=0.85 used (slice matched to total transmission).*

*Optimized sqz-coh sup.: loss_ch state at eta=0.85 used (slice matched to total transmission).*

