# Back-action only, ⟨n⟩ = 5, 10 Hz – 1 kHz

Channel: unitary, `H = √2 ε_a x + (κ/2) x²`, no loss, no phase noise. κ(Ω) from
`quantum_sensing.rp_kappa` (package LIGO defaults). QFI for `ε_a` at `ε_a = 0`.

## Result

The QFI is **frequency independent** for every probe: `x²` commutes with the signal
generator `x`, so `F = 8 Var(x)` of the input, exactly. The Fock track reproduces this
at every frequency where its cutoff ladder converged. The strain bound
`σ_h/h_SQL = √(2/(κF))` is frequency dependent through κ(Ω) only.

| probe | ⟨n⟩ | F exact (8 Var x) | Fock-track range over f | worst rel. err. | converged points | largest cutoff | σ_h/h_SQL at 30 Hz |
|---|---|---|---|---|---|---|---|
| squeezed vacuum | 5.000 | 87.818 | 87.818 – 87.818 | 2.9e-06 | 23/25 (2 beyond budget) | 219 | 0.155 |
| coherent | 5.000 | 4.000 | 4.000 – 4.001 | 1.6e-04 | 23/25 (2 beyond budget) | 219 | 0.726 |
| even cat | 5.000 | 84.004 | 84.004 – 84.004 | 7.2e-08 | 23/25 (2 beyond budget) | 219 | 0.158 |
| squeezed cat | 5.000 | 87.813 | 87.813 – 87.813 | 4.4e-06 | 23/25 (2 beyond budget) | 219 | 0.155 |
| Fock |5⟩ | 5.000 | 44.000 | 44.000 – 44.000 | 5.4e-09 | 23/25 (2 beyond budget) | 219 | 0.219 |
| optimised Fock superposition | 5.000 | 81.929 | 81.929 – 81.929 | 4.5e-09 | 23/25 (2 beyond budget) | 219 | 0.160 |
| optimised sqz-vac superposition | 5.004 | 87.821 | 87.815 – 87.821 | 6.2e-05 | 23/25 (2 beyond budget) | 219 | 0.155 |
| optimised sqz-coh superposition | 5.015 | 88.021 | 88.021 – 88.023 | 2.1e-05 | 23/25 (2 beyond budget) | 219 | 0.155 |

## Photon-number scaling (exact, lossless)

| probe | ⟨n⟩ values | F at ⟨n⟩ = 5 | F at largest ⟨n⟩ | local slope at top |
|---|---|---|---|---|
| squeezed vacuum | 1, 2, 3, 5, 8, 12, 20 | 87.82 | 327.88 | 0.97 |
| coherent | 1, 2, 3, 5, 8, 12, 20 | 4.00 | 4.00 | -0.00 |
| even cat | 1, 2, 3, 5, 8, 12, 20 | 84.00 | 324.00 | 0.98 |
| squeezed cat | 1, 2, 3, 5, 8, 12, 20 | 87.81 | 327.95 | 0.97 |
| Fock |5⟩ | 1, 2, 3, 5, 8, 12, 20 | 44.00 | 164.00 | 0.97 |
| optimised Fock superposition | 1, 2, 5, 10, 30, 35, 40 | 81.93 | 553.09 | 0.99 |
| optimised sqz-vac superposition | 0.5, 1, 2, 3, 5, 8, 10 | 87.82 | 169.02 | 0.95 |
| optimised sqz-coh superposition | 1, 2, 5, 10 | 88.02 | 168.30 | 0.94 |
