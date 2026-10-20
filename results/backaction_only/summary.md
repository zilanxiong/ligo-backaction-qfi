# Back-action only, ⟨n⟩ = 5, 10 Hz – 1 kHz

Channel: unitary, `H = √2 ε_a x + (κ/2) x²`, no loss, no phase noise. κ(Ω) from
`quantum_sensing.rp_kappa` (package LIGO defaults). QFI for `ε_a` at `ε_a = 0`.

## Result

The QFI is **frequency independent** for every probe: `x²` commutes with the signal
generator `x`, so `F = 8 Var(x)` of the input, exactly. The Fock track reproduces this
at every frequency.

| probe | ⟨n⟩ | F exact (8 Var x) | Fock-track range over f | worst rel. err. | converged points | largest cutoff |
|---|---|---|---|---|---|---|
| squeezed vacuum | 5.000 | 87.818 | 87.818 – 87.818 | 2.9e-06 | 25/25 (2 via pure-state track) | 560 |
| coherent | 5.000 | 4.000 | 4.000 – 4.001 | 1.6e-04 | 25/25 (2 via pure-state track) | 560 |
| even cat | 5.000 | 84.004 | 84.004 – 84.004 | 7.2e-08 | 25/25 (2 via pure-state track) | 560 |
| squeezed cat | 5.000 | 87.813 | 87.813 – 87.813 | 4.4e-06 | 25/25 (2 via pure-state track) | 560 |
| Fock |5⟩ | 5.000 | 44.000 | 44.000 – 44.000 | 5.4e-09 | 25/25 (2 via pure-state track) | 560 |
| optimised Fock superposition | 5.000 | 81.929 | 81.929 – 81.929 | 4.5e-09 | 25/25 (2 via pure-state track) | 560 |
| optimised sqz-vac superposition | 5.004 | 87.821 | 87.815 – 87.821 | 6.2e-05 | 25/25 (2 via pure-state track) | 560 |
| optimised sqz-coh superposition | 5.015 | 88.021 | 88.021 – 88.023 | 2.1e-05 | 25/25 (2 via pure-state track) | 560 |
