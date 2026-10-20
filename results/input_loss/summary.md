# Input loss + back-action, ⟨n⟩ = 5, 10 Hz – 1 kHz

Input loss acts before the sensing unitary, so the QFI must be frequency independent and
equal to its no-back-action value. Reference: the stored QFI of the paper's optimised states
(loss_in slice), closed forms for squeezed vacuum (paper Eq. 17) and coherent light (4), and the exact
no-back-action Fock value for the rest.

## η_in = 0.95

| probe | reference | source | Fock track over f (min – max) | worst rel. dev. | converged | tracks | largest cutoff |
|---|---|---|---|---|---|---|---|
| Squeezed vac. (p) | 42.8856 | closed form (paper Eq. 17) | 42.8829 – 42.8856 | 6.3e-05 | 25/25 | mesolve×23, branches×2 | 1120 |
| Coherent | 4.0000 | closed form: 4, independent of input loss | 4.0000 – 4.0000 | 7.0e-08 | 25/25 | mesolve×23, branches×2 | 1120 |
| Even cat | 31.9597 | no-back-action Fock value | 31.9597 – 31.9597 | 9.3e-08 | 25/25 | mesolve×23, branches×2 | 1120 |
| Squeezed cat | 42.7811 | no-back-action Fock value | 42.7810 – 42.7811 | 2.2e-06 | 25/25 | mesolve×23, branches×2 | 1120 |
| Fock | 27.8101 | no-back-action Fock value | 27.8101 – 27.8101 | 6.9e-09 | 25/25 | mesolve×23, branches×2 | 1120 |
| Optimized Fock sup. | 36.1954 | stored QFI (paper data, eta=0.95) | 36.1954 – 36.1954 | 2.5e-08 | 25/25 | mesolve×23, branches×2 | 1120 |
| Optimized sqz-vac sup. | 42.6405 | stored QFI (paper data, eta=0.95) | 42.7322 – 42.7398 | 2.3e-03 | 25/25 | mesolve×23, branches×2 | 1120 |
| Optimized sqz-coh sup. | 42.7022 | stored QFI (paper data, eta=0.95) | 43.3681 – 43.3688 | 1.6e-02 | 25/25 | mesolve×23, branches×2 | 1120 |

