# Numerical checks on sweep points

## (a) finite-difference step and (b) exact derivative

| point | probe | h = 1e-4 | h = 1e-5 (default) | h = 1e-6 | exact | (h=1e-5 − exact)/exact | time FD / exact |
|---|---|---|---|---|---|---|---|
| channel loss + phase noise, 30 Hz | squeezed | 14.96237473 | 14.96237707 | 14.96237709 | 14.96237709 | -1.9e-09 | 3 s / 16 s |
| channel loss + phase noise, 30 Hz | opt_fock_sup | 25.38919874 | 25.38920294 | 25.38920298 | 25.38920301 | -3.1e-09 | 2 s / 1 s |
| channel loss + phase noise, 10 Hz | fock | 2.55694500 | 2.55694501 | 2.55694501 | 2.55694501 | -8.8e-11 | 4 s / 3 s |
| three losses, 10 Hz | squeezed | 0.60197987 | 0.60197987 | 0.60197987 | 0.60197987 | -3.9e-11 | 59 s / 34 s |
| three losses, 30 Hz | cat | 4.97975646 | 4.97975675 | 4.97975675 | 4.97975674 | +3.1e-10 | 2 s / 1 s |

## (c) extended cutoff ladder at 10 Hz of the phase-noise sweep (exact derivative)

| probe | N=160 | N=220 (sweep) | N=300 | N=400 | N=560 | (560 − 220)/220 |
|---|---|---|---|---|---|---|
| squeezed | 1.799085 | 1.799085 | 1.799085 | 1.799085 | 1.799085 | +4.9e-12 |
| fock | 2.556945 | 2.556945 | 2.556945 | 2.556945 | 2.556945 | +3.4e-14 |
| opt_fock_sup | 1.662326 | 1.662326 | 1.662326 | 1.662326 | 1.662326 | +1.7e-09 |

## (d) eigenvalue cutoff in the QFI sum

| point | probe | eig_tol = 1e-6 | 1e-8 (default) | 1e-10 | (1e-10 − 1e-8)/1e-8 |
|---|---|---|---|---|---|
| channel loss + phase noise, 30 Hz | squeezed | 14.96235555 | 14.96237709 | 14.96237733 | +1.6e-08 |
| channel loss + phase noise, 30 Hz | opt_fock_sup | 25.38917624 | 25.38920301 | 25.38920339 | +1.5e-08 |
| channel loss + phase noise, 10 Hz | fock | 2.55692890 | 2.55694501 | 2.55694515 | +5.5e-08 |
