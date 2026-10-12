"""Probe states at a given Fock cutoff.

Every builder takes the cutoff ``N`` first, so ``lambda N: squeezed(N, 2.0)``
is a valid ``state_builder`` for :func:`backaction_qfi.convergence.converged_qfi`.

Quadrature convention (matches ``quantum_sensing.dynamics``):
``x = (a + a†)/√2``, ``p = -i(a - a†)/√2``, vacuum ``Var(x) = Var(p) = 1/2``.
The GW-frame signal ``epsilon_a (a + a†)`` has generator ``√2 x``, so for a
pure probe through a lossless channel the QFI is ``8 Var(x)``: the useful
squeezing for this signal is *p-squeezed* (``x`` anti-squeezed).
"""

import numpy as np
import qutip as qt

__all__ = ["coherent", "squeezed", "cat_even", "fock", "mean_photon_number", "var_x"]


def coherent(N, nbar):
    """Coherent state with mean photon number ``nbar`` (real amplitude)."""
    return qt.coherent(N, np.sqrt(nbar))


def squeezed(N, nbar, squeeze_p=True):
    """Squeezed vacuum with mean photon number ``nbar = sinh²r``.

    ``squeeze_p=True`` (default) squeezes ``p`` and anti-squeezes ``x``: the
    optimal fixed-``nbar`` Gaussian probe for the ``epsilon_a`` signal.
    """
    r = np.arcsinh(np.sqrt(nbar))
    # qutip.squeeze(N, z) = exp((z* a² - z a†²)/2): real z > 0 squeezes x.
    z = -r if squeeze_p else r
    return (qt.squeeze(N, z) * qt.fock(N, 0)).unit()


def cat_even(N, alpha):
    """Even cat ``|α⟩ + |-α⟩`` (normalised); mean photon number ``|α|² tanh|α|²``."""
    return (qt.coherent(N, alpha) + qt.coherent(N, -alpha)).unit()


def fock(N, n):
    """Number state ``|n⟩``."""
    return qt.fock(N, n)


def mean_photon_number(state):
    N = state.shape[0]
    return float(np.real(qt.expect(qt.num(N), state)))


def var_x(state):
    """Variance of ``x = (a + a†)/√2`` (vacuum: 1/2)."""
    N = state.shape[0]
    a = qt.destroy(N)
    x = (a + a.dag()) / np.sqrt(2)
    return float(np.real(qt.expect(x * x, state) - qt.expect(x, state) ** 2))
