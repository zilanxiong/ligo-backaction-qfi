"""Exact Fock-track QFI for *input noise followed by the lossless sensing stage*.

Input loss/dephasing makes the probe mixed, so the pure-state shortcut no longer
applies; but the sensing stage itself (shear + signal, no channel loss) is a
unitary ``U``.  A mixed state ``ρ = Σ λ_i |φ_i⟩⟨φ_i|`` therefore maps to
``Σ λ_i U|φ_i⟩⟨φ_i|U†``: the spectrum is unchanged and only the eigenvectors
evolve.  So the noise stage is run through the *package's own channel* at a
modest cutoff ``N_noise`` (density matrix, cheap), the state is diagonalised,
each eigen-branch is embedded in a large cutoff ``N_unitary`` and propagated
through ``U`` as a ket, and the QFI is assembled exactly::

    F = Σ_{a,b ∈ S} 2 |⟨ψ_a|∂ρ|ψ_b⟩|² / (λ_a + λ_b)  +  4 Σ_{a ∈ S} λ_a ‖P_⊥ ∂ψ_a‖²

where ``S`` is the support of ρ, ``ψ_a = U φ_a``, ``∂ψ_a`` the central-difference
derivative in ``epsilon_a``, and ``P_⊥`` projects off the support (the second term
is the contribution of pairs with one index in the kernel of ρ).

Valid only when stage 2 carries no loss or dephasing and stage 3 is off.
Used to reach the low-frequency points (kappa_ba ~ 6-10) where the
density-matrix track would need a cutoff of ~1000.
"""

import numpy as np
import qutip

from quantum_sensing.conversions import ba_to_g
from quantum_sensing.dynamics import get_state_single_mode

from .solver import tight_solver

__all__ = ["rp_qfi_input_noise_branches"]


def _unitaries(N, kappa_ba, prec):
    """U(ε) = exp(-i[(κ/2)x² + ε(a+a†)]) for ε = -h, 0, +h, by exact diagonalisation."""
    a = qutip.destroy(N).full()
    x2 = (a + a.conj().T) @ (a + a.conj().T) / 2
    H0 = ba_to_g(kappa_ba, 1.0) * x2
    D = a + a.conj().T
    out = []
    for eps in (-prec, 0.0, prec):
        H = H0 + eps * D
        w, v = np.linalg.eigh(H)
        out.append((v * np.exp(-1j * w)) @ v.conj().T)
    return out


def rp_qfi_input_noise_branches(psi, *, kappa_ba=0.0, eta_in=1.0, pn_in=0.0, N_unitary=None,
                                N_noise=None, prec=1e-5, lam_min=1e-13, unitaries=None):
    """QFI for ``epsilon_a`` of: input noise (``eta_in``, ``pn_in``) → shear ``kappa_ba`` + signal.

    ``psi``: probe ket at cutoff ``N_noise`` (defaults to its own size).
    ``N_unitary``: cutoff for the sensing unitary (default ``4 * N_noise``).
    ``unitaries``: optional precomputed ``_unitaries(N_unitary, kappa_ba, prec)``
    so several probes can share them.
    """
    N0 = int(N_noise or psi.shape[0])
    N = int(N_unitary or 4 * N0)
    if N < N0:
        raise ValueError("N_unitary must be at least the noise-stage cutoff")
    psi0 = qutip.Qobj(psi.full().ravel()[:N0]).unit() if psi.shape[0] != N0 else psi

    # stage 1 through the package: noise only (stage 2 is the identity here)
    with tight_solver():
        rho = get_state_single_mode(rho=psi0, N_basis=N0, eta_in=eta_in, pn_in=pn_in)
    lam, vec = np.linalg.eigh(rho.full())
    keep = lam > lam_min
    lam, vec = lam[keep], vec[:, keep]
    phi = np.zeros((N, lam.size), dtype=complex)
    phi[:N0, :] = vec                                   # embed branches in the big cutoff

    Um, U0, Up = unitaries if unitaries is not None else _unitaries(N, kappa_ba, prec)
    psi_out = U0 @ phi                                  # eigenvectors of ρ_out (orthonormal)
    dpsi = (Up @ phi - Um @ phi) / (2 * prec)           # their derivatives

    # ∂ρ = Σ λ_i (|∂ψ_i⟩⟨ψ_i| + |ψ_i⟩⟨∂ψ_i|); matrix elements on the support
    G = psi_out.conj().T @ dpsi                         # G[a,i] = ⟨ψ_a|∂ψ_i⟩
    drho_S = (G * lam[None, :]) + (G * lam[None, :]).conj().T   # ⟨ψ_a|∂ρ|ψ_b⟩ = λ_b G[a,b] + λ_a conj(G[b,a])
    denom = lam[:, None] + lam[None, :]
    F_support = float(np.sum(2 * np.abs(drho_S) ** 2 / denom))

    # pairs with one index off the support: 4 Σ_a λ_a ‖P_⊥ ∂ψ_a‖²
    perp = dpsi - psi_out @ G
    F_kernel = float(4 * np.sum(lam * np.sum(np.abs(perp) ** 2, axis=0)))
    return F_support + F_kernel
