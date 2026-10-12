"""Drive the package's channel with tight ODE tolerances, without editing it.

``quantum_sensing.dynamics`` calls ``qt.mesolve`` with QuTiP's default
tolerances (atol 1e-8, rtol 1e-6).  ``quantum_sensing.sld`` then differentiates
the channel by central differences with step 1e-5, so solver noise of order
1e-6 corrupts the derivative, and the back-action term ``g x^2`` is stiff at
large cutoffs.  The package exposes no way to pass solver options, so
:func:`tight_solver` temporarily replaces the ``qt`` name *inside the
dynamics module* with a proxy whose ``mesolve`` injects the options.  QuTiP
itself and the package source are never modified.
"""

import contextlib

import qutip

import quantum_sensing.dynamics as _dynamics
from quantum_sensing.sld import calculate_qfi
from quantum_sensing.dynamics import get_state_single_mode_rp

__all__ = ["SOLVER_OPTIONS", "tight_solver", "rp_channel", "rp_qfi"]

#: Tolerances used for every ``mesolve`` call made through this module.
SOLVER_OPTIONS = {"atol": 1e-12, "rtol": 1e-10, "nsteps": 100_000}


class _QutipProxy:
    """Looks like the ``qutip`` module, except ``mesolve`` carries fixed options."""

    def __init__(self, options):
        self._options = dict(options)

    def __getattr__(self, name):
        return getattr(qutip, name)

    def mesolve(self, H, rho0, tlist, c_ops=None, **kwargs):
        options = {**self._options, **(kwargs.pop("options", None) or {})}
        return qutip.mesolve(H, rho0, tlist, c_ops, options=options, **kwargs)


@contextlib.contextmanager
def tight_solver(options=None):
    """Within the block, the package's channels solve with tight tolerances.

    ``options`` overrides individual entries of :data:`SOLVER_OPTIONS`.
    Nested use is safe; the previous binding is always restored.
    """
    merged = {**SOLVER_OPTIONS, **(options or {})}
    previous = _dynamics.qt
    _dynamics.qt = _QutipProxy(merged)
    try:
        yield merged
    finally:
        _dynamics.qt = previous


def rp_channel(**kwargs):
    """``quantum_sensing.get_state_single_mode_rp`` solved with tight tolerances.

    Same keyword arguments as the package function (``kappa_ba``, ``rho``,
    ``N_basis``, ``epsilon_a``, ``eta_in``, ``eta_ch``, ``eta_out``, ...).
    """
    with tight_solver():
        return get_state_single_mode_rp(**kwargs)


def rp_qfi(rho, *, param_type="epsilon_a", param_value=0.0, **channel_kwargs):
    """QFI of the back-action channel's output for ``param_type`` (default the
    GW-frame signal ``epsilon_a``), for probe ``rho`` (ket or density matrix).

    ``N_basis`` is taken from ``rho`` unless given explicitly.  Returns a float.
    """
    channel_kwargs.setdefault("N_basis", rho.shape[0])
    return float(calculate_qfi(rp_channel, param_value=param_value,
                               param_type=param_type, rho=rho, **channel_kwargs))
