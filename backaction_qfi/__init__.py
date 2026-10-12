"""Study layer for QFI under radiation-pressure back-action.

Everything here sits *outside* the group's ``quantum_sensing`` package, which
is vendored unmodified.  The package is driven through :mod:`.solver` (tight
ODE tolerances injected from outside), probe states come from :mod:`.probes`,
and every user-facing number goes through :mod:`.convergence`.
"""

from .solver import SOLVER_OPTIONS, tight_solver, rp_channel, rp_qfi
from .probes import (coherent, squeezed, cat_even, cat, squeezed_cat, fock,
                     mean_photon_number, var_x, FAMILIES, make_probe)
from .gaussian import gaussian_channel, gaussian_qfi, homodyne_noise_spectrum, fd_squeeze_angle, strain_noise_from_qfi
from .optimised import load_probe, quarter_turn, rotate
from .convergence import (
    ConvergenceResult, converge, converged_qfi, tail_population, suggested_cutoff,
)

__all__ = [
    "SOLVER_OPTIONS", "tight_solver", "rp_channel", "rp_qfi",
    "coherent", "squeezed", "cat_even", "cat", "squeezed_cat", "fock", "mean_photon_number", "var_x",
    "FAMILIES", "make_probe",
    "ConvergenceResult", "converge", "converged_qfi", "tail_population", "suggested_cutoff",
    "gaussian_channel", "gaussian_qfi", "homodyne_noise_spectrum", "fd_squeeze_angle", "strain_noise_from_qfi",
    "load_probe", "quarter_turn", "rotate",
]
