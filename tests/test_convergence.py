"""Fock-cutoff convergence for the back-action channel.

Why these exist: the package default ``N_basis=20`` is fine for displacement
and rotation but silently wrong once radiation pressure shears the state.
Nothing in the package checks this, so the study layer must.
"""

import numpy as np
import pytest
import qutip as qt

import backaction_qfi as b
import quantum_sensing.dynamics as dynamics


# --- the ladder itself (no physics) ---------------------------------------

def test_converge_reports_success_and_failure():
    good = b.converge(lambda N: 1.0 + 1.0 / N**2, cutoffs=[20, 40, 80], rtol=1e-2)
    assert good.converged and good.cutoff == 40
    assert good.history == [(20, 1.0 + 1 / 400), (40, 1.0 + 1 / 1600)]

    bad = b.converge(lambda N: float(N), cutoffs=[20, 40, 80], rtol=1e-6)
    assert not bad.converged and bad.cutoff == 80 and bad.value == 80.0


def test_converge_can_demand_several_stable_rungs():
    # 1/N² steps: 20->40 and 40->80 both pass rtol=1e-2, so n_stable=2 lands on 80.
    res = b.converge(lambda N: 1.0 + 1.0 / N**2, cutoffs=[20, 40, 80], rtol=1e-2, n_stable=2)
    assert res.converged and res.cutoff == 80


def test_converge_needs_at_least_two_cutoffs():
    with pytest.raises(ValueError):
        b.converge(lambda N: 1.0, cutoffs=[20])


def test_suggested_cutoff_grows_with_kappa_and_nbar():
    assert b.suggested_cutoff(0.0, 1.0) < b.suggested_cutoff(0.0, 3.0)
    assert b.suggested_cutoff(1.0, 1.0) < b.suggested_cutoff(5.0, 1.0)
    assert b.suggested_cutoff(0.0, 0.0) >= 8


# --- truncation alarm ------------------------------------------------------

def test_tail_population_flags_truncation():
    assert b.tail_population(qt.coherent(12, 2.5)) > 1e-3     # badly truncated
    assert b.tail_population(qt.coherent(60, 2.5)) < 1e-10    # comfortably held
    cat = b.cat_even(40, 2.0)
    sheared = b.rp_channel(kappa_ba=3.0, N_basis=40, rho=cat)
    assert b.tail_population(sheared) > 10 * b.tail_population(cat)


# --- solver injection ------------------------------------------------------

def test_tight_solver_is_scoped_and_restored():
    assert dynamics.qt is qt
    with b.tight_solver() as opts:
        assert dynamics.qt is not qt
        assert opts["atol"] == b.SOLVER_OPTIONS["atol"]
        assert dynamics.qt.destroy(3) == qt.destroy(3)    # everything else passes through
    assert dynamics.qt is qt


def test_rp_channel_matches_package_at_default_tolerances():
    # Same physics as the raw package call; tight tolerances only sharpen it.
    # The two outputs differ by the default solver's own error (~1e-6 relative,
    # measured ~4e-6 here), which is the same order as the 1e-5 finite-difference
    # step used for the QFI: the reason this module injects tight tolerances.
    psi = b.cat_even(40, 2.0)
    ours = b.rp_channel(kappa_ba=1.0, N_basis=40, rho=psi)
    theirs = dynamics.get_state_single_mode_rp(kappa_ba=1.0, N_basis=40, rho=psi)
    assert (ours - theirs).norm() < 1e-4
    assert abs(ours.tr() - 1) < 1e-9


# --- probes ----------------------------------------------------------------

@pytest.mark.parametrize("nbar", [0.5, 2.0, 5.0])
def test_probes_hit_their_photon_number(nbar):
    N = 120
    assert b.mean_photon_number(b.coherent(N, nbar)) == pytest.approx(nbar, rel=1e-6)
    assert b.mean_photon_number(b.squeezed(N, nbar)) == pytest.approx(nbar, rel=1e-3)
    assert b.mean_photon_number(b.fock(N, 3)) == pytest.approx(3.0)


def test_squeezed_probe_anti_squeezes_x():
    # QFI for the epsilon_a signal is 8 Var(x), so the useful probe has large Var(x).
    r = np.arcsinh(np.sqrt(2.0))
    assert b.var_x(b.squeezed(120, 2.0)) == pytest.approx(np.exp(2 * r) / 2, rel=1e-3)
    assert b.var_x(b.squeezed(120, 2.0, squeeze_p=False)) == pytest.approx(np.exp(-2 * r) / 2, rel=1e-3)


# --- convergence on the real channel --------------------------------------

def test_vacuum_qfi_is_four_and_converges_immediately():
    res = b.converged_qfi(lambda N: b.fock(N, 0), cutoffs=[10, 15, 20], kappa_ba=1.0)
    assert res.converged and res.cutoff == 15
    assert res.value == pytest.approx(4.0, rel=1e-6)          # 8 Var(x) = 8 * 1/2


@pytest.mark.parametrize("kappa", [0.0, 1.0, 2.0])
def test_cat_converges_to_the_lossless_analytic_qfi(kappa):
    # Lossless: QFI(epsilon_a) = 8 Var(x) for a pure probe, independent of
    # kappa_ba because the shear generator x² commutes with the signal
    # generator x.  Agreement with the closed form is the real test that the
    # accepted cutoff is converged, not merely self-consistent.
    res = b.converged_qfi(lambda N: b.cat_even(N, 2.0),
                          cutoffs=[20, 30, 40, 60, 80], kappa_ba=kappa, rtol=1e-5)
    assert res.converged
    assert res.value == pytest.approx(8 * b.var_x(b.cat_even(res.cutoff, 2.0)), rel=1e-5)


def test_undersized_ladder_is_rejected():
    # A cat with alpha=2 under shear needs far more than 16 levels.
    res = b.converged_qfi(lambda N: b.cat_even(N, 2.0),
                          cutoffs=[12, 14, 16], kappa_ba=1.0, rtol=1e-3)
    assert not res.converged


def test_default_cutoff_is_not_trustworthy_under_backaction():
    # Regression guard for the trap: a p-squeezed probe at nbar=2, shear 1,
    # 10% detection loss.  N_basis=20 is fine for the input but the sheared
    # output has nbar ~ 4.5 and the QFI at 20 levels is >10% off.
    coarse = b.rp_qfi(b.squeezed(20, 2.0), kappa_ba=1.0, eta_out=0.9)
    res = b.converged_qfi(lambda N: b.squeezed(N, 2.0),
                          cutoffs=[20, 30, 40, 60, 80, 120], kappa_ba=1.0, eta_out=0.9,
                          rtol=1e-3)
    assert res.converged and res.cutoff >= 60
    assert abs(coarse - res.value) / res.value > 0.10
    # and the ladder is monotone: truncation over-estimates this QFI
    values = [v for _, v in res.history]
    assert all(a > c for a, c in zip(values, values[1:]))
