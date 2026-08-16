"""Fock-sector weights and the correlated binding-energy contract."""
import numpy as np
import pytest

from dlcq.fock_weights import binding_series, dominant_sector, fock_weights


def test_weights_sum_to_one(python_K21):
    w = fock_weights(python_K21, state_idx=0)
    assert sum(w.values()) == pytest.approx(1.0, abs=1e-10)
    assert all(v > -1e-12 for v in w.values())


def test_sectors_are_odd_baryon_counts(python_K21):
    # N=3, B=1: parton counts are 3, 5, 7, ... — valence plus whole pairs.
    w = fock_weights(python_K21, state_idx=0)
    assert all(L % 2 == 1 and L >= 3 for L in w)


def test_ground_state_is_valence_dominated(python_K21):
    # At m/g = 1.6 the five-quark content is suppressed by 10^2-10^3
    # (Fig. 6 rescalings), so valence must dominate outright.
    L, weight = dominant_sector(python_K21, state_idx=0)
    assert L == 3
    assert weight > 0.9


def test_weights_agree_with_structure_function(python_K21):
    # integral (q + qbar) dx counts partons: it must equal sum_L L * w_L.
    from dlcq.observables import structure_function

    w = fock_weights(python_K21, state_idx=0)
    x, q, qbar = structure_function(python_K21, state_idx=0)
    dx = 2.0 / python_K21.K_code
    n_partons = float(np.sum(q + qbar) * dx)
    assert n_partons == pytest.approx(sum(L * v for L, v in w.items()),
                                      rel=1e-10)


def test_binding_series_difference():
    K, delta = binding_series([10, 12, 14], [5.0, 4.8, 4.7],
                              [10, 12, 14], [2.0, 1.9, 1.85], n_thresh=2)
    assert np.array_equal(K, [10, 12, 14])
    assert np.allclose(delta, [1.0, 1.0, 1.0])


def test_binding_series_rejects_mismatched_grids():
    # The N=3 hexaquark trap: B=2 grids are even, B=1 grids odd — a
    # correlated difference must refuse, not interpolate.
    with pytest.raises(ValueError, match="K grids differ"):
        binding_series([10, 12], [5.0, 4.8], [11, 13], [2.0, 1.9])


def test_binding_series_rejects_hamiltonian_mix():
    with pytest.raises(ValueError, match="Hamiltonian mismatch"):
        binding_series([10, 12], [5.0, 4.8], [10, 12], [2.0, 1.9],
                       hamiltonian_state="improved",
                       hamiltonian_thresh="standard")


def test_pair_correlation_counts_pairs(python_K21):
    # Total correlation mass equals the mean pair count 2 * <L(L-1)/2>.
    import numpy as np
    from dlcq.fock_weights import fock_weights, pair_correlation

    ks, C = pair_correlation(python_K21, 0)
    w = fock_weights(python_K21, 0)
    expected = sum(v * L * (L - 1) for L, v in w.items())
    assert np.isclose(C.sum(), expected, rtol=1e-10)


def test_pair_correlation_momentum_constraint(python_K21):
    # No pair can exceed the total momentum: C[k1,k2]=0 for k1+k2 > K_code.
    import numpy as np
    from dlcq.fock_weights import pair_correlation

    ks, C = pair_correlation(python_K21, 0)
    K = python_K21.K_code
    for i, k1 in enumerate(ks):
        for j, k2 in enumerate(ks):
            if k1 + k2 > K:
                assert C[i, j] == 0.0
