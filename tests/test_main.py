"""Tests for TopologicalInsulator-Lab. Every test grades a computed
invariant or spectrum against a closed form, a known spectrum, or an
independent referee inside main.py."""
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import main as ti


def test_jacobi_known_spectra():
    vals, _ = ti.jacobi_eigh([[0.0, 1.0], [1.0, 0.0]])
    assert abs(vals[0] + 1.0) < 1e-12 and abs(vals[1] - 1.0) < 1e-12
    tri = [[2.0, 1.0, 0.0], [1.0, 2.0, 1.0], [0.0, 1.0, 2.0]]
    vals, V = ti.jacobi_eigh(tri)
    expect = (2.0 - math.sqrt(2.0), 2.0, 2.0 + math.sqrt(2.0))
    assert all(abs(a - b) < 1e-12 for a, b in zip(vals, expect))
    assert ti.eig_residual_max(tri, vals, V) < 1e-12


def test_hermitian_known_spectra():
    assert all(abs(a - b) < 1e-12 for a, b in
               zip(ti.hermitian_eigenvalues(ti.PAULI_Y), (-1.0, 1.0)))
    known = [[1 + 0j, 1j], [-1j, 2 + 0j]]
    expect = ((3.0 - math.sqrt(5.0)) / 2.0, (3.0 + math.sqrt(5.0)) / 2.0)
    assert all(abs(a - b) < 1e-12
               for a, b in zip(ti.hermitian_eigenvalues(known), expect))


def test_winding_exact_integers():
    for v, w, want in ((2.0, 1.0, 0), (1.5, 1.0, 0), (1.1, 1.0, 0),
                       (0.9, 1.0, 1), (0.5, 1.0, 1), (0.2, 1.0, 1),
                       (1.0, 0.5, 0), (1.0, 2.0, 1)):
        assert abs(ti.ssh_winding(v, w) - want) < 1e-9


def test_zak_phase():
    assert abs(ti.ssh_zak_phase(0.5, 1.0) - math.pi) < 1e-9
    assert abs(ti.ssh_zak_phase(1.5, 1.0)) < 1e-9


def test_open_chain_topological_edge_modes():
    rep = ti.ssh_open_report(0.5, 1.0)
    assert rep["residual_max"] < 1e-10
    assert rep["n_in_gap_below_1pct"] == 2
    assert rep["min_abs_energy"] < 0.01 * rep["half_gap"]
    assert len(rep["in_gap_modes"]) == 2
    for m in rep["in_gap_modes"]:
        assert m["end_cell_weight"] > 0.5
        assert m["end_over_mid"] > 100.0


def test_open_chain_trivial_clean():
    rep = ti.ssh_open_report(1.5, 1.0)
    assert rep["n_in_gap_below_1pct"] == 0
    assert rep["in_gap_modes"] == []
    assert rep["min_abs_energy"] >= rep["half_gap"]


def test_chern_three_regimes():
    assert abs(ti.qwz_chern(3.0, 12)) < 1e-6
    assert abs(ti.qwz_chern(1.0, 12) - 1.0) < 1e-6
    assert abs(ti.qwz_chern(-1.0, 12) + 1.0) < 1e-6
    assert abs(ti.qwz_chern(1.5, 12) - 1.0) < 1e-6
    assert abs(ti.qwz_chern(-1.5, 12) + 1.0) < 1e-6


def test_chern_convergence():
    for u in ti.QWZ_MASSES:
        assert abs(ti.qwz_chern(u, 8) - ti.qwz_chern(u, 24)) < 1e-6


def test_qwz_eigenvector_solves_eigen_equation():
    assert ti._qwz_eigvec_residual_ok()


def test_ribbon_chern_phases_cross():
    for u in (1.0, -1.0):
        rep = ti.ribbon_report(u)
        assert rep["crossing_slices"] >= 3
        assert rep["min_abs_energy"] < 0.02 * rep["bulk_half_gap"]


def test_ribbon_trivial_stays_out_of_gap():
    rep = ti.ribbon_report(3.0)
    assert rep["crossing_slices"] == 0
    assert rep["in_gap_eigenvalues"] == 0
    assert rep["min_abs_energy"] > 0.9 * rep["bulk_half_gap"]


def test_ribbon_trace_invariants():
    H = ti.ribbon_hamiltonian(0.7, 1.0)
    ev = ti.hermitian_eigenvalues(H)
    n = len(H)
    assert abs(sum(ev) - sum(H[i][i].real for i in range(n))) < 1e-8
    assert abs(sum(e * e for e in ev)
               - sum(abs(H[i][j]) ** 2 for i in range(n)
                     for j in range(n))) < 1e-6


def test_classifier_scenarios():
    rep = ti.classifier_report()
    assert rep["ssh_trivial_v1.5_w1"]["knn"] == "trivial-insulator"
    assert rep["ssh_topological_v0.5_w1"]["knn"] == "ssh-topological"
    assert rep["ssh_critical_v1_w1"]["knn"] == "gapless-critical"
    assert rep["qwz_trivial_u3"]["knn"] == "trivial-insulator"
    assert rep["qwz_chern_plus_u1"]["knn"] == "chern-insulator"
    assert rep["qwz_chern_minus_u-1"]["knn"] == "chern-insulator"
    assert all(row["knn"] == row["policy"] for row in rep.values())


def test_classifier_grid_consistency():
    for feats, label in ti.TRAINING_GRID:
        assert ti.classify_features(feats) == label


def test_hold_gates():
    assert len(ti.HOLD_GATES) == 14
    assert ti.gate("real_material_claim") == "HOLD"
    assert ti.gate("topological_qubit_claim") == "HOLD"
    assert ti.gate("run_simulation") == "PROCEED"


def test_seats_never_override():
    brains = ti.brains_report()
    assert brains["decider"]["can_override_hold"] is False
    assert brains["decider"]["model_id"] == "gemini-argon-4"
    assert brains["advisor"]["can_override_hold"] is False
    assert brains["flow"]["rendered"] is False
    assert brains["flow"]["live"] == 0


def test_checks_all_pass():
    checks = ti.run_checks()
    assert all(ok for _, ok in checks), [n for n, ok in checks if not ok]


def test_cli_hold_demo_exit_2():
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent.parent
                             / "main.py"), "hold-demo"],
        capture_output=True, text=True)
    assert proc.returncode == 2
    assert "HOLD" in proc.stdout


def test_cli_checks_command():
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent.parent
                             / "main.py"), "checks"],
        capture_output=True, text=True)
    assert proc.returncode == 0
    assert "CHECKS 48/48" in proc.stdout


def test_demo_deterministic():
    main_py = str(Path(__file__).resolve().parent.parent / "main.py")
    a = subprocess.run([sys.executable, main_py, "demo"],
                       capture_output=True, text=True)
    b = subprocess.run([sys.executable, main_py, "demo"],
                       capture_output=True, text=True)
    assert a.stdout == b.stdout and a.returncode == 0
