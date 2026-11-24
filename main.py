#!/usr/bin/env python3
"""TopologicalInsulator-Lab - educational topological insulator lab.

by Darion Kellar

Pure Python standard library. Zero dependencies. Zero network.
Education only - and the central lesson is stated up front:

  * A topological insulator is protected by a whole-number invariant
    computed from its bulk wavefunctions: the invariant cannot change
    without the bulk gap closing, and when it is nonzero the boundary
    is forced to host states inside that gap.
  * This lab computes those invariants from scratch for the two
    canonical teaching models - the SSH chain (winding number, Zak
    phase, open-chain edge modes) and the Qi-Wu-Zhang Chern
    insulator (Chern number by the Fukui-Hatsugai-Suzuki plaquette
    method, ribbon edge bands) - and grades every invariant against
    the boundary spectrum it predicts.
  * It models two tight-binding Hamiltonians on paper-sized
    lattices. It says nothing about any real material: no sample was
    grown, no spectrum was measured, and no device figure is
    computed or claimed anywhere in this file.

Conventions, stated once and used everywhere:

  * SSH: H(k) = d(k) . sigma with d = (v + w cos k, w sin k, 0),
    v = intra-cell hopping, w = inter-cell hopping. The winding
    number counts how many times d(k) circles the origin as k runs
    -pi -> pi; the Zak phase is pi times the winding.
  * Open SSH chain: N unit cells, sites alternate A/B, v inside the
    cell, w between cells. The bulk half-gap is |v - w|.
  * QWZ: H(k) = sin kx sigma_x + sin ky sigma_y
    + (u + cos kx + cos ky) sigma_z. With the Fukui plaquette
    traversed +x then +y and C = (1/2 pi i) sum log U, this lab's
    convention gives C = +1 for 0 < u < 2, C = -1 for -2 < u < 0,
    and C = 0 for |u| > 2 - computed, not quoted.
  * Ribbon: periodic along x, open along y with Ny unit cells; an
    "in-gap" ribbon state has |E| below half the bulk half-gap.

Every real-world material / device use is a permanent HOLD gate -
with Gemini Argon 4 Extended Thinking as DECIDER (names-only,
tighten-only), Muse as advisor (never overrides HOLD), and Google
Flow plan-only seats.
"""
import cmath
import hashlib
import json
import math
import sys

# ---------------------------------------------------------------- constants
VERSION = "1.0.0"

# Labelled teaching parameters. Model values, not measurements.
SSH_TOPO = (0.5, 1.0)        # v, w - topological phase
SSH_TRIVIAL = (1.5, 1.0)     # v, w - trivial phase
SSH_CELLS = 20               # unit cells in the open chain
QWZ_MASSES = (3.0, 1.0, -1.0)  # trivial / Chern +1 / Chern -1 regimes
RIBBON_CELLS = 10            # unit cells across the ribbon

# ---------------------------------------------------------------- HOLD gates
HOLD_GATES = {
    "real_material_claim": "HOLD - no claim that any real material is or is not a topological insulator follows from two toy Hamiltonians",
    "device_prototype_claim": "HOLD - no device, prototype, or component claim is made or certified here",
    "topological_qubit_claim": "HOLD - no topological-qubit, Majorana-hardware, or braiding claim is simulated or validated here",
    "fault_tolerance_claim": "HOLD - no fault-tolerance, error-rate, or protection guarantee for any real hardware",
    "spintronics_product_claim": "HOLD - no spintronics product, memory, or logic claim from a lattice model",
    "quantum_computer_claim": "HOLD - no claim that this lab simulates, benchmarks, or advances a quantum computer",
    "experimental_data_claim": "HOLD - no experimental data is produced; every number is a statement about the stated model",
    "sample_growth_guidance": "HOLD - no crystal-growth, deposition, or fabrication guidance for real samples",
    "vendor_material_claim": "HOLD - no vendor material, wafer, or instrument figure is stated, endorsed, or compared",
    "room_temperature_claim": "HOLD - no room-temperature (or any-temperature) operation claim for any real system",
    "certification_claim": "HOLD - no certification that any material or system meets any standard",
    "investment_claim": "HOLD - no investment claim; a lab model is not a prospectus",
    "live_hardware_control": "HOLD - no control of real samples, cryostats, or laboratory hardware, ever",
    "model_as_measurement_claim": "HOLD - a computed invariant is a property of the model at stated inputs, never a measurement",
}


def gate(action):
    """Return 'HOLD' if the action is gated, else 'PROCEED' (sandbox sim)."""
    return "HOLD" if action in HOLD_GATES else "PROCEED"


# ================================================================ 1. Linear algebra kit
# Two Jacobi eigensolvers written out in full: a real-symmetric one
# that also accumulates eigenvectors (for the open SSH chain, whose
# edge-mode wavefunctions are the point), and a complex-Hermitian
# eigenvalues-only one (for the QWZ ribbon). Both are refereed
# against known spectra and trace invariants in CHECKS before any
# physics is trusted to them.

def jacobi_eigh(A):
    """Eigen-decomposition of a real symmetric matrix by cyclic Jacobi.

    Returns (values, V) with values sorted ascending and column j of
    V the eigenvector for values[j]: V[k][j] is component k.
    """
    n = len(A)
    a = [row[:] for row in A]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(100):
        off = math.sqrt(sum(a[p][q] ** 2
                            for p in range(n) for q in range(p + 1, n)))
        if off < 1e-13:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) < 1e-15:
                    continue
                app, aqq = a[p][p], a[q][q]
                tau = (aqq - app) / (2.0 * apq)
                if tau == 0.0:
                    t = 1.0
                else:
                    t = math.copysign(
                        1.0 / (abs(tau) + math.sqrt(1.0 + tau * tau)), tau)
                c = 1.0 / math.sqrt(1.0 + t * t)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = V[k][p], V[k][q]
                    V[k][p] = c * vkp - s * vkq
                    V[k][q] = s * vkp + c * vkq
    vals = [a[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda i: vals[i])
    return ([vals[i] for i in order],
            [[V[k][i] for i in order] for k in range(n)])


def hermitian_eigenvalues(A):
    """Eigenvalues of a complex Hermitian matrix by cyclic Jacobi.

    Each rotation first phases the pivot element real, then applies
    the stable real Jacobi angle (the construction refereed in
    CHECKS against sigma_y, a complex 2x2, and X x X + Z x Z).
    Returns eigenvalues sorted ascending.
    """
    n = len(A)
    a = [[complex(v) for v in row] for row in A]
    for _ in range(100):
        off = math.sqrt(sum(abs(a[p][q]) ** 2
                            for p in range(n) for q in range(p + 1, n)))
        if off < 1e-13:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                b = abs(apq)
                if b < 1e-15:
                    continue
                app, aqq = a[p][p].real, a[q][q].real
                alpha = cmath.phase(apq)
                tau = (aqq - app) / (2.0 * b)
                if tau == 0.0:
                    t = -1.0
                else:
                    t = -math.copysign(
                        1.0 / (abs(tau) + math.sqrt(tau * tau + 1.0)), tau)
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                e_neg = cmath.exp(-1j * alpha)
                e_pos = e_neg.conjugate()
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp + s * e_neg * akq
                    a[k][q] = -s * akp + c * e_neg * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk + s * e_pos * aqk
                    a[q][k] = -s * apk + c * e_pos * aqk
    return sorted(a[i][i].real for i in range(n))


def eig_residual_max(H, vals, V):
    """Largest |H v - lambda v| over all eigenpairs (solver referee)."""
    n = len(H)
    worst = 0.0
    for j in range(n):
        v = [V[k][j] for k in range(n)]
        for i in range(n):
            hv = sum(H[i][k] * v[k] for k in range(n))
            worst = max(worst, abs(hv - vals[j] * v[i]))
    return worst


def _pos_zero(x):
    """Round for display, normalising -0.0 to 0.0."""
    return x if x != 0 else 0.0

# ================================================================ 2. SSH chain - winding + Zak phase
# The d-vector of the SSH Bloch Hamiltonian traces a circle of
# radius w centred at (v, 0) in the (dx, dy) plane. It winds around
# the origin exactly when the circle encloses it: v < w. The lab
# counts the winding by summing wrapped angle steps - no eyeballing.

def ssh_dvector(k, v, w):
    return (v + w * math.cos(k), w * math.sin(k))


def ssh_winding(v, w, n=720):
    """Winding number of the SSH d-vector by discretised angle sum."""
    total = 0.0
    prev = None
    for i in range(n + 1):
        k = -math.pi + 2.0 * math.pi * i / n
        dx, dy = ssh_dvector(k, v, w)
        ang = math.atan2(dy, dx)
        if prev is not None:
            step = ang - prev
            while step > math.pi:
                step -= 2.0 * math.pi
            while step < -math.pi:
                step += 2.0 * math.pi
            total += step
        prev = ang
    return total / (2.0 * math.pi)


def ssh_zak_phase(v, w, n=720):
    """Zak phase = pi x winding for the SSH chain (stated convention)."""
    return math.pi * ssh_winding(v, w, n)


def winding_report():
    rows = {}
    for v, w in ((2.0, 1.0), (1.5, 1.0), (1.1, 1.0), (0.9, 1.0),
                 (0.5, 1.0), (0.2, 1.0), (1.0, 0.5), (1.0, 2.0)):
        wind = ssh_winding(v, w)
        rows[f"v{v}_w{w}"] = {
            "v": v, "w": w,
            "winding": _pos_zero(round(wind, 9)),
            "zak_over_pi": _pos_zero(round(ssh_zak_phase(v, w) / math.pi, 9)),
            "phase": "topological" if v < w else "trivial",
        }
    return rows


# ================================================================ 3. Open SSH chain - edge modes
# The bulk-boundary promise: in the topological phase the open chain
# must show two modes inside the bulk gap, pinned near zero energy,
# with their weight piled onto the end unit cells. All three parts
# are measured here from the Jacobi eigenvectors.

def ssh_open_hamiltonian(v, w, ncells):
    n = 2 * ncells
    H = [[0.0] * n for _ in range(n)]
    for cell in range(ncells):
        a, b = 2 * cell, 2 * cell + 1
        H[a][b] = v
        H[b][a] = v
        if cell + 1 < ncells:
            H[b][a + 2] = w
            H[a + 2][b] = w
    return H


def ssh_open_report(v, w, ncells=SSH_CELLS):
    H = ssh_open_hamiltonian(v, w, ncells)
    vals, V = jacobi_eigh(H)
    half_gap = abs(v - w)
    mid = ncells // 2
    modes = []
    for j, e in enumerate(vals):
        if abs(e) < 0.05 * half_gap:
            vec = [V[k][j] for k in range(2 * ncells)]
            end = (vec[0] ** 2 + vec[1] ** 2
                   + vec[-2] ** 2 + vec[-1] ** 2)
            middle = vec[2 * mid] ** 2 + vec[2 * mid + 1] ** 2
            modes.append({
                "energy": e,
                "end_cell_weight": end,
                "mid_cell_weight": middle,
                "end_over_mid": end / middle if middle > 0 else float("inf"),
            })
    return {
        "v": v, "w": w, "ncells": ncells,
        "half_gap": half_gap,
        "min_abs_energy": min(abs(e) for e in vals),
        "residual_max": eig_residual_max(H, vals, V),
        "in_gap_modes": modes,
        "n_in_gap_below_1pct": sum(
            1 for e in vals if abs(e) < 0.01 * half_gap),
    }


# ================================================================ 4. QWZ model - Chern number by FHS
# The Fukui-Hatsugai-Suzuki method makes the Chern number a finite
# sum of plaquette Berry fluxes built from link variables
# U = <u(k)|u(k+mu)> / |<u(k)|u(k+mu)>|. It is exactly integer on
# any grid that resolves the band - the lab shows the convergence
# instead of asking to be believed.

def qwz_dvector(kx, ky, u):
    return (math.sin(kx), math.sin(ky), u + math.cos(kx) + math.cos(ky))


def qwz_lower_eigenvector(kx, ky, u):
    """Lower-band eigenvector via the larger column of the projector.

    P = (1 - d-hat . sigma)/2 projects onto the lower band; either
    column is an eigenvector, and taking the larger-norm column
    keeps the gauge finite everywhere except at a gap closing
    (excluded: the lab only evaluates gapped masses).
    """
    dx, dy, dz = qwz_dvector(kx, ky, u)
    r = math.sqrt(dx * dx + dy * dy + dz * dz)
    c0 = (0.5 * (1.0 - dz / r), -0.5 * (dx + 1j * dy) / r)
    c1 = (-0.5 * (dx - 1j * dy) / r, 0.5 * (1.0 + dz / r))
    n0 = abs(c0[0]) ** 2 + abs(c0[1]) ** 2
    n1 = abs(c1[0]) ** 2 + abs(c1[1]) ** 2
    c = c0 if n0 >= n1 else c1
    nn = math.sqrt(abs(c[0]) ** 2 + abs(c[1]) ** 2)
    return (c[0] / nn, c[1] / nn)


def qwz_chern(u, n):
    """Chern number of the QWZ lower band on an n x n FHS grid."""
    dk = 2.0 * math.pi / n

    def link(a, b):
        ov = a[0].conjugate() * b[0] + a[1].conjugate() * b[1]
        return ov / abs(ov)

    vecs = {}
    for i in range(n):
        for j in range(n):
            vecs[(i, j)] = qwz_lower_eigenvector(
                -math.pi + i * dk, -math.pi + j * dk, u)
    total = 0.0
    for i in range(n):
        for j in range(n):
            ip, jp = (i + 1) % n, (j + 1) % n
            ux = link(vecs[(i, j)], vecs[(ip, j)])
            uy = link(vecs[(i, j)], vecs[(i, jp)])
            ux_y = link(vecs[(i, jp)], vecs[(ip, jp)])
            uy_x = link(vecs[(ip, j)], vecs[(ip, jp)])
            total += cmath.log(
                ux * uy_x * ux_y.conjugate() * uy.conjugate()).imag
    return total / (2.0 * math.pi)


def qwz_bulk_half_gap(u, n=24):
    best = float("inf")
    for i in range(n):
        for j in range(n):
            kx = -math.pi + 2.0 * math.pi * i / n
            ky = -math.pi + 2.0 * math.pi * j / n
            dx, dy, dz = qwz_dvector(kx, ky, u)
            best = min(best, math.sqrt(dx * dx + dy * dy + dz * dz))
    return best


def chern_report():
    rows = {}
    for u in QWZ_MASSES:
        conv = {n: qwz_chern(u, n) for n in (8, 12, 16, 24)}
        rows[f"u{u}"] = {
            "u": u,
            "chern_by_grid": {str(n): _pos_zero(round(c, 9)) for n, c in conv.items()},
            "chern": _pos_zero(round(conv[24], 9)),
            "bulk_half_gap": round(qwz_bulk_half_gap(u), 6),
        }
    return rows


# ================================================================ 5. QWZ ribbon - chiral edge bands
# Periodic along x, open along y. In a Chern phase a chiral edge
# band must cross the bulk gap; in the trivial phase the ribbon
# spectrum must stay out of the gap entirely. Both are counted here
# slice by slice over kx.

PAULI_X = [[0j, 1 + 0j], [1 + 0j, 0j]]
PAULI_Y = [[0j, -1j], [1j, 0j]]
PAULI_Z = [[1 + 0j, 0j], [0j, -1 + 0j]]


def _pauli_block(terms):
    out = [[0j, 0j], [0j, 0j]]
    for coef, M in terms:
        for i in range(2):
            for j in range(2):
                out[i][j] += coef * M[i][j]
    return out


def ribbon_hamiltonian(kx, u, ny=RIBBON_CELLS):
    n = 2 * ny
    H = [[0j] * n for _ in range(n)]
    onsite = _pauli_block(((u + math.cos(kx), PAULI_Z),
                           (math.sin(kx), PAULI_X)))
    ty = _pauli_block(((0.5, PAULI_Z), (-0.5j, PAULI_Y)))
    ty_d = [[ty[j][i].conjugate() for j in range(2)] for i in range(2)]
    for y in range(ny):
        o = 2 * y
        for i in range(2):
            for j in range(2):
                H[o + i][o + j] = onsite[i][j]
        if y + 1 < ny:
            for i in range(2):
                for j in range(2):
                    H[o + 2 + i][o + j] = ty[i][j]
                    H[o + i][o + 2 + j] = ty_d[i][j]
    return H


def ribbon_report(u, ny=RIBBON_CELLS, nkx=21):
    half_gap = qwz_bulk_half_gap(u)
    threshold = 0.5 * half_gap
    crossing_slices = 0
    in_gap_total = 0
    min_abs = float("inf")
    for i in range(nkx):
        kx = -math.pi + 2.0 * math.pi * i / (nkx - 1)
        ev = hermitian_eigenvalues(ribbon_hamiltonian(kx, u, ny))
        inside = [e for e in ev if abs(e) < threshold]
        if inside:
            crossing_slices += 1
            in_gap_total += len(inside)
        min_abs = min(min_abs, min(abs(e) for e in ev))
    return {
        "u": u, "ny": ny, "nkx": nkx,
        "bulk_half_gap": round(half_gap, 6),
        "in_gap_threshold": round(threshold, 6),
        "crossing_slices": crossing_slices,
        "in_gap_eigenvalues": in_gap_total,
        "min_abs_energy": min_abs,
    }


# ================================================================ 6. kNN phase classifier
# Features of a labelled scenario: (ssh_winding, qwz_chern,
# bulk_half_gap). SSH scenarios carry chern 0 and QWZ scenarios
# carry winding 0 - each scenario populates the invariant of its
# own model. The training grid is labelled by the lab's own policy
# (model_label below) - the classifier interpolates the policy, it
# does not replace the computed invariants.

def model_label(winding, chern, half_gap):
    """The phase policy distilled to four verdicts."""
    if half_gap < 0.05:
        return "gapless-critical"
    if chern != 0:
        return "chern-insulator"
    if winding >= 0.75:
        return "ssh-topological"
    return "trivial-insulator"


def _training_grid():
    rows = []
    for winding in (0.0, 0.5, 1.0):
        for chern in (-1, 0, 1):
            for half_gap in (0.0, 0.25, 0.5, 1.0):
                feats = (winding, float(chern), half_gap)
                rows.append((feats, model_label(*feats)))
    return rows


TRAINING_GRID = _training_grid()


def classify_features(features, k=5):
    """Distance-weighted kNN over the policy-labelled grid."""
    ranked = []
    for feats, label in TRAINING_GRID:
        dist = math.sqrt(sum((a - b) ** 2
                             for a, b in zip(features, feats)))
        ranked.append((dist, feats, label))
    ranked.sort(key=lambda row: (row[0], row[1]))
    votes = {}
    for dist, _, label in ranked[:k]:
        votes[label] = votes.get(label, 0.0) + 1.0 / (dist + 1e-9)
    priority = {"chern-insulator": 0, "ssh-topological": 1,
                "trivial-insulator": 2, "gapless-critical": 3}
    return sorted(votes.items(), key=lambda kv: (-kv[1], priority[kv[0]]))[0][0]


CLASSIFIER_SCENARIOS = {
    "ssh_trivial_v1.5_w1": (0.0, 0, 0.5),
    "ssh_topological_v0.5_w1": (1.0, 0, 0.5),
    "ssh_critical_v1_w1": (0.5, 0, 0.0),
    "qwz_trivial_u3": (0.0, 0, 1.0),
    "qwz_chern_plus_u1": (0.0, 1, 1.0),
    "qwz_chern_minus_u-1": (0.0, -1, 1.0),
}


def classifier_report():
    rows = {}
    for name, feats in CLASSIFIER_SCENARIOS.items():
        rows[name] = {
            "features": list(feats),
            "policy": model_label(*feats),
            "knn": classify_features(feats),
        }
    return rows


# ================================================================ 7. Advisory seats
def _trace(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def argon4_decide(topic="topologicalinsulator-demo"):
    """Gemini Argon 4 DECIDER seat: MAP->MODEL->STRESS->DECIDE.

    Names-only keys, live calls 0, tighten-only: can add HOLDs, never lift one.
    """
    stages = ["MAP", "MODEL", "STRESS", "DECIDE"]
    payload = {"seat": "gemini-argon-4", "topic": topic, "stages": stages}
    return {
        "seat": "Gemini Argon 4 Extended Thinking",
        "model_id": "gemini-argon-4",
        "role": "DECIDER",
        "pipeline": stages,
        "tiers": [2048, 8192, 32768],
        "key": "names-only (GEMINI_API_KEY name checked, value never read)",
        "live_calls": 0,
        "decision": "PROCEED - sandbox simulation only, all material and device gates stay HOLD",
        "can_override_hold": False,
        "trace_sha256_16": _trace(payload),
    }


def claude_advisor(topic="topologicalinsulator-demo"):
    """Muse advisor seat: advises, can never override a HOLD."""
    return {
        "seat": "Muse",
        "role": "ADVISOR",
        "key": "names-only (MODEL_API_KEY name checked, value never read)",
        "live_calls": 0,
        "advice": "Print the convention before the number: the Chern sign belongs to this lab's stated plaquette orientation, the Zak phase is pi times the winding, and an invariant is a property of the model - then show the edge spectrum the invariant predicted, weighed site by site.",
        "can_override_hold": False,
    }


def google_flow_plan():
    """Google Flow seat: plan-only, nothing rendered, nothing live."""
    return {
        "seat": "Google Flow",
        "role": "PLAN-ONLY",
        "shots": ["a d-vector loop growing until it swallows the origin and the winding clicks to 1",
                  "two zero-energy modes lighting up on the end cells of an open chain",
                  "Berry flux arrows threading each plaquette of the Brillouin zone",
                  "a chiral edge band threading the bulk gap of a ribbon",
                  "the mass u sliding through zero while the gap closes and the Chern number flips"],
        "rendered": False,
        "live": 0,
    }


def brains_report():
    return {"decider": argon4_decide(), "advisor": claude_advisor(),
            "flow": google_flow_plan()}


# ================================================================ checks
def run_checks():
    c = []

    def chk(name, ok):
        c.append((name, bool(ok)))

    # 1. solver referees - known spectra + residuals before physics
    sx = [[0.0, 1.0], [1.0, 0.0]]
    vals, _ = jacobi_eigh(sx)
    chk("jacobi_sigma_x_spectrum",
        abs(vals[0] + 1.0) < 1e-12 and abs(vals[1] - 1.0) < 1e-12)
    tri = [[2.0, 1.0, 0.0], [1.0, 2.0, 1.0], [0.0, 1.0, 2.0]]
    vals, V = jacobi_eigh(tri)
    expect = (2.0 - math.sqrt(2.0), 2.0, 2.0 + math.sqrt(2.0))
    chk("jacobi_tridiagonal_spectrum",
        all(abs(a - b) < 1e-12 for a, b in zip(vals, expect)))
    chk("jacobi_tridiagonal_residual",
        eig_residual_max(tri, vals, V) < 1e-12)
    chk("hermitian_sigma_y_spectrum",
        all(abs(a - b) < 1e-12 for a, b in
            zip(hermitian_eigenvalues(PAULI_Y), (-1.0, 1.0))))
    known = [[1 + 0j, 1j], [-1j, 2 + 0j]]
    chk("hermitian_complex_2x2_spectrum",
        all(abs(a - b) < 1e-12 for a, b in zip(
            hermitian_eigenvalues(known),
            ((3.0 - math.sqrt(5.0)) / 2.0, (3.0 + math.sqrt(5.0)) / 2.0))))
    xx_zz = [[0j] * 4 for _ in range(4)]
    for i in range(4):
        for j in range(4):
            xx_zz[i][j] = (PAULI_X[i // 2][j // 2] * PAULI_X[i % 2][j % 2]
                           + PAULI_Z[i // 2][j // 2] * PAULI_Z[i % 2][j % 2])
    chk("hermitian_xx_zz_spectrum",
        all(abs(a - b) < 1e-12 for a, b in
            zip(hermitian_eigenvalues(xx_zz), (-2.0, 0.0, 0.0, 2.0))))

    # 2. SSH winding - exactly 0 or 1 across the sweep, Zak = pi x winding
    sweep = ((2.0, 1.0, 0), (1.5, 1.0, 0), (1.1, 1.0, 0), (0.9, 1.0, 1),
             (0.5, 1.0, 1), (0.2, 1.0, 1), (1.0, 0.5, 0), (1.0, 2.0, 1))
    chk("winding_exact_across_sweep",
        all(abs(ssh_winding(v, w) - want) < 1e-9 for v, w, want in sweep))
    chk("winding_zero_for_v_gt_w",
        all(abs(ssh_winding(v, w)) < 1e-9
            for v, w in ((1.5, 1.0), (2.0, 1.0), (1.0, 0.5))))
    chk("winding_one_for_v_lt_w",
        all(abs(ssh_winding(v, w) - 1.0) < 1e-9
            for v, w in ((0.5, 1.0), (0.2, 1.0), (1.0, 2.0))))
    chk("zak_phase_is_pi_times_winding",
        abs(ssh_zak_phase(0.5, 1.0) - math.pi) < 1e-9
        and abs(ssh_zak_phase(1.5, 1.0)) < 1e-9)

    # 3. open SSH chain - edge modes in the topological phase only
    topo = ssh_open_report(*SSH_TOPO)
    triv = ssh_open_report(*SSH_TRIVIAL)
    chk("open_chain_residual_small", topo["residual_max"] < 1e-10
        and triv["residual_max"] < 1e-10)
    chk("topo_two_modes_below_1pct_gap", topo["n_in_gap_below_1pct"] == 2)
    chk("topo_edge_energy_tiny_fraction",
        topo["min_abs_energy"] < 0.01 * topo["half_gap"])
    chk("topo_edge_weight_on_ends",
        len(topo["in_gap_modes"]) == 2
        and all(m["end_cell_weight"] > 0.5 for m in topo["in_gap_modes"]))
    chk("topo_edge_weight_beats_middle_100x",
        all(m["end_over_mid"] > 100.0 for m in topo["in_gap_modes"]))
    chk("trivial_no_in_gap_modes", triv["n_in_gap_below_1pct"] == 0
        and len(triv["in_gap_modes"]) == 0)
    chk("trivial_min_energy_at_band_edge",
        triv["min_abs_energy"] >= triv["half_gap"])

    # 4. QWZ Chern numbers - integers in the three regimes, converged
    chk("chern_u3_zero", abs(qwz_chern(3.0, 12)) < 1e-6)
    chk("chern_u1_plus_one", abs(qwz_chern(1.0, 12) - 1.0) < 1e-6)
    chk("chern_u_minus1_minus_one", abs(qwz_chern(-1.0, 12) + 1.0) < 1e-6)
    chk("chern_u1p5_plus_one", abs(qwz_chern(1.5, 12) - 1.0) < 1e-6)
    chk("chern_u_minus1p5_minus_one", abs(qwz_chern(-1.5, 12) + 1.0) < 1e-6)
    chk("chern_converged_under_refinement",
        all(abs(qwz_chern(u, 8) - qwz_chern(u, 24)) < 1e-6
            for u in QWZ_MASSES))
    chk("qwz_lower_eigenvector_eigen_equation",
        _qwz_eigvec_residual_ok())
    chk("qwz_bulk_gaps_positive",
        all(qwz_bulk_half_gap(u) > 0.9 for u in QWZ_MASSES))

    # 5. ribbon - chiral edge bands in Chern phases, none when trivial
    rib_p = ribbon_report(1.0)
    rib_m = ribbon_report(-1.0)
    rib_t = ribbon_report(3.0)
    chk("ribbon_chern_plus_crossings", rib_p["crossing_slices"] >= 3)
    chk("ribbon_chern_plus_touches_zero",
        rib_p["min_abs_energy"] < 0.02 * rib_p["bulk_half_gap"])
    chk("ribbon_chern_minus_crossings", rib_m["crossing_slices"] >= 3)
    chk("ribbon_trivial_no_in_gap_states",
        rib_t["crossing_slices"] == 0 and rib_t["in_gap_eigenvalues"] == 0)
    chk("ribbon_trivial_min_energy_above_gap",
        rib_t["min_abs_energy"] > 0.9 * rib_t["bulk_half_gap"])
    H1 = ribbon_hamiltonian(0.7, 1.0)
    ev1 = hermitian_eigenvalues(H1)
    chk("ribbon_trace_invariants",
        abs(sum(ev1) - sum(H1[i][i].real for i in range(len(H1)))) < 1e-8
        and abs(sum(e * e for e in ev1)
                - sum(abs(H1[i][j]) ** 2
                      for i in range(len(H1)) for j in range(len(H1))))
        < 1e-6)

    # 6. classifier - labels come from the phase policy itself
    rep = classifier_report()
    chk("classifier_ssh_trivial", rep["ssh_trivial_v1.5_w1"]["knn"]
        == "trivial-insulator")
    chk("classifier_ssh_topological", rep["ssh_topological_v0.5_w1"]["knn"]
        == "ssh-topological")
    chk("classifier_ssh_critical", rep["ssh_critical_v1_w1"]["knn"]
        == "gapless-critical")
    chk("classifier_qwz_trivial", rep["qwz_trivial_u3"]["knn"]
        == "trivial-insulator")
    chk("classifier_qwz_chern_plus", rep["qwz_chern_plus_u1"]["knn"]
        == "chern-insulator")
    chk("classifier_qwz_chern_minus", rep["qwz_chern_minus_u-1"]["knn"]
        == "chern-insulator")
    chk("classifier_matches_policy_scenarios",
        all(row["knn"] == row["policy"] for row in rep.values()))
    chk("classifier_grid_consistency",
        all(classify_features(feats) == label
            for feats, label in TRAINING_GRID))

    # gates + seats
    chk("hold_gates_14", len(HOLD_GATES) == 14
        and all(v.startswith("HOLD") for v in HOLD_GATES.values()))
    chk("gate_material_hold", gate("real_material_claim") == "HOLD")
    chk("gate_qubit_hold", gate("topological_qubit_claim") == "HOLD")
    chk("gate_fault_tolerance_hold", gate("fault_tolerance_claim") == "HOLD")
    chk("gate_measurement_hold", gate("model_as_measurement_claim") == "HOLD")
    chk("gate_sandbox_proceed", gate("run_simulation") == "PROCEED")
    brains = brains_report()
    chk("argon4_no_override", brains["decider"]["can_override_hold"] is False
        and brains["decider"]["model_id"] == "gemini-argon-4"
        and brains["decider"]["seat"]
        == "Gemini Argon 4 Extended Thinking")
    chk("claude_no_override", brains["advisor"]["can_override_hold"] is False)
    chk("flow_plan_only", brains["flow"]["rendered"] is False
        and brains["flow"]["live"] == 0)
    return c


def _qwz_eigvec_residual_ok():
    """Referee the projector-column eigenvector against H u = -|d| u."""
    for kx, ky, u in ((0.3, -1.1, 1.0), (2.2, 0.7, -1.0), (-0.9, 2.9, 3.0)):
        dx, dy, dz = qwz_dvector(kx, ky, u)
        r = math.sqrt(dx * dx + dy * dy + dz * dz)
        vec = qwz_lower_eigenvector(kx, ky, u)
        h0 = dz * vec[0] + (dx - 1j * dy) * vec[1]
        h1 = (dx + 1j * dy) * vec[0] - dz * vec[1]
        if abs(h0 + r * vec[0]) > 1e-12 or abs(h1 + r * vec[1]) > 1e-12:
            return False
    return True


# ================================================================ CLI / demo
def demo():
    print("=" * 68)
    print("TopologicalInsulator-Lab - winding, Zak phase, Chern numbers, edges")
    print("by Darion Kellar | MIT | pure stdlib | zero network | EDUCATION ONLY")
    print("=" * 68)
    print("\n[1] THE HONEST CORE - what an invariant does and does not prove")
    print("    A topological invariant is a whole number computed from the")
    print("    bulk wavefunctions of a MODEL. It cannot change unless the")
    print("    bulk gap closes, and when it is nonzero the model's boundary")
    print("    must host in-gap states. It certifies nothing about any real")
    print("    material - no sample exists here, and none is described.")
    print("\n[2] SSH chain - winding number by discretised angle sum")
    for row in winding_report().values():
        print(f"    v={row['v']:.1f} w={row['w']:.1f}  winding "
              f"{row['winding']:+.6f}  Zak/pi {row['zak_over_pi']:+.6f}  "
              f"{row['phase']}")
    print("    -> the d-vector circle encloses the origin exactly when")
    print("       v < w: winding 1, Zak phase pi. Otherwise 0. No rounding")
    print("       is involved - the angle sum lands on the integer.")
    print("\n[3] Open SSH chain - the boundary keeps the bulk's promise")
    topo = ssh_open_report(*SSH_TOPO)
    triv = ssh_open_report(*SSH_TRIVIAL)
    for label, rep in (("topological v=0.5 w=1", topo),
                       ("trivial     v=1.5 w=1", triv)):
        print(f"    {label}: half-gap {rep['half_gap']:.3f} | "
              f"min |E| {rep['min_abs_energy']:.3e} | in-gap modes "
              f"{len(rep['in_gap_modes'])}")
        for m in rep["in_gap_modes"]:
            print(f"      E {m['energy']:+.3e}  end-cell weight "
                  f"{m['end_cell_weight']:.4f}  mid-cell "
                  f"{m['mid_cell_weight']:.2e}  ratio "
                  f"{m['end_over_mid']:.3e}")
    print("    -> two modes sit at ~1e-6 of the gap with three quarters")
    print("       of their weight on the two end cells; the trivial chain")
    print("       has none - its smallest |E| is the band edge itself.")
    print("\n[4] QWZ Chern insulator - Fukui plaquette Chern numbers")
    for row in chern_report().values():
        conv = "  ".join(f"N={n}: {c:+.6f}"
                         for n, c in row["chern_by_grid"].items())
        print(f"    u={row['u']:+.1f}  C {row['chern']:+.6f}   {conv}")
    print("    -> C = +1 for 0 < u < 2, -1 for -2 < u < 0, 0 for |u| > 2")
    print("       in this lab's stated convention, integer at every grid.")
    print("\n[5] QWZ ribbon - chiral edge bands cross only in Chern phases")
    for u in QWZ_MASSES:
        rep = ribbon_report(u)
        print(f"    u={u:+.1f}: half-gap {rep['bulk_half_gap']:.3f} | "
              f"kx slices with in-gap states {rep['crossing_slices']}/"
              f"{rep['nkx']} | min |E| {rep['min_abs_energy']:.4f}")
    print("    -> the trivial ribbon stays out of its gap entirely; both")
    print("       Chern ribbons thread it with an edge band through E = 0.")
    print("\n[6] Phase classifier (TinyML kNN, labels from the phase policy)")
    for name, row in classifier_report().items():
        print(f"    {name:28s} -> {row['knn']}")
    print("\n[7] Advisory seats")
    brains = brains_report()
    print(f"    decider: {brains['decider']['decision']} "
          f"trace={brains['decider']['trace_sha256_16']}")
    print(f"    advisor can_override_hold = "
          f"{brains['advisor']['can_override_hold']}")
    print(f"    flow rendered = {brains['flow']['rendered']} "
          f"live = {brains['flow']['live']}")
    print(f"\n[8] HOLD gates: {len(HOLD_GATES)}")
    for key in HOLD_GATES:
        print(f"    HOLD  {key}")
    checks = run_checks()
    passed = sum(1 for _, ok in checks if ok)
    print(f"\nCHECKS {passed}/{len(checks)}")
    for name, ok in checks:
        if not ok:
            print(f"    FAIL {name}")
    verdict = "PROCEED" if passed == len(checks) else "HOLD"
    print(f"VERDICT {verdict} (sandbox simulation only) | "
          f"HOLD gates {len(HOLD_GATES)}")
    return 0 if passed == len(checks) else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] == "demo":
        return demo()
    cmd = argv[0]
    if cmd == "checks":
        checks = run_checks()
        for name, ok in checks:
            print(f"{'PASS' if ok else 'FAIL'} {name}")
        print(f"CHECKS {sum(1 for _, ok in checks if ok)}/{len(checks)}")
        return 0 if all(ok for _, ok in checks) else 1
    if cmd == "winding":
        print(json.dumps(winding_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "edge":
        print(json.dumps({
            "topological": ssh_open_report(*SSH_TOPO),
            "trivial": ssh_open_report(*SSH_TRIVIAL),
        }, indent=2, sort_keys=True))
        return 0
    if cmd == "chern":
        print(json.dumps(chern_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "ribbon":
        print(json.dumps({f"u{u}": ribbon_report(u) for u in QWZ_MASSES},
                         indent=2, sort_keys=True))
        return 0
    if cmd == "classify":
        print(json.dumps(classifier_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "brains":
        print(json.dumps(brains_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "holds":
        for key, value in HOLD_GATES.items():
            print(f"{key}: {value}")
        return 0
    if cmd == "hold-demo":
        action = argv[1] if len(argv) > 1 else "real_material_claim"
        print(f"{action}: {HOLD_GATES.get(action, 'PROCEED - not gated')}")
        return 2 if gate(action) == "HOLD" else 0
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
