"""Finite-chamber BPS spectra for quivers.

The chamber is encoded by a mutation sequence.  For a green mutation sequence,
the positive c-vector at each mutation is the charge/dimension vector of a
stable hypermultiplet in that finite BPS chamber.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
import re


MAX_NODES = 20


@dataclass
class BpsInput:
    matrix: list[list[int]]
    labels_are_zero_based: bool


def _as_int_matrix(data) -> list[list[int]]:
    if not isinstance(data, list) or not data:
        raise ValueError("exchange matrix must be a non-empty square array")
    matrix: list[list[int]] = []
    for row in data:
        if not isinstance(row, list):
            raise ValueError("exchange matrix rows must be arrays")
        try:
            matrix.append([int(v) for v in row])
        except (TypeError, ValueError) as exc:
            raise ValueError("exchange matrix entries must be integers") from exc
    return matrix


def parse_matrix(text: str) -> list[list[int]]:
    """Parse a skew-symmetric integer exchange matrix."""
    text = (text or "").strip()
    if not text:
        raise ValueError("empty exchange matrix")
    if text.startswith("["):
        matrix = _as_int_matrix(json.loads(text))
    else:
        rows = [r.strip() for r in re.split(r"[\n;]+", text) if r.strip()]
        matrix = []
        for row in rows:
            vals = [v for v in re.split(r"[\s,]+", row) if v]
            matrix.append([int(v) for v in vals])
    _validate_exchange_matrix(matrix)
    return matrix


def parse_arrows(text: str) -> BpsInput:
    """Parse arrows like ``1->2``, ``1->2 x2`` or ``0->1:3``."""
    text = (text or "").strip()
    if not text:
        raise ValueError("empty arrow list")

    chunks = [c.strip() for c in re.split(r"[\n,;]+", text) if c.strip()]
    arrows: list[tuple[int, int, int]] = []
    directed: set[tuple[int, int]] = set()
    labels: list[int] = []
    pat = re.compile(r"^(\d+)\s*->\s*(\d+)(?:\s*(?:x|\*|:)\s*(\d+))?$")
    for chunk in chunks:
        m = pat.match(chunk)
        if not m:
            raise ValueError(
                f"bad arrow {chunk!r}; expected forms like '1->2' or '1->2 x2'"
            )
        a, b = int(m.group(1)), int(m.group(2))
        mult = int(m.group(3) or 1)
        if a == b:
            raise ValueError("BPS quivers should not contain self-arrows")
        if mult <= 0:
            raise ValueError("arrow multiplicities must be positive")
        if (b, a) in directed:
            raise ValueError("arrow list contains a 2-cycle; BPS quivers must be 2-acyclic")
        directed.add((a, b))
        arrows.append((a, b, mult))
        labels.extend([a, b])

    zero_based = 0 in labels
    offset = 0 if zero_based else 1
    n = max(labels) - offset + 1
    if n <= 0:
        raise ValueError("could not infer nodes from arrow list")
    matrix = [[0 for _ in range(n)] for _ in range(n)]
    for a, b, mult in arrows:
        i, j = a - offset, b - offset
        if not (0 <= i < n and 0 <= j < n):
            raise ValueError("mixed 0-based and 1-based node labels")
        matrix[i][j] += mult
        matrix[j][i] -= mult

    _validate_exchange_matrix(matrix)
    return BpsInput(matrix=matrix, labels_are_zero_based=zero_based)


def parse_quiver(text: str, kind: str = "matrix") -> BpsInput:
    kind = (kind or "matrix").lower()
    if kind == "auto":
        kind = "arrows" if "->" in (text or "") else "matrix"
    if kind == "arrows":
        return parse_arrows(text)
    if kind != "matrix":
        raise ValueError("quiver kind must be 'matrix', 'arrows', or 'auto'")
    return BpsInput(matrix=parse_matrix(text), labels_are_zero_based=False)


def _validate_exchange_matrix(matrix: list[list[int]]) -> None:
    n = len(matrix)
    if n > MAX_NODES:
        raise ValueError(f"BPS quiver UI is limited to {MAX_NODES} nodes")
    if any(len(row) != n for row in matrix):
        raise ValueError("exchange matrix must be square")
    for i in range(n):
        if matrix[i][i] != 0:
            raise ValueError("exchange matrix diagonal must vanish")
        for j in range(n):
            if matrix[i][j] != -matrix[j][i]:
                raise ValueError("exchange matrix must be skew-symmetric")


def _parse_number_list(text: str, n: int, name: str) -> list[float] | None:
    text = (text or "").strip()
    if not text:
        return None
    if text.startswith("["):
        data = json.loads(text)
        vals = [float(v) for v in data]
    else:
        vals = [float(v) for v in re.split(r"[\s,;]+", text) if v]
    if len(vals) != n:
        raise ValueError(f"{name} must contain exactly {n} numbers")
    if name == "masses" and any(v <= 0 for v in vals):
        raise ValueError("masses must be positive")
    return vals


def parse_sequence(text: str, n: int, zero_based_hint: bool = False) -> list[int]:
    text = (text or "").strip()
    if not text:
        raise ValueError("empty chamber mutation sequence")
    vals = [int(v) for v in re.findall(r"\d+", text)]
    if not vals:
        raise ValueError("mutation sequence must contain node labels")
    zero_based = zero_based_hint or 0 in vals
    offset = 0 if zero_based else 1
    seq = [v - offset for v in vals]
    bad = [vals[i] for i, k in enumerate(seq) if k < 0 or k >= n]
    if bad:
        lo, hi = (0, n - 1) if zero_based else (1, n)
        raise ValueError(f"mutation node {bad[0]} is outside the range {lo}..{hi}")
    return seq


def _mutate_exchange(B: list[list[int]], k: int) -> list[list[int]]:
    n = len(B)
    out = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == k or j == k:
                out[i][j] = -B[i][j]
            else:
                out[i][j] = (
                    B[i][j]
                    + max(B[i][k], 0) * max(B[k][j], 0)
                    - max(B[j][k], 0) * max(B[k][i], 0)
                )
    return out


def _vector_sign(v: list[int]) -> str:
    if all(x >= 0 for x in v) and any(x > 0 for x in v):
        return "green"
    if all(x <= 0 for x in v) and any(x < 0 for x in v):
        return "red"
    return "mixed"


def _mutate_c_matrix(C: list[list[int]], B: list[list[int]], k: int, sign: str):
    n = len(C)
    ck = C[k][:]
    sgn = 1 if sign == "green" else (-1 if sign == "red" else 1)
    out = [c[:] for c in C]
    for j in range(n):
        if j == k:
            out[j] = [-x for x in ck]
        else:
            factor = max(sgn * B[k][j], 0)
            if factor:
                out[j] = [out[j][a] + factor * ck[a] for a in range(n)]
    return out


def _basis_label(v: list[int], label_base: int = 1) -> str:
    terms = []
    for i, coeff in enumerate(v, start=label_base):
        if coeff == 0:
            continue
        base = f"gamma_{i}"
        terms.append(base if coeff == 1 else f"{coeff} {base}")
    return " + ".join(terms) if terms else "0"


def _central_charge(v: list[int], phases: list[float], masses: list[float]) -> dict:
    z = 0j
    for coeff, phase, mass in zip(v, phases, masses):
        z += coeff * mass * complex(math.cos(math.pi * phase),
                                    math.sin(math.pi * phase))
    mag = abs(z)
    arg = math.atan2(z.imag, z.real) / math.pi
    if arg < 0:
        arg += 2.0
    return {
        "re": round(z.real, 8),
        "im": round(z.imag, 8),
        "mass": round(mag, 8),
        "phase": round(arg, 8),
    }


def finite_chamber_spectrum(matrix: list[list[int]], sequence: list[int],
                            phases: list[float] | None = None,
                            masses: list[float] | None = None,
                            label_base: int = 1) -> dict:
    n = len(matrix)
    B = [row[:] for row in matrix]
    C = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    stable = []
    mutations = []

    for step, k in enumerate(sequence, start=1):
        cvec = C[k][:]
        sign = _vector_sign(cvec)
        central = None
        if phases is not None and masses is not None:
            central = _central_charge(cvec, phases, masses)
        mutation = {
            "step": step,
            "node": k + label_base,
            "status": sign,
            "c_vector": cvec,
            "charge": _basis_label(cvec, label_base),
        }
        if central is not None:
            mutation["central_charge"] = central
        mutations.append(mutation)
        if sign == "green":
            stable.append({
                "index": len(stable) + 1,
                "step": step,
                "node": k + label_base,
                "dimension_vector": cvec,
                "charge": _basis_label(cvec, label_base),
                "kind": "hypermultiplet",
                **({"central_charge": central} if central is not None else {}),
            })
        C = _mutate_c_matrix(C, B, k, sign)
        B = _mutate_exchange(B, k)

    green = all(m["status"] == "green" for m in mutations)
    final_signs = [_vector_sign(c) for c in C]
    maximal_green = green and all(s == "red" for s in final_signs)

    phase_order_ok = None
    if stable and phases is not None and masses is not None:
        seq_phases = [s["central_charge"]["phase"] for s in stable]
        phase_order_ok = all(
            seq_phases[i] + 1e-8 >= seq_phases[i + 1]
            for i in range(len(seq_phases) - 1)
        )

    return {
        "stable": stable,
        "mutations": mutations,
        "final_c_matrix": C,
        "final_exchange_matrix": B,
        "green_sequence": green,
        "maximal_green": maximal_green,
        "final_c_vector_signs": final_signs,
        "phase_order_ok": phase_order_ok,
    }




def _validate_adjacency(adjacency) -> list[list[int]]:
    matrix = _as_int_matrix(adjacency)
    n = len(matrix)
    if n > MAX_NODES:
        raise ValueError(f"BPS quiver UI is limited to {MAX_NODES} nodes")
    if any(len(row) != n for row in matrix):
        raise ValueError("quiver adjacency must be square")
    for i, row in enumerate(matrix):
        for j, val in enumerate(row):
            if val < 0:
                raise ValueError("quiver adjacency entries must be non-negative")
    return matrix


def exchange_from_adjacency(adjacency) -> list[list[int]]:
    """Return the skew exchange matrix B_ij = A_ij - A_ji.

    Geometry builders return an oriented adjacency matrix.  BPS mutation uses the
    antisymmetric exchange matrix, so reciprocal arrows are netted and loops are
    ignored on the diagonal.
    """
    A = _validate_adjacency(adjacency)
    n = len(A)
    B = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                B[i][j] = A[i][j] - A[j][i]
    _validate_exchange_matrix(B)
    return B


def bps_quiver_from_adjacency_json(adjacency, sequence_text: str,
                                   phases_text: str = "", masses_text: str = "",
                                   label_base: int = 0,
                                   source: dict | None = None) -> dict:
    """Run BPS chamber analysis on an existing geometry-built quiver."""
    A = _validate_adjacency(adjacency)
    B = exchange_from_adjacency(A)
    n = len(B)
    sequence = parse_sequence(sequence_text, n, zero_based_hint=(label_base == 0))
    phases = _parse_number_list(phases_text, n, "phases")
    masses = _parse_number_list(masses_text, n, "masses")
    if masses is None and phases is not None:
        masses = [1.0] * n
    if phases is None and masses is not None:
        raise ValueError("masses require phases; phases may be given without masses")

    spectrum = finite_chamber_spectrum(B, sequence, phases, masses, label_base)
    warnings = []
    if any(A[i][i] for i in range(n)):
        warnings.append("Geometry quiver has self-arrows; loops are ignored in the BPS exchange matrix.")
    if any(A[i][j] and A[j][i] for i in range(n) for j in range(i + 1, n)):
        warnings.append("Geometry quiver has reciprocal arrows; BPS exchange matrix uses the net orientation A_ij - A_ji.")
    if all(B[i][j] == 0 for i in range(n) for j in range(n)):
        warnings.append("The antisymmetrized BPS exchange matrix is zero for this geometry quiver.")
    if not spectrum["green_sequence"]:
        warnings.append(
            "The mutation sequence is not green throughout; stable objects are "
            "reported only for green mutations before/inside that path."
        )
    if spectrum["green_sequence"] and not spectrum["maximal_green"]:
        warnings.append(
            "The sequence is green but not maximal green; this is a partial "
            "finite chamber, not a complete finite BPS spectrum."
        )
    if spectrum["phase_order_ok"] is False:
        warnings.append(
            "The supplied central-charge phases are not decreasing along the "
            "mutation order; check whether they describe the intended chamber."
        )

    return {
        "source": source or {"kind": "geometry"},
        "quiver": {
            "num_nodes": n,
            "exchange_matrix": B,
            "adjacency": A,
            "num_arrows": sum(sum(row) for row in A),
            "dims": [1] * n,
            "node_labels": [f"gamma_{i}" for i in range(label_base, label_base + n)],
            "label_base": label_base,
        },
        "chamber": {
            "mutation_sequence": [k + label_base for k in sequence],
            "num_mutations": len(sequence),
            "green_sequence": spectrum["green_sequence"],
            "maximal_green": spectrum["maximal_green"],
            "stable_count": len(spectrum["stable"]),
            "phase_order_ok": spectrum["phase_order_ok"],
            "warnings": warnings,
        },
        "stable_objects": spectrum["stable"],
        "mutations": spectrum["mutations"],
        "final": {
            "c_matrix": spectrum["final_c_matrix"],
            "exchange_matrix": spectrum["final_exchange_matrix"],
            "c_vector_signs": spectrum["final_c_vector_signs"],
        },
    }

def bps_quiver_json(quiver_text: str, sequence_text: str, kind: str = "matrix",
                    phases_text: str = "", masses_text: str = "") -> dict:
    parsed = parse_quiver(quiver_text, kind)
    matrix = parsed.matrix
    n = len(matrix)
    raw_sequence_labels = [int(v) for v in re.findall(r"\d+", sequence_text or "")]
    labels_are_zero_based = parsed.labels_are_zero_based or 0 in raw_sequence_labels
    sequence = parse_sequence(sequence_text, n, labels_are_zero_based)
    label_base = 0 if labels_are_zero_based else 1
    phases = _parse_number_list(phases_text, n, "phases")
    masses = _parse_number_list(masses_text, n, "masses")
    if masses is None and phases is not None:
        masses = [1.0] * n
    if phases is None and masses is not None:
        raise ValueError("masses require phases; phases may be given without masses")

    spectrum = finite_chamber_spectrum(matrix, sequence, phases, masses, label_base)
    adjacency = [[max(matrix[i][j], 0) for j in range(n)] for i in range(n)]
    warnings = []
    if not spectrum["green_sequence"]:
        warnings.append(
            "The mutation sequence is not green throughout; stable objects are "
            "reported only for green mutations before/inside that path."
        )
    if spectrum["green_sequence"] and not spectrum["maximal_green"]:
        warnings.append(
            "The sequence is green but not maximal green; this is a partial "
            "finite chamber, not a complete finite BPS spectrum."
        )
    if spectrum["phase_order_ok"] is False:
        warnings.append(
            "The supplied central-charge phases are not decreasing along the "
            "mutation order; check whether they describe the intended chamber."
        )

    return {
        "quiver": {
            "num_nodes": n,
            "exchange_matrix": matrix,
            "adjacency": adjacency,
            "num_arrows": sum(sum(row) for row in adjacency),
            "dims": [1] * n,
            "node_labels": [f"gamma_{i}" for i in range(label_base, label_base + n)],
            "label_base": label_base,
        },
        "chamber": {
            "mutation_sequence": [k + label_base for k in sequence],
            "num_mutations": len(sequence),
            "green_sequence": spectrum["green_sequence"],
            "maximal_green": spectrum["maximal_green"],
            "stable_count": len(spectrum["stable"]),
            "phase_order_ok": spectrum["phase_order_ok"],
            "warnings": warnings,
        },
        "stable_objects": spectrum["stable"],
        "mutations": spectrum["mutations"],
        "final": {
            "c_matrix": spectrum["final_c_matrix"],
            "exchange_matrix": spectrum["final_exchange_matrix"],
            "c_vector_signs": spectrum["final_c_vector_signs"],
        },
    }
