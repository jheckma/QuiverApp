"""Freeform group input -> MatrixGroup.

The user types either
  * a constructor expression, e.g.   cyclic(7, 1, 2, 4)   /   delta(6, 4)
  * explicit SU(3) generators, e.g.
        w = exp(2j*pi/7)
        [ diag(w, w**2, w**4),
          [[0,1,0],[0,0,1],[1,0,0]] ]
  * or a mix (assignments followed by a final expression).

The final expression (or a variable named group/G/Gamma/generators) is taken as
the result: a MatrixGroup is used directly; a single 3x3 matrix or a list of
3x3 matrices is fed to closure().  Evaluation happens in a restricted namespace
(numpy-backed helpers, no builtins beyond a safe handful) -- intended for local,
single-user use.
"""

from __future__ import annotations

import ast
import math

import numpy as np

from .groups import (MatrixGroup, closure, cyclic, delta_3n2, delta_6n2,
                     icosahedral_A5, octahedral_S4, tetrahedral_A4)


def _cyclic(n, a, b, c):
    return cyclic(int(n), (int(a), int(b), int(c)))


def _delta(k, n):
    if int(k) == 3:
        return delta_3n2(int(n))
    if int(k) == 6:
        return delta_6n2(int(n))
    raise ValueError("delta(k, n): the first argument k must be 3 or 6 "
                     "(Delta(3 n^2) or Delta(6 n^2)).")


def _diag(*xs):
    return np.diag(np.array(xs, dtype=complex))


_SAFE_BUILTINS = {
    "range": range, "len": len, "abs": abs, "complex": complex,
    "float": float, "int": int, "list": list, "tuple": tuple, "round": round,
}


def _namespace() -> dict:
    return {
        "__builtins__": _SAFE_BUILTINS,
        # group constructors
        "cyclic": _cyclic,
        "delta": _delta,
        "A4": tetrahedral_A4,
        "S4": octahedral_S4,
        "A5": icosahedral_A5,
        "closure": lambda gens, name="Gamma (freeform)": closure(
            [_to_matrix(g) for g in gens], name=name,
            description="user-supplied generators"),
        # math helpers (numpy-backed so they act on scalars and arrays)
        "exp": np.exp, "sqrt": np.sqrt, "cos": np.cos, "sin": np.sin,
        "pi": np.pi, "I": 1j, "phi": (1 + math.sqrt(5)) / 2,
        "diag": _diag, "eye": lambda: np.eye(3, dtype=complex),
        "array": lambda x: np.array(x, dtype=complex),
        "np": np,
    }


def _to_matrix(m) -> np.ndarray:
    arr = np.array(m, dtype=complex)
    if arr.shape != (3, 3):
        raise ValueError(f"each generator must be a 3x3 matrix; got shape {arr.shape}.")
    return arr


def build_from_expr(src: str, max_order: int = 5000) -> MatrixGroup:
    src = (src or "").strip()
    if not src:
        raise ValueError("empty input.")

    try:
        tree = ast.parse(src, mode="exec")
    except SyntaxError as exc:
        raise ValueError(f"syntax error: {exc.msg} (line {exc.lineno}).") from exc

    # turn a trailing expression into  __result__ = <expr>
    if tree.body and isinstance(tree.body[-1], ast.Expr):
        last = tree.body.pop()
        assign = ast.Assign(
            targets=[ast.Name(id="__result__", ctx=ast.Store())],
            value=last.value)
        ast.copy_location(assign, last)
        tree.body.append(assign)
    ast.fix_missing_locations(tree)

    ns = _namespace()
    try:
        exec(compile(tree, "<freeform>", "exec"), ns)  # noqa: S102 (local tool)
    except ValueError:
        raise
    except Exception as exc:  # noqa: BLE001
        raise ValueError(f"could not evaluate input: {exc}") from exc

    result = ns.get("__result__")
    if result is None:
        for key in ("group", "G", "Gamma", "generators", "gens"):
            if ns.get(key) is not None:
                result = ns[key]
                break
    if result is None:
        raise ValueError(
            "input produced no group. End with an expression that is a "
            "MatrixGroup (e.g. cyclic(7,1,2,4)) or a list of 3x3 generators, "
            "or assign it to `generators`.")

    if isinstance(result, MatrixGroup):
        return result

    arr = np.array(result, dtype=complex)
    if arr.shape == (3, 3):
        gens = [arr]
    elif arr.ndim == 3 and arr.shape[1:] == (3, 3):
        gens = [arr[k] for k in range(arr.shape[0])]
    else:
        raise ValueError(
            f"expected a 3x3 matrix or a list of 3x3 matrices; got shape "
            f"{arr.shape}.")

    return closure(gens, name="Gamma (freeform)",
                   description="user-supplied SU(3) generators",
                   max_order=max_order)
