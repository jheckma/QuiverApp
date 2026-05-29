"""Command-line / interactive entry point.

    python -m conformalmanifold                 # interactive menu
    python -m conformalmanifold "A4 = Delta(12)"
    python -m conformalmanifold --cyclic 10 2 3 5
    python -m conformalmanifold --list
"""

from __future__ import annotations

import sys

from .groups import cyclic, list_groups, make_group
from .pipeline import run


def _interactive():
    names = list_groups()
    print("Choose a finite subgroup Gamma < SU(3) for C^3/Gamma:\n")
    for i, n in enumerate(names):
        print(f"  [{i}] {n}")
    print("  [c] custom cyclic Z_n with weights (a,b,c)")
    choice = input("\nselection> ").strip()
    if choice.lower() == "c":
        n = int(input("  n = "))
        a = int(input("  a = "))
        b = int(input("  b = "))
        c = int(input("  c = "))
        return cyclic(n, (a, b, c))
    return make_group(names[int(choice)])


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    if argv and argv[0] == "--list":
        for n in list_groups():
            print(n)
        return 0

    if argv and argv[0] == "--serve":
        from .webapp import serve
        port = int(argv[argv.index("--port") + 1]) if "--port" in argv else 8000
        serve(port)
        return 0

    if argv and argv[0] == "--cyclic":
        n, a, b, c = (int(x) for x in argv[1:5])
        run(cyclic(n, (a, b, c)))
        return 0

    if argv:
        run(make_group(argv[0]))
        return 0

    try:
        group = _interactive()
    except (EOFError, KeyboardInterrupt):
        print("\naborted.")
        return 1
    run(group)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
