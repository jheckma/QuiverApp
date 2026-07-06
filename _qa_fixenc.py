"""Diagnose + repair PS5.1 double-encoded UTF-8 (mojibake) in the repo's
uncommitted files, per the documented recipe: for each non-ASCII run, try
encode cp1252 (passing through the 5 undefined bytes 81 8D 8F 90 9D) then
decode utf-8; repair the run only if the round trip succeeds.  Runs that fail
are already-clean single-encoded UTF-8 and are left alone."""
import io, re, sys, unicodedata
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

PASSTHRU = {0x81, 0x8D, 0x8F, 0x90, 0x9D}

def de_mojibake_run(s):
    b = bytearray()
    for ch in s:
        o = ord(ch)
        if o in PASSTHRU:
            b.append(o)
            continue
        try:
            b += ch.encode("cp1252")
        except UnicodeEncodeError:
            return None
    try:
        out = b.decode("utf-8")
    except UnicodeDecodeError:
        return None
    return out

def scan(path, apply=False):
    raw = open(path, "rb").read()
    txt = raw.decode("utf-8")
    total = 0
    # iterate to fixpoint: a file rewritten k times needs k peeling passes
    # (index.html was hit twice: '✓' -> 'âœ“' -> 'Ã¢Å“â€œ')
    for _ in range(6):
        runs = sorted(set(re.findall(r"[^\x00-\x7f]+", txt)))
        fixed, clean = [], []
        for r in runs:
            rep = de_mojibake_run(r)
            # a genuine mojibake run repairs to something SHORTER; a clean run
            # (e.g. a lone arrow) can't cp1252-encode + utf-8-decode.
            if rep is not None and len(rep) < len(r):
                fixed.append((r, rep))
            else:
                clean.append(r)
        if not fixed:
            break
        total += len(fixed)
        for r, rep in sorted(fixed, key=lambda t: -len(t[0])):
            txt = txt.replace(r, rep)
    print(f"{path}: {total} mojibake runs peeled; "
          f"{len(clean)} clean runs remain: {[c[:6] for c in clean[:12]]!r}")
    if apply and total:
        crlf = b"\r\n" in raw
        if crlf:
            txt = txt.replace("\r\n", "\n")   # the rewrite also flipped LF->CRLF
        open(path, "wb").write(txt.encode("utf-8"))
        print(f"   -> written back (utf-8, no BOM{', CRLF->LF' if crlf else ''})")

if __name__ == "__main__":
    apply = "--apply" in sys.argv
    for p in ["conformalmanifold/static/index.html",
              "conformalmanifold/api.py",
              "conformalmanifold/webapp.py",
              "conformalmanifold/inverse.py",
              "conformalmanifold/bps.py",
              "tests/test_bps.py",
              "tests/test_inverse.py"]:
        try:
            scan(p, apply=apply)
        except FileNotFoundError:
            print(p, ": missing")
