"""A tiny dependency-free web frontend for the C^3/Gamma pipeline.

    python -m conformalmanifold.webapp           # serve on http://localhost:8000
    python -m conformalmanifold.webapp --port 9000

Endpoints
    GET /                       the single-page UI
    GET /api/groups             {"groups": [...names...]}
    GET /api/compute?name=...   full pipeline summary as JSON
    GET /api/compute?name=__cyclic__&n=10&a=2&b=3&c=5
    GET /api/toric_web?pts=x0,y0;x1,y1;...   toric diagram -> (p,q) web + quiver
"""

from __future__ import annotations

import json
import os
import sys
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs, urlparse

from . import api
from .groups import list_groups

_HERE = os.path.dirname(os.path.abspath(__file__))
_INDEX = os.path.join(_HERE, "static", "index.html")


def _parse_points(s: str):
    """'x0,y0;x1,y1;...' -> [(x0,y0), (x1,y1), ...]  (integer lattice points)."""
    pts = []
    for chunk in s.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        xy = chunk.split(",")
        if len(xy) != 2:
            raise ValueError(f"bad point {chunk!r}; expected 'x,y'")
        pts.append((int(xy[0]), int(xy[1])))
    if not pts:
        raise ValueError("no toric-diagram points supplied")
    return pts


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):  # quiet
        pass

    def _send(self, code: int, body: bytes, ctype: str):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _json(self, obj, code: int = 200):
        self._send(code, json.dumps(obj).encode("utf-8"), "application/json")

    def do_GET(self):
        parsed = urlparse(self.path)
        path = parsed.path

        if path in ("/", "/index.html"):
            with open(_INDEX, "rb") as fh:
                self._send(200, fh.read(), "text/html; charset=utf-8")
            return

        if path == "/api/groups":
            self._json({"groups": list_groups()})
            return

        if path == "/api/compute":
            q = parse_qs(parsed.query)
            name = q.get("name", [""])[0]
            try:
                if name == "__cyclic__":
                    n = int(q.get("n", ["0"])[0])
                    a = int(q.get("a", ["0"])[0])
                    b = int(q.get("b", ["0"])[0])
                    c = int(q.get("c", ["0"])[0])
                    self._json(api.summarize_cyclic(n, a, b, c))
                elif name == "__freeform__":
                    self._json(api.summarize_freeform(q.get("expr", [""])[0]))
                else:
                    self._json(api.summarize_named(name))
            except (ValueError, KeyError) as exc:
                self._json({"error": str(exc)}, code=400)
            except Exception as exc:  # noqa: BLE001
                self._json({"error": f"internal error: {exc}"}, code=500)
            return

        if path == "/api/toric_web":
            q = parse_qs(parsed.query)
            try:
                pts = _parse_points(q.get("pts", [""])[0])
                self._json(api.summarize_toric_web(pts))
            except (ValueError, KeyError) as exc:
                self._json({"error": str(exc)}, code=400)
            except Exception as exc:  # noqa: BLE001
                self._json({"error": f"internal error: {exc}"}, code=500)
            return

        self._json({"error": "not found"}, code=404)


def serve(port: int = 8000):
    httpd = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    url = f"http://localhost:{port}"
    print(f"conformalmanifold UI serving at {url}  (Ctrl-C to stop)")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped.")


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    port = 8000
    if "--port" in argv:
        port = int(argv[argv.index("--port") + 1])
    serve(port)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
