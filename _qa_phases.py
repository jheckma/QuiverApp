"""Throwaway visual-QA: dP2/dP3 Seiberg-phase cycling (mass-integration feature).

Drives the toric tab for dP2 and dP3, waits for the async inverse fetch, and
clicks through every toric phase, checking the phase label / field counts match
the literature (dP2: 2 phases 11+13 fields; dP3: 4 phases 12+14+14+18).
Output -> _qa/ (untracked).  Requires the webapp on :8011 and Chrome.
"""
import base64, json, os, subprocess, time, urllib.request
import websocket

CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
PORT = 9349
QADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_qa")
USERDIR = os.path.join(QADIR, f"chrome-profile-{os.getpid()}")
BASE = "http://localhost:8011"
os.makedirs(QADIR, exist_ok=True)

proc = subprocess.Popen(
    [CHROME, "--headless=new", "--disable-gpu", "--no-sandbox",
     "--hide-scrollbars", "--remote-allow-origins=*",
     f"--remote-debugging-port={PORT}",
     f"--user-data-dir={USERDIR}", "--window-size=1500,2600",
     "about:blank"],
    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def wait_devtools():
    for _ in range(100):
        try:
            with urllib.request.urlopen(f"http://localhost:{PORT}/json/version", timeout=1) as r:
                json.load(r); return
        except Exception:
            time.sleep(0.2)
    raise RuntimeError("devtools never came up")

wait_devtools()
with urllib.request.urlopen(f"http://localhost:{PORT}/json/list") as r:
    targets = json.load(r)
page = next(t for t in targets if t.get("type") == "page")
ws = websocket.create_connection(page["webSocketDebuggerUrl"], max_size=64*1024*1024)

_id = 0
def send(method, **params):
    global _id
    _id += 1
    ws.send(json.dumps({"id": _id, "method": method, "params": params}))
    while True:
        msg = json.loads(ws.recv())
        if msg.get("id") == _id:
            return msg

console = []
def drain(seconds):
    end = time.time() + seconds
    ws.settimeout(0.4)
    while time.time() < end:
        try:
            msg = json.loads(ws.recv())
        except Exception:
            continue
        m = msg.get("method")
        if m == "Runtime.consoleAPICalled":
            p = msg["params"]
            txt = " ".join(str(a.get("value", a.get("description",""))) for a in p.get("args",[]))
            console.append((p.get("type","log"), txt))
        elif m == "Runtime.exceptionThrown":
            d = msg["params"]["exceptionDetails"]
            console.append(("EXCEPTION", d.get("text","") + " " +
                            (d.get("exception",{}).get("description","") or "")))
    ws.settimeout(None)

send("Page.enable"); send("Runtime.enable"); send("Log.enable")

def evaluate(expr):
    r = send("Runtime.evaluate", expression=expr, awaitPromise=True, returnByValue=True)
    return r.get("result", {}).get("result", {}).get("value")

def navigate(url, settle=4.0):
    send("Page.navigate", url="about:blank")
    drain(0.6)
    send("Page.navigate", url=url)
    drain(settle)

def shot(name):
    r = send("Page.captureScreenshot", format="png", captureBeyondViewport=True)
    with open(os.path.join(QADIR, name), "wb") as f:
        f.write(base64.b64decode(r["result"]["data"]))

def wait_phases(timeout=30):
    """Wait until the async inverse fetch fills _ivPhases."""
    end = time.time() + timeout
    while time.time() < end:
        n = evaluate("typeof _ivPhases !== 'undefined' ? _ivPhases.length : -1")
        if n and n > 0:
            return n
        drain(0.7)
    return evaluate("typeof _ivPhases !== 'undefined' ? _ivPhases.length : -1")

report = {}
for label, pts, want in [
    ("dP2", "1,0;0,1;-1,0;-1,-1;0,-1",      [11, 13]),
    ("dP3", "1,0;0,1;-1,1;-1,0;0,-1;1,-1",  [12, 14, 14, 18]),
]:
    console.clear()
    navigate(f"{BASE}/#pts={pts}", settle=4)
    n = wait_phases()
    fields, labels, svg_ok = [], [], []
    for i in range(max(n, 0)):
        evaluate(f"showPhase({i}); null")
        drain(0.6)
        fields.append(evaluate("_ivPhases[_ivPhaseIdx].num_fields"))
        labels.append(evaluate("(document.getElementById('iv-phase-label')||{}).textContent"))
        svg_ok.append(evaluate("var e=document.getElementById('iv-tiling'); e?e.childElementCount:-1"))
        shot(f"phase_{label}_{i}.png")
    report[label] = {
        "n_phases": n,
        "fields": fields,
        "fields_ok": sorted(f for f in fields if f is not None) == want,
        "nav_visible": evaluate("!document.getElementById('iv-phase-nav').classList.contains('hide')"),
        "labels": labels,
        "tiling_svg_children": svg_ok,
        "hash": evaluate("location.hash"),
        "console_errors": [c for c in console if c[0] in ("error", "EXCEPTION")],
    }

print(json.dumps(report, indent=2))
ws.close(); proc.kill()
