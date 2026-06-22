"""Throwaway visual-QA driver: headless Chrome over CDP (websocket-client).

Drives the QuiverApp, captures full-page screenshots + every console
message / JS exception for each render path. Output -> _qa/ (untracked).
"""
import base64, json, os, subprocess, time, urllib.request
import websocket  # websocket-client 1.8.0

CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
PORT = 9347
QADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_qa")
USERDIR = os.path.join(QADIR, f"chrome-profile-{os.getpid()}")  # unique per run
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
# use the initial about:blank page target
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

console = []  # collected across the run
def drain(seconds):
    """Pump events for a while, collecting console + exceptions."""
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

send("Page.enable")
send("Runtime.enable")
send("Log.enable")

def evaluate(expr):
    r = send("Runtime.evaluate", expression=expr, awaitPromise=True, returnByValue=True)
    return r.get("result", {}).get("result", {}).get("value")

def navigate(url, settle=4.0):
    # reset to a blank document first so hash-based deep-link init re-runs
    send("Page.navigate", url="about:blank")
    drain(0.6)
    send("Page.navigate", url=url)
    drain(settle)

def shot(name):
    r = send("Page.captureScreenshot", format="png", captureBeyondViewport=True)
    data = r["result"]["data"]
    path = os.path.join(QADIR, name)
    with open(path, "wb") as f:
        f.write(base64.b64decode(data))
    return path

report = {}

# ---- 1. Toric tab: conifold (square) — both new cards + dimer ----
console.clear()
navigate(f"{BASE}/#pts=0,0;1,0;1,1;0,1", settle=5)
report["toric_conifold"] = {
    "shot": shot("01_toric_conifold.png"),
    "has_scft_card": evaluate("!!document.getElementById('ts-a') && document.getElementById('ts-a').textContent.length>0"),
    "ts_a": evaluate("(document.getElementById('ts-a')||{}).textContent"),
    "has_inverse_card": evaluate("!!document.getElementById('iv-tiling')"),
    "iv_badge": evaluate("(document.getElementById('iv-badge')||{}).textContent"),
    "tiling_svg_children": evaluate("var e=document.getElementById('iv-tiling'); e?e.childElementCount:-1"),
    "fieldR_shown": evaluate("!document.getElementById('ts-fieldR-wrap').classList.contains('hide')"),
    "fieldR_badge": evaluate("(document.getElementById('ts-fieldR-badge')||{}).textContent"),
    "fieldR_distinct": evaluate("(document.getElementById('ts-fieldR-distinct')||{}).textContent"),
    "fieldR_ntags": evaluate("var e=document.getElementById('ts-fieldR-tags'); e?e.childElementCount:-1"),
    "console": list(console),
}

# ---- 2. Toric tab: dP3 hexagon — richer dimer ----
console.clear()
navigate(f"{BASE}/#pts=1,0;1,-1;0,-1;-1,0;-1,1;0,1", settle=5)
report["toric_dP3"] = {
    "shot": shot("02_toric_dP3_hex.png"),
    "ts_a": evaluate("(document.getElementById('ts-a')||{}).textContent"),
    "iv_badge": evaluate("(document.getElementById('iv-badge')||{}).textContent"),
    "tiling_svg_children": evaluate("var e=document.getElementById('iv-tiling'); e?e.childElementCount:-1"),
    "fieldR_shown": evaluate("!document.getElementById('ts-fieldR-wrap').classList.contains('hide')"),
    "fieldR_badge": evaluate("(document.getElementById('ts-fieldR-badge')||{}).textContent"),
    "fieldR_distinct": evaluate("(document.getElementById('ts-fieldR-distinct')||{}).textContent"),
    "fieldR_ntags": evaluate("var e=document.getElementById('ts-fieldR-tags'); e?e.childElementCount:-1"),
    "console": list(console),
}

# ---- 3. Orbifold tab: select a group + compute() ----
console.clear()
navigate(f"{BASE}/", settle=2)
# pick Delta(27) (non-abelian, interesting) then trigger compute
evaluate("document.getElementById('group').value='Delta(27)'; null")
evaluate("compute(); null")
drain(4)
report["orbifold_delta27"] = {
    "shot": shot("03_orbifold_delta27.png"),
    "scft_card_a": evaluate("(document.getElementById('s-a')||{}).textContent"),
    "scft_trR": evaluate("(document.getElementById('s-trR')||{}).textContent"),
    "scft_badge": evaluate("(document.getElementById('s-badge')||{}).textContent"),
    "console": list(console),
}

# ---- 4. Orbifold tab: an abelian group Z3(1,1,1) = dP0 cross-check ----
console.clear()
evaluate("document.getElementById('group').value='Z3(1,1,1)'; null")
evaluate("compute(); null")
drain(4)
report["orbifold_z3"] = {
    "shot": shot("04_orbifold_z3.png"),
    "scft_card_a": evaluate("(document.getElementById('s-a')||{}).textContent"),
    "scft_trR": evaluate("(document.getElementById('s-trR')||{}).textContent"),
    "console": list(console),
}

print(json.dumps(report, indent=2, default=str))

ws.close()
proc.terminate()
