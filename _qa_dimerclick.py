"""Throwaway QA: clickable square faces in the dimer (Seiberg duality by
clicking the brane tiling).  Headless Chrome over CDP, same pattern as
_qa_driver.py.  Server must already be up on :8011."""
import base64, json, os, subprocess, time, urllib.request
import websocket

CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
PORT = 9351
QADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_qa")
USERDIR = os.path.join(QADIR, f"chrome-profile-{os.getpid()}")
BASE = "http://localhost:8000"
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

def navigate(url, settle=5.0):
    send("Page.navigate", url="about:blank")
    drain(0.6)
    send("Page.navigate", url=url)
    drain(settle)

def shot(name):
    r = send("Page.captureScreenshot", format="png", captureBeyondViewport=True)
    with open(os.path.join(QADIR, name), "wb") as f:
        f.write(base64.b64decode(r["result"]["data"]))
    return name

report = {}

# ---- F0: all 4 faces square -> 4 glowing clickable faces in the tiling ----
console.clear()
navigate(f"{BASE}/#pts=-1,0;1,0;0,1;0,-1", settle=6)
report["f0_before"] = {
    "shot": shot("dc_01_f0.png"),
    "iv_badge": evaluate("(document.getElementById('iv-badge')||{}).textContent"),
    "n_dface": evaluate("document.querySelectorAll('#iv-tiling .dface').length"),
    "n_dface_sq": evaluate("document.querySelectorAll('#iv-tiling .dface.sq').length"),
    "n_clickable": evaluate("document.querySelectorAll('#iv-tiling polygon[data-qnode]').length"),
    "labels": evaluate("[...new Set([...document.querySelectorAll('#iv-tiling .dface text')].map(t=>t.textContent))].sort()"),
    "fields": evaluate("(document.getElementById('iv-fields')||{}).textContent"),
    "console": list(console),
}

# ---- click a square face (node 0) in the TILING -> urban renewal ----
console.clear()
evaluate("document.querySelector('#iv-tiling polygon[data-qnode=\"0\"]').dispatchEvent(new MouseEvent('click',{bubbles:true})); null")
drain(6)
report["f0_after_face_click"] = {
    "shot": shot("dc_02_f0_dualized.png"),
    "hash": evaluate("location.hash"),
    "dual_path_badge": evaluate("(document.getElementById('iv-dual-path')||{}).textContent"),
    "iv_badge": evaluate("(document.getElementById('iv-badge')||{}).textContent"),
    "fields": evaluate("(document.getElementById('iv-fields')||{}).textContent"),
    "n_dface_sq": evaluate("document.querySelectorAll('#iv-tiling .dface.sq').length"),
    "labels_sq": evaluate("[...new Set([...document.querySelectorAll('#iv-tiling .dface.sq text')].map(t=>t.textContent))].sort()"),
    "console": list(console),
}

# ---- click AGAIN on the same face (involution: back to 8 fields) ----
console.clear()
evaluate("document.querySelector('#iv-tiling polygon[data-qnode=\"0\"]').dispatchEvent(new MouseEvent('click',{bubbles:true})); null")
drain(6)
report["f0_after_second_click"] = {
    "shot": shot("dc_03_f0_backagain.png"),
    "hash": evaluate("location.hash"),
    "dual_path_badge": evaluate("(document.getElementById('iv-dual-path')||{}).textContent"),
    "fields": evaluate("(document.getElementById('iv-fields')||{}).textContent"),
    "console": list(console),
}

# ---- C3: hexagonal face, NO square faces -> labels but nothing clickable ----
console.clear()
navigate(f"{BASE}/#pts=0,0;1,0;0,1", settle=6)
report["c3"] = {
    "shot": shot("dc_04_c3.png"),
    "n_dface": evaluate("document.querySelectorAll('#iv-tiling .dface').length"),
    "n_dface_sq": evaluate("document.querySelectorAll('#iv-tiling .dface.sq').length"),
    "n_clickable": evaluate("document.querySelectorAll('#iv-tiling polygon[data-qnode]').length"),
    "console": list(console),
}

# ---- deep-link: &dualize= restores the face-clicked state ----
console.clear()
navigate(f"{BASE}/#pts=-1,0;1,0;0,1;0,-1&dualize=0", settle=6)
report["deeplink"] = {
    "shot": shot("dc_05_deeplink.png"),
    "fields": evaluate("(document.getElementById('iv-fields')||{}).textContent"),
    "dual_path_badge": evaluate("(document.getElementById('iv-dual-path')||{}).textContent"),
    "console": list(console),
}

print(json.dumps(report, indent=2, default=str))
ws.close()
proc.terminate()
