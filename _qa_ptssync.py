"""QA: pts= hash sync — replay the user's flow that produced the chimera URL.
Load #pts=F0 -> click dP3 preset -> click a glowing dimer face -> the hash must
carry dP3's pts (not F0's), and reloading that hash must reproduce the display."""
import base64, json, os, subprocess, time, urllib.request
import websocket

CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
PORT = 9353
QADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_qa")
USERDIR = os.path.join(QADIR, f"chrome-profile-{os.getpid()}")
BASE = "http://localhost:8000"
os.makedirs(QADIR, exist_ok=True)

proc = subprocess.Popen(
    [CHROME, "--headless=new", "--disable-gpu", "--no-sandbox",
     "--hide-scrollbars", "--remote-allow-origins=*",
     f"--remote-debugging-port={PORT}",
     f"--user-data-dir={USERDIR}", "--window-size=1500,2600", "about:blank"],
    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

for _ in range(100):
    try:
        with urllib.request.urlopen(f"http://localhost:{PORT}/json/version", timeout=1) as r:
            json.load(r); break
    except Exception:
        time.sleep(0.2)
with urllib.request.urlopen(f"http://localhost:{PORT}/json/list") as r:
    page = next(t for t in json.load(r) if t.get("type") == "page")
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

def drain(seconds):
    end = time.time() + seconds
    ws.settimeout(0.4)
    while time.time() < end:
        try:
            ws.recv()
        except Exception:
            pass
    ws.settimeout(None)

send("Page.enable"); send("Runtime.enable")

def evaluate(expr):
    r = send("Runtime.evaluate", expression=expr, awaitPromise=True, returnByValue=True)
    return r.get("result", {}).get("result", {}).get("value")

def navigate(url, settle=6.0):
    send("Page.navigate", url="about:blank"); drain(0.6)
    send("Page.navigate", url=url); drain(settle)

report = {}

# 1. load F0 deep-link (what the opened window started from)
navigate(f"{BASE}/#pts=-1,0;1,0;0,1;0,-1", settle=6)
# 2. click the dP3 preset (diagram changes; hash must follow now)
evaluate("[...document.querySelectorAll('button.preset[data-pts]')].find(b=>b.dataset.pts==='1,0;0,1;-1,1;-1,0;0,-1;1,-1').click(); null")
drain(8)
report["after_preset"] = {
    "hash": evaluate("location.hash"),
    "num_gauge": evaluate("(document.getElementById('iv-gauge')||{}).textContent"),
}
# 3. click a glowing square face in the dimer
evaluate("document.querySelector('#iv-tiling polygon[data-qnode]').dispatchEvent(new MouseEvent('click',{bubbles:true})); null")
drain(8)
report["after_face_click"] = {
    "hash": evaluate("location.hash"),
    "num_gauge": evaluate("(document.getElementById('iv-gauge')||{}).textContent"),
    "fields": evaluate("(document.getElementById('iv-fields')||{}).textContent"),
    "badge": evaluate("(document.getElementById('iv-badge')||{}).textContent"),
}
# 4. cold-reload the resulting hash: must reproduce the SAME state
h = report["after_face_click"]["hash"]
navigate(f"{BASE}/{h}", settle=8)
report["reloaded"] = {
    "hash": evaluate("location.hash"),
    "num_gauge": evaluate("(document.getElementById('iv-gauge')||{}).textContent"),
    "fields": evaluate("(document.getElementById('iv-fields')||{}).textContent"),
    "badge": evaluate("(document.getElementById('iv-badge')||{}).textContent"),
}
# 5. lattice-dot edit must also sync pts= (add a point to F0)
navigate(f"{BASE}/#pts=-1,0;1,0;0,1;0,-1", settle=6)
evaluate("document.querySelector('[data-gx=\"1\"][data-gy=\"-1\"]') ? (document.querySelector('[data-gx=\"1\"][data-gy=\"-1\"]').dispatchEvent(new MouseEvent('click',{bubbles:true})), 1) : 0")
drain(8)
report["after_dot_click"] = {"hash": evaluate("location.hash")}

print(json.dumps(report, indent=2))
ws.close(); proc.terminate()
