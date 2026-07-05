"""Throwaway QA: verify blow-up/blow-down labelling directions.

P^2 (&blow) -> click a + handle -> must identify dP1 (grow, 3->4 gauge).
dP1 (&blowdown) -> click a - handle -> must identify P^2/dP0 (shrink).
"""
import base64, json, os, subprocess, time, urllib.request
import websocket

CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
PORT = 9351
QADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_qa")
USERDIR = os.path.join(QADIR, f"chrome-profile-{os.getpid()}")
BASE = "http://localhost:8765"
os.makedirs(QADIR, exist_ok=True)

proc = subprocess.Popen(
    [CHROME, "--headless=new", "--disable-gpu", "--no-sandbox",
     "--hide-scrollbars", "--remote-allow-origins=*",
     f"--remote-debugging-port={PORT}",
     f"--user-data-dir={USERDIR}", "--window-size=1500,1400",
     "about:blank"],
    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

for _ in range(100):
    try:
        with urllib.request.urlopen(f"http://localhost:{PORT}/json/version", timeout=1) as r:
            json.load(r); break
    except Exception:
        time.sleep(0.2)
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

def navigate(url, settle=3.5):
    send("Page.navigate", url="about:blank"); drain(0.5)
    send("Page.navigate", url=url); drain(settle)

def shot(name):
    r = send("Page.captureScreenshot", format="png", captureBeyondViewport=True)
    with open(os.path.join(QADIR, name), "wb") as f:
        f.write(base64.b64decode(r["result"]["data"]))

def state():
    return {
        "identified": evaluate("(document.getElementById('t-geom')||{}).textContent"
                               " || Array.from(document.querySelectorAll('td,div'))"
                               ".map(e=>e.textContent).find(t=>/Cone over|del Pezzo|P\\^2|projective/i.test(t)) || null"),
        "gauge": evaluate("(Array.from(document.querySelectorAll('tr')).map(r=>r.textContent)"
                          ".find(t=>/gauge groups/i.test(t))||'').trim()"),
        "pts_hash": evaluate("location.hash"),
        "n_handles_up": evaluate("document.querySelectorAll('.blowhandle:not(.down)').length"),
        "n_handles_down": evaluate("document.querySelectorAll('.blowhandle.down').length"),
    }

report = {}

# --- 1. P^2 in blow-UP mode: click the first + handle -> expect dP1 ---
navigate(f"{BASE}/#pts=1,0;0,1;-1,-1&blow", settle=4)
report["P2_before"] = state()
shot("blow_P2_before.png")
evaluate("var h=document.querySelector('.blowhandle:not(.down)');"
         "h && h.dispatchEvent(new MouseEvent('click',{bubbles:true})); null")
drain(3.5)
report["P2_after_blowup"] = state()
shot("blow_P2_after_up.png")

# --- 2. dP1 in blow-DOWN mode: click a - handle -> expect P^2 ---
navigate(f"{BASE}/#pts=1,0;0,1;-1,-1;0,-1&blowdown", settle=4)
report["dP1_before"] = state()
shot("blow_dP1_before.png")
evaluate("var h=document.querySelector('.blowhandle.down');"
         "h && h.dispatchEvent(new MouseEvent('click',{bubbles:true})); null")
drain(3.5)
report["dP1_after_blowdown"] = state()
shot("blow_dP1_after_down.png")

print(json.dumps(report, indent=2))
ws.close(); proc.kill()
