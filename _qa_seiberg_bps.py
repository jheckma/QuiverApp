"""Throwaway QA: multi-duality Seiberg clicks + adjoint refusal + BPS tab.
Headless Chrome over CDP (same pattern as _qa_dimerclick.py).
Usage: start the server, then `python _qa_seiberg_bps.py <port>`."""
import base64, json, os, subprocess, sys, time, urllib.request
import websocket

CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
PORT = 9361
QADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_qa")
USERDIR = os.path.join(QADIR, f"chrome-profile-{os.getpid()}")
BASE = f"http://localhost:{sys.argv[1] if len(sys.argv) > 1 else '8000'}"
os.makedirs(QADIR, exist_ok=True)

proc = subprocess.Popen(
    [CHROME, "--headless=new", "--disable-gpu", "--no-sandbox",
     "--hide-scrollbars", "--remote-allow-origins=*",
     f"--remote-debugging-port={PORT}",
     f"--user-data-dir={USERDIR}", "--window-size=1500,2800", "about:blank"],
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
            txt = " ".join(str(a.get("value", a.get("description", ""))) for a in p.get("args", []))
            console.append((p.get("type", "log"), txt))
        elif m == "Runtime.exceptionThrown":
            d = msg["params"]["exceptionDetails"]
            console.append(("EXCEPTION", d.get("text", "") + " " +
                            (d.get("exception", {}).get("description", "") or "")))
    ws.settimeout(None)


send("Page.enable"); send("Runtime.enable"); send("Log.enable")

def ev(expr):
    r = send("Runtime.evaluate", expression=expr, awaitPromise=True, returnByValue=True)
    return r.get("result", {}).get("result", {}).get("value")

def navigate(url, settle=6.0):
    send("Page.navigate", url="about:blank"); drain(0.6)
    send("Page.navigate", url=url); drain(settle)

def click_qnode(node, settle=6):
    ev(f"(document.querySelector('#iv-quiver [data-qnode=\"{node}\"]')||{{dispatchEvent(){{}}}}).dispatchEvent(new MouseEvent('click',{{bubbles:true}})); null")
    drain(settle)

def errs():
    return [c for c in console if c[0] in ("error", "EXCEPTION")]

report = {}

# ===== 1. F0: multiple successive dualities on quiver nodes =====
console.clear()
navigate(f"{BASE}/#pts=-1,0;1,0;0,1;0,-1", settle=7)
seq = []
for k in [0, 1, 2, 0]:
    console.clear()
    click_qnode(k)
    seq.append({
        "clicked": k,
        "fields": ev("(document.getElementById('iv-fields')||{}).textContent"),
        "badge": ev("(document.getElementById('iv-badge')||{}).textContent"),
        "path": ev("(document.getElementById('iv-dual-path')||{}).textContent"),
        "errors": errs(),
    })
report["f0_multi_duality"] = seq

# ===== 2. Adjoint refusal: L131 node 2 (seed adjoint) should be refused =====
console.clear()
navigate(f"{BASE}/#pts=0,0;1,0;2,1;0,1", settle=7)
before_fields = ev("(document.getElementById('iv-fields')||{}).textContent")
console.clear()
click_qnode(2)                      # node 2 has an adjoint -> refuse, re-render valid
report["adjoint_refusal"] = {
    "before_fields": before_fields,
    "after_fields": ev("(document.getElementById('iv-fields')||{}).textContent"),
    "blowhint": ev("(document.getElementById('blowhint')||{}).textContent"),
    "iv_quiver_has_nodes": ev("document.querySelectorAll('#iv-quiver [data-qnode]').length"),
    "errors": errs(),
}

# ===== 3. BPS tab: canonical chamber shows for dP0 =====
console.clear()
navigate(f"{BASE}/#bps", settle=7)
drain(3)
report["bps_tab_dp0"] = {
    "stable_count": ev("(document.getElementById('bps-stable-count')||{}).textContent"),
    "canon": ev("(document.getElementById('bps-canon')||{}).textContent"),
    "kk": ev("(document.getElementById('bps-kk')||{}).textContent"),
    "flavor": ev("(document.getElementById('bps-flavor')||{}).textContent"),
    "mgs": ev("(document.getElementById('bps-mgs')||{}).textContent"),
    "sn_legs": ev("(document.getElementById('bps-sn-legs')||{}).textContent"),
    "sn_rank": ev("(document.getElementById('bps-sn-rank')||{}).textContent"),
    "errors": errs(),
}
# click a BPS quiver node -> appends a mutation and recomputes
console.clear()
ev("(document.querySelector('#bps-quiver-svg [data-qnode]')||{dispatchEvent(){}}).dispatchEvent(new MouseEvent('click',{bubbles:true})); null")
drain(5)
report["bps_after_click"] = {
    "seq": ev("(document.getElementById('bps-seq')||{}).textContent"),
    "stable_count": ev("(document.getElementById('bps-stable-count')||{}).textContent"),
    "errors": errs(),
}

print(json.dumps(report, indent=2, default=str))
ws.close(); proc.terminate()
