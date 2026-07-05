"""Throwaway QA: interactive resolution blow-up on the Z3xZ3 singular cone.

Load the corners-only (singular) state, click one glowing + on the interior
divisor, and check the subdivision refines + residual labels update.
"""
import base64, json, os, subprocess, time, urllib.request
import websocket

CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
PORT = 9353
QADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_qa")
USERDIR = os.path.join(QADIR, f"chrome-profile-{os.getpid()}")
BASE = "http://localhost:8765"
os.makedirs(QADIR, exist_ok=True)

proc = subprocess.Popen(
    [CHROME, "--headless=new", "--disable-gpu", "--no-sandbox",
     "--hide-scrollbars", "--remote-allow-origins=*",
     f"--remote-debugging-port={PORT}",
     f"--user-data-dir={USERDIR}", "--window-size=1400,1300", "about:blank"],
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

def drain(sec):
    end = time.time() + sec
    ws.settimeout(0.4)
    while time.time() < end:
        try: ws.recv()
        except Exception: pass
    ws.settimeout(None)

send("Page.enable"); send("Runtime.enable")
def ev(expr):
    r = send("Runtime.evaluate", expression=expr, awaitPromise=True, returnByValue=True)
    return r.get("result", {}).get("result", {}).get("value")

def shot(name):
    r = send("Page.captureScreenshot", format="png", captureBeyondViewport=True)
    with open(os.path.join(QADIR, name), "wb") as f:
        f.write(base64.b64decode(r["result"]["data"]))

send("Page.navigate", url=f"{BASE}/#pts=0,0;3,0;0,3&active=0,0;3,0;0,3&blow")
drain(4)

def state():
    return {
        "cells": ev("(document.getElementById('t-tris')||{}).textContent"),
        "resid": ev("(document.getElementById('t-resid')||{}).textContent"),
        "hash": ev("location.hash"),
        "n_resolve_handles": ev("document.querySelectorAll('[data-resolve]').length"),
        "n_sing_cells": ev("document.querySelectorAll('.singcell').length"),
    }

out = {"cone": state()}
# click the + on the interior divisor (1,1)
ev("var h=document.querySelector('[data-resolve=\"1,1\"]');"
   "h && h.dispatchEvent(new MouseEvent('click',{bubbles:true})); null")
drain(3)
out["after_resolve_1_1"] = state()
shot("resolve_after_1_1.png")

# blow it back down (switch to down mode, click the - on 1,1)
ev("setBlowMode('down'); null")
drain(2.5)
out["down_mode_handles"] = ev("document.querySelectorAll('[data-unresolve]').length")
ev("var h=document.querySelector('[data-unresolve=\"1,1\"]');"
   "h && h.dispatchEvent(new MouseEvent('click',{bubbles:true})); null")
drain(3)
out["after_unresolve_1_1"] = state()

print(json.dumps(out, indent=2))
ws.close(); proc.kill()
