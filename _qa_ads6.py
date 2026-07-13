"""Headless check of the AdS6/CFT5 tab. Usage: python _qa_ads6.py <port>"""
import json, os, subprocess, sys, time, urllib.request
import websocket
CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
PORT = 9381
QADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "_qa")
USERDIR = os.path.join(QADIR, f"chrome-profile-{os.getpid()}")
BASE = f"http://localhost:{sys.argv[1] if len(sys.argv) > 1 else '8080'}"
os.makedirs(QADIR, exist_ok=True)
proc = subprocess.Popen([CHROME, "--headless=new", "--disable-gpu", "--no-sandbox",
    "--remote-allow-origins=*", f"--remote-debugging-port={PORT}",
    f"--user-data-dir={USERDIR}", "--window-size=1500,2600", "about:blank"],
    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
for _ in range(100):
    try: urllib.request.urlopen(f"http://localhost:{PORT}/json/version", timeout=1); break
    except Exception: time.sleep(0.2)
page = next(t for t in json.load(urllib.request.urlopen(f"http://localhost:{PORT}/json/list")) if t.get("type")=="page")
ws = websocket.create_connection(page["webSocketDebuggerUrl"], max_size=64*1024*1024)
_id=0
def send(m,**p):
    global _id;_id+=1;ws.send(json.dumps({"id":_id,"method":m,"params":p}))
    while True:
        r=json.loads(ws.recv())
        if r.get("id")==_id:return r
console=[]
def drain(s):
    end=time.time()+s;ws.settimeout(0.3)
    while time.time()<end:
        try:
            m=json.loads(ws.recv())
            if m.get("method")=="Runtime.consoleAPICalled":
                p=m["params"];console.append((p.get("type"),"".join(str(a.get("value","")) for a in p.get("args",[]))))
            elif m.get("method")=="Runtime.exceptionThrown":
                console.append(("EXCEPTION",m["params"]["exceptionDetails"].get("text","")+" "+(m["params"]["exceptionDetails"].get("exception",{}).get("description","") or "")))
        except Exception: pass
    ws.settimeout(None)
send("Page.enable");send("Runtime.enable")
def ev(e): return send("Runtime.evaluate",expression=e,awaitPromise=True,returnByValue=True).get("result",{}).get("result",{}).get("value")
send("Page.navigate",url=f"{BASE}/");drain(4)
ev("document.getElementById('tab-ads6').click(); null");drain(6)
rep={
 "tab_active": ev("document.getElementById('tab-ads6').classList.contains('on')"),
 "grid_dots": ev("document.querySelectorAll('#a6-grid .gridpt').length"),
 "grid_sel": ev("document.querySelectorAll('#a6-grid .gridpt.sel').length"),
 "col_rank": ev("(document.getElementById('a6-rank2')||{}).textContent"),
 "rank": ev("(document.getElementById('a6-rank')||{}).textContent"),
 "oneform": ev("(document.getElementById('a6-oneform')||{}).textContent"),
 "defect": ev("(document.getElementById('a6-defect')||{}).textContent"),
 "cubic": ev("(document.getElementById('a6-cubic')||{}).textContent"),
 "ansatz": ev("(document.getElementById('a6-ansatz')||{}).textContent"),
 "sigma": ev("(document.getElementById('a6-sigma')||{}).textContent"),
 "stacks": ev("(document.getElementById('a6-stacks')||{}).textContent"),
 "charge": ev("(document.getElementById('a6-charge')||{}).textContent"),
 "holo": ev("(document.getElementById('a6-holo')||{}).textContent"),
 "n_pole_rows": ev("document.querySelectorAll('#a6-poles tr').length"),
 "n_sigma_dots": ev("document.querySelectorAll('#a6-sigma-svg circle[fill]').length"),
 "web_junctions": ev("document.querySelectorAll('#a6-web-svg .wjunc').length"),
 "web_legs": ev("document.querySelectorAll('#a6-web-svg .wleg').length"),
 "web_leg_labels": ev("[...document.querySelectorAll('#a6-web-svg .wleglabel')].map(t=>t.textContent)"),
 "errors": [c for c in console if c[0] in ('error','EXCEPTION')],
}
# click the dP3 preset -> live recompute
ev("[...document.querySelectorAll('button.a6preset')].find(b=>/dP3/.test(b.textContent)).click(); null");drain(5)
rep["dp3_stacks"]=ev("(document.getElementById('a6-stacks')||{}).textContent")
rep["dp3_flavor"]=ev("(document.getElementById('a6-flavor')||{}).textContent")
rep["dp3_grid_sel"]=ev("document.querySelectorAll('#a6-grid .gridpt.sel').length")
# now toggle one lattice point on the grid -> live recompute to a new geometry
before_stacks = ev("(document.getElementById('a6-stacks')||{}).textContent")
ev("var d=document.querySelector('#a6-grid .gridpt:not(.sel)'); if(d) d.dispatchEvent(new MouseEvent('click',{bubbles:true})); null");drain(5)
rep["after_toggle_stacks"]=ev("(document.getElementById('a6-stacks')||{}).textContent")
rep["after_toggle_sel"]=ev("document.querySelectorAll('#a6-grid .gridpt.sel').length")
rep["errors_final"]=[c for c in console if c[0] in ('error','EXCEPTION')]
print(json.dumps(rep,indent=2,default=str))
ws.close();proc.terminate()
