set pagination off
set confirm off
set debuginfod enabled off
set auto-load safe-path /dev/null
set print thread-events off
set sysroot /
handle SIGTERM nostop noprint pass
handle SIGINT nostop noprint pass
handle SIGSEGV stop print nopass
handle SIGBUS stop print nopass
python
import gdb,json,pathlib,time
case=pathlib.Path('/nfs/lambda_stor_01/homes/rzhu/0_redox/results/validation/response_diagnostic_2026-10-31/C_watchdog_control')
deadline=time.monotonic()+60
while not (case/"debugger_permitted.json").exists():
 if time.monotonic()>deadline:raise RuntimeError("Target did not grant debugger permission")
 time.sleep(.1)
gdb.execute("file "+'/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python3.11')
gdb.execute("attach 1390721")
(case/"debugger_attached.json").write_text(json.dumps({"pid":1390721}))
last_signal=[]
def stopped(event):
 if isinstance(event,gdb.SignalEvent):last_signal.append(event.stop_signal)
gdb.events.stop.connect(stopped)
gdb.execute("continue")
active=gdb.selected_inferior().pid
if active:
 gdb.execute("thread apply all bt 40")
 gdb.execute("info sharedlibrary")
(case/"debugger_stop.json").write_text(json.dumps({"signals":last_signal,"inferior_pid":active,"native_G_qualified":False}))
if active:gdb.execute("kill")
end
quit
