"""Cooperative, PID-specific debugging of a new owned diagnostic process.

The target opts in to one verified debugger through the Linux Yama API. Shared
ptrace policy stays unchanged; existing scientific workers are never attached.
"""
import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from redox.core.io import atomic_write_json


def sha256_file(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def process_start(pid):return Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()[19]


def permit(case,executable_hash,timeout=60):
    pid=os.getpid();owner=dict(pid=pid,process_start=process_start(pid),host=os.uname().nodename,
        executable=sys.executable,executable_sha256=sha256_file('/proc/self/exe'))
    atomic_write_json(case/'debugger_wait.json',owner);deadline=time.monotonic()+timeout
    path=case/'debugger_owner.json'
    while not path.exists():
        if time.monotonic()>deadline:raise TimeoutError('Owned debugger did not arrive')
        time.sleep(.1)
    debugger=json.loads(path.read_text());proc=Path('/proc')/str(debugger['pid'])
    # A sandboxed inferior cannot inspect the debugger's host-namespace exe.
    # Its controller verifies that executable before publishing this receipt;
    # the inferior checks PID/start plus the bound executable/script identities.
    if (proc.stat().st_uid!=os.getuid() or process_start(debugger['pid'])!=debugger['process_start']
            or debugger.get('executable_sha256')!=executable_hash
            or debugger.get('controller_verified_executable') is not True
            or sha256_file(case/'debugger.gdb')!=debugger.get('script_sha256')):
        raise ValueError('Debugger PID/start, executable or exact command identity differs')
    libc=ctypes.CDLL(None,use_errno=True);prctl=libc.prctl
    prctl.argtypes=[ctypes.c_int,ctypes.c_ulong,ctypes.c_ulong,ctypes.c_ulong,ctypes.c_ulong];prctl.restype=ctypes.c_int
    dumpable_before=prctl(3,0,0,0,0)
    if prctl(4,1,0,0,0):raise OSError(ctypes.get_errno(),'Owned diagnostic dumpability setup failed')
    if prctl(0x59616d61,debugger['pid'],0,0,0):raise OSError(ctypes.get_errno(),'PID-specific debugger permission failed')
    try:
        atomic_write_json(case/'debugger_permitted.json',dict(target=owner,debugger=debugger,dumpable_before=dumpable_before,global_ptrace_policy_changed=False))
        while not (case/'debugger_attached.json').exists():
            if time.monotonic()>deadline:raise TimeoutError('Owned debugger did not attach')
            time.sleep(.1)
    finally:prctl(0x59616d61,0,0,0,0)


def configuration(executable):
    executable=executable.resolve();prefix=executable.parent.parent
    data=prefix/'share/gdb';library=prefix/'lib/x86_64-linux-gnu';env=dict(os.environ,LD_LIBRARY_PATH=str(library))
    run=subprocess.run([str(executable),'--nx','--batch','--data-directory='+str(data),'--version'],env=env,capture_output=True,text=True,check=True,timeout=15)
    return dict(executable=str(executable),sha256=sha256_file(executable),data_directory=str(data),library_path=str(library),version=run.stdout.splitlines()[0])


def start(case,config,command_fragment):
    owner=json.loads((case/'debugger_wait.json').read_text());pid=owner['pid'];proc=Path('/proc')/str(pid)
    if (proc.stat().st_uid!=os.getuid() or process_start(pid)!=owner['process_start']
            or command_fragment.encode() not in (proc/'cmdline').read_bytes()):raise ValueError('Only this exact waiting diagnostic can be debugged')
    executable=Path(config['executable'])
    if sha256_file(executable)!=config['sha256']:raise ValueError('Debugger executable changed')
    target_executable=Path(owner['executable']).resolve()
    if sha256_file(target_executable)!=owner['executable_sha256']:raise ValueError('Host and namespace Python executable identities differ')
    commands=('set pagination off\nset confirm off\nset debuginfod enabled off\nset auto-load safe-path /dev/null\nset print thread-events off\nset sysroot /\n'
        'handle SIGTERM nostop noprint pass\nhandle SIGINT nostop noprint pass\nhandle SIGSEGV stop print nopass\nhandle SIGBUS stop print nopass\npython\n'
        'import gdb,json,pathlib,time\ncase=pathlib.Path('+repr(str(case))+')\ndeadline=time.monotonic()+60\n'
        'while not (case/"debugger_permitted.json").exists():\n'
        ' if time.monotonic()>deadline:raise RuntimeError("Target did not grant debugger permission")\n time.sleep(.1)\n'
        'gdb.execute("file "+'+repr(str(target_executable))+')\ngdb.execute("attach '+str(pid)+'")\n'
        '(case/"debugger_attached.json").write_text(json.dumps({"pid":'+str(pid)+'}))\n'
        'last_signal=[]\n'
        'def stopped(event):\n if isinstance(event,gdb.SignalEvent):last_signal.append(event.stop_signal)\n'
        'gdb.events.stop.connect(stopped)\ngdb.execute("continue")\n'
        'active=gdb.selected_inferior().pid\n'
        'if active:\n gdb.execute("thread apply all bt 40")\n gdb.execute("info sharedlibrary")\n'
        '(case/"debugger_stop.json").write_text(json.dumps({"signals":last_signal,"inferior_pid":active,"native_G_qualified":False}))\n'
        'if active:gdb.execute("kill")\nend\nquit\n')
    script=case/'debugger.gdb';script.write_text(commands)
    env=dict(os.environ,LD_LIBRARY_PATH=config['library_path'],DEBUGINFOD_URLS='')
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
    with (case/'debugger.log').open('x') as log:
        child=subprocess.Popen([str(executable),'--nx','--batch','--data-directory='+config['data_directory'],'--command',str(script)],env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    if sha256_file(Path('/proc')/str(child.pid)/'exe')!=config['sha256']:
        child.terminate();child.wait(timeout=3);raise ValueError('Spawned debugger executable differs')
    atomic_write_json(case/'debugger_owner.json',dict(pid=child.pid,process_start=process_start(child.pid),executable_sha256=config['sha256'],controller_verified_executable=True,script_sha256=sha256_file(script),target=owner))
    return child


def control(executable,output,namespace_command=None,helper_path=None,exercise='null_pointer'):
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=2);os.sched_setaffinity(0,cpu_allocation(1))
    if output.exists():raise FileExistsError('Preserve earlier debugger controls')
    output.mkdir(parents=True);config=configuration(executable);target=output/'control_target.py'
    importer=('from redox.workflows.owned_debugger import permit\n' if helper_path is None else
        'import importlib.util\nspec=importlib.util.spec_from_file_location("checked_debug_helper",'+repr(str(helper_path))+')\nhelper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper);permit=helper.permit\n')
    if exercise=='null_pointer':body='ctypes.string_at(0)\n'
    elif exercise in ('C_watchdog','GIL_held'):
        body=('import faulthandler,threading,time\nstream=open("/dev/null","w")\nstopped=threading.Event()\n'
            'def sample():\n while not stopped.wait(.001):faulthandler.dump_traceback(file=stream,all_threads=True)\n'
            +('faulthandler.dump_traceback_later(.001,repeat=True,file=stream)\n' if exercise=='C_watchdog' else 'sampler=threading.Thread(target=sample);sampler.start()\n')+
            'for i in range(5000):\n thread=threading.Thread(target=lambda:None);thread.start();thread.join()\n'
            +('faulthandler.cancel_dump_traceback_later()\n' if exercise=='C_watchdog' else 'stopped.set();sampler.join()\n')+'stream.close()\n')
    else:raise ValueError('Unknown owned CPU debugger exercise')
    target.write_text('import ctypes,resource\nfrom pathlib import Path\n'+importer+
        'resource.setrlimit(resource.RLIMIT_CORE,(0,0))\npermit(Path('+repr(str(output))+'),'+repr(config['sha256'])+')\n'+body)
    with (output/'target.log').open('x') as log:child=subprocess.Popen([*(namespace_command or []),sys.executable,str(target)],stdout=log,stderr=subprocess.STDOUT)
    debugger=None
    try:
        deadline=time.monotonic()+20
        while not (output/'debugger_wait.json').exists():
            if child.poll() is not None or time.monotonic()>deadline:raise RuntimeError('Owned CPU debug target did not start')
            time.sleep(.1)
        debugger=start(output,config,str(target));debugger_code=debugger.wait(timeout=60);target_code=child.wait(timeout=10)
    finally:
        for process in (debugger,child):
            if process is not None and process.poll() is None:
                process.terminate()
                try:process.wait(timeout=3)
                except subprocess.TimeoutExpired:process.kill();process.wait(timeout=3)
    stop=json.loads((output/'debugger_stop.json').read_text())
    expected_codes=(-9,137) if namespace_command is not None else (-9,)
    fault=stop['signals']==['SIGSEGV'] and target_code in expected_codes
    normal=not stop['signals'] and target_code==0
    if debugger_code!=0 or not (fault or (normal and exercise!='null_pointer')):raise ValueError('Owned CPU exercise did not capture a fault or normal completion')
    result=dict(status='owned_CPU_debugger_exercise_checked',exercise=exercise,debugger=config,debugger_returncode=debugger_code,target_returncode=target_code,signals=stop['signals'],intentional_test_fault=exercise=='null_pointer',read_only_namespace_tested=namespace_command is not None,GPU_used=False,global_ptrace_policy_changed=False,native_G_qualified=False)
    atomic_write_json(output/'control.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--debugger',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--exercise',choices=('null_pointer','C_watchdog','GIL_held'),default='null_pointer');a=p.parse_args()
    result=control(a.debugger.resolve(),a.out.resolve(),exercise=a.exercise);print(json.dumps(result))
