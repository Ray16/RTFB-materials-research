import ctypes,resource
from pathlib import Path
from redox.workflows.owned_debugger import permit
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
permit(Path('/nfs/lambda_stor_01/homes/rzhu/0_redox/results/validation/response_diagnostic_2026-10-31/GIL_held_control'),'a6222ff0f450ca040e3479c04fdb7d0f1e0cabb099688270e720b7115d933c3d')
import faulthandler,threading,time
stream=open("/dev/null","w")
stopped=threading.Event()
def sample():
 while not stopped.wait(.001):faulthandler.dump_traceback(file=stream,all_threads=True)
sampler=threading.Thread(target=sample);sampler.start()
for i in range(5000):
 thread=threading.Thread(target=lambda:None);thread.start();thread.join()
stopped.set();sampler.join()
stream.close()
