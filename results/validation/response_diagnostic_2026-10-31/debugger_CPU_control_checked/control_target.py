import ctypes,resource
from pathlib import Path
from redox.workflows.owned_debugger import permit
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
permit(Path('/nfs/lambda_stor_01/homes/rzhu/0_redox/results/validation/response_diagnostic_2026-10-31/debugger_CPU_control_checked'),'a6222ff0f450ca040e3479c04fdb7d0f1e0cabb099688270e720b7115d933c3d')
ctypes.string_at(0)
