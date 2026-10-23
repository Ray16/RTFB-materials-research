"""Use the installed bounded gradient contractions with the direct LU solver.

Only grad_solver's local contraction branch changes. The solvent object's
lowmem_intermediate_storage flag is left unchanged, preserving direct energies.
"""
from contextlib import contextmanager
import hashlib
import importlib.util
import inspect
from pathlib import Path
import sys
import textwrap


def contraction_source(source):
    old='lowmem_mode = getattr(pcmobj, "lowmem_intermediate_storage", False)'
    if source.count(old)!=1:
        raise ValueError('Installed solvent-gradient branch differs from the reviewed template')
    return source.replace(old,'lowmem_mode = True')


@contextmanager
def compact_smd_execution(directory,budget_GiB):
    from gpu4pyscf.solvent.grad import pcm
    from gpu4pyscf.solvent import smd
    source=Path(__file__).with_name('smd_workspace_execution.py')
    spec=importlib.util.spec_from_file_location('qualified_workspace_underlay',source)
    workspace=importlib.util.module_from_spec(spec);spec.loader.exec_module(workspace)
    with workspace.compact_smd_execution(directory,budget_GiB=budget_GiB) as proof:
        original=pcm.grad_solver;before=inspect.getsource(original);after=contraction_source(before)
        path=Path(directory)/'bounded_gradient_solver.py'
        header='\n'.join(Path(pcm.__file__).read_text().splitlines()[:13])+'\n\n'
        path.write_text(header+after)
        namespace={};exec(compile(textwrap.dedent(after),str(path),'exec'),original.__globals__,namespace)
        replacement=namespace['grad_solver']
        proof.update(gradient_solver_calls=0,direct_LU_solver_retained=True,gradient_contraction_branch_only=True)
        proof['source_sha256'].update({str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest()
            for p in (Path(__file__).resolve(),source,Path(pcm.__file__),path)})
        def solver(pcmobj,*args,**kwargs):
            if not isinstance(pcmobj,smd.SMD):return original(pcmobj,*args,**kwargs)
            if getattr(pcmobj,'lowmem_intermediate_storage',False):
                raise ValueError('Direct-LU gradient profile cannot use iterative energy storage')
            proof['gradient_solver_calls']+=1
            return replacement(pcmobj,*args,**kwargs)
        for name,module in list(sys.modules.items()):
            if name.startswith('gpu4pyscf') and module is not None:
                for key,value in list(vars(module).items()):
                    if value is original:setattr(module,key,solver)
        try:yield proof
        finally:
            for name,module in list(sys.modules.items()):
                if name.startswith('gpu4pyscf') and module is not None:
                    for key,value in list(vars(module).items()):
                        if value is solver:setattr(module,key,original)
