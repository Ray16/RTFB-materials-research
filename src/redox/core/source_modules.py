"""Load explicitly selected source snapshots without importing campaign orchestration.

Callers must verify provenance before executing a trusted Python snapshot.
Loading a module does not qualify its scientific results.
"""
import importlib.util


def load_source_module(path, name):
    """Execute a trusted source file under the supplied isolated module name."""
    definition = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(module)
    return module
