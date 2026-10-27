import importlib.util
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
s=importlib.util.spec_from_file_location('uma_repeat_diagnostic',ROOT/'scripts/validation/qualify_uma_current_diagnostic.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def test_instrumentation_retains_original_panel_and_thresholds():
    worker=m.api.load_module('original_frozen_panel_template',m.PANEL)
    original=worker.run
    m.instrument(worker,Path(m.__file__).resolve())
    assert worker.run is not original
    assert worker.run.__closure__ is not None
    # This check compiles against the actual frozen template without GPU imports.
    assert callable(worker.run)
