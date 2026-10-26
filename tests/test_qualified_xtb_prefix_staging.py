"""Safety checks for the local relocation-test prefix stager."""
import importlib.util
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/release/prepare_qualified_xtb_prefix.py'
spec = importlib.util.spec_from_file_location('prefix_stager', SCRIPT)
stager = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stager)


def test_ldd_distinguishes_source_and_system_and_resolves_links(tmp_path):
    source = tmp_path / 'source'
    source.mkdir()
    lib = source / 'liba.so.1'
    lib.write_bytes(b'qualified dependency')
    alias = source / 'liba.so'
    alias.symlink_to(lib)
    external = tmp_path / 'libsystem.so'
    external.write_bytes(b'system dependency')
    deps, system = stager.parse_ldd(
        f'liba.so => {alias} (0x1)\nlibsystem.so => {external} (0x2)', source
    )
    assert deps == {'liba.so': lib}
    assert system == {'libsystem.so': external}


def test_ldd_rejects_unresolved_and_conflicting_libraries(tmp_path):
    with pytest.raises(ValueError, match='Unresolved'):
        stager.parse_ldd('libmissing.so => not found', tmp_path)
    one, two = tmp_path / 'one', tmp_path / 'two'
    one.touch(); two.touch()
    with pytest.raises(ValueError, match='Conflicting'):
        stager.parse_ldd(f'lib.so => {one} (0x1)\nlib.so => {two} (0x2)', tmp_path)


def test_staging_refuses_existing_or_inside_source(tmp_path):
    source = tmp_path / 'source'
    source.mkdir()
    with pytest.raises(ValueError, match='outside'):
        stager.prepare(source, source / 'staged')
    with pytest.raises(FileExistsError, match='Preserve'):
        stager.prepare(source, tmp_path)
    assert list(source.iterdir()) == []


def test_relocated_probe_requires_no_site_and_rejects_original_project(tmp_path):
    script = SCRIPT.with_name('probe_relocated_xtb.py')
    args = [str(script), '--project', str(tmp_path), '--forbid-root', str(tmp_path),
            '--xtb-prefix', str(tmp_path), '--site-packages', str(tmp_path),
            '--out', str(tmp_path / 'output.json')]
    default = subprocess.run([sys.executable, *args], capture_output=True, text=True)
    assert default.returncode != 0 and 'with -S' in default.stderr
    isolated = subprocess.run([sys.executable, '-S', *args], capture_output=True, text=True)
    assert isolated.returncode != 0 and 'outside forbidden root' in isolated.stderr
    assert not (tmp_path / 'output.json').exists()
