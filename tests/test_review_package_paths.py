"""A relocated review package must never follow paths back outside itself."""
import importlib.util
import io
from pathlib import Path
import tarfile

import pytest


ROOT = Path(__file__).resolve().parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts/release" / (name + ".py"))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


@pytest.mark.parametrize("relative", ["../outside.json", "/outside.json"])
def test_replayer_refuses_external_paths(tmp_path, relative):
    with pytest.raises(ValueError, match="escapes"):
        module("reproduce_review_package").local(tmp_path, relative)


def test_replayer_refuses_symlink_escape(tmp_path):
    (tmp_path / "link").symlink_to(tmp_path.parent, target_is_directory=True)
    with pytest.raises(ValueError, match="escapes"):
        module("reproduce_review_package").local(tmp_path, "link/unrelated")


def archive(name, symlink=False):
    data = io.BytesIO()
    with tarfile.open(fileobj=data, mode="w") as stream:
        entry = tarfile.TarInfo(name)
        if symlink:
            entry.type = tarfile.SYMTYPE
            entry.linkname = "/unrelated"
            stream.addfile(entry)
        else:
            entry.size = 4
            stream.addfile(entry, io.BytesIO(b"text"))
    return data.getvalue()


@pytest.mark.parametrize("name,symlink", [("../escape", False), ("/escape", False), ("link", True)])
def test_code_export_refuses_unsafe_members(tmp_path, name, symlink):
    with pytest.raises(ValueError, match="Unsafe"):
        module("prepare_review_package").extract_source(archive(name, symlink), tmp_path)


def test_code_export_preserves_regular_file_bytes(tmp_path):
    module("prepare_review_package").extract_source(archive("src/example.py"), tmp_path)
    assert (tmp_path / "src/example.py").read_bytes() == b"text"
