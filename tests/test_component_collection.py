import importlib.util
import json
from pathlib import Path
import pytest
pytest.importorskip('h5py')
from redox.core.hashing import sha256_file


def test_cached_review_cannot_qualify_another_namespace_or_state(tmp_path):
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('component_collector_test',root/'scripts/validation/collect_component_screening.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    namespace=tmp_path/'namespace.json';namespace.write_text('{}')
    source=dict(id='one',charge=0,mult=1);manifest=tmp_path/'manifest.json'
    manifest.write_text(json.dumps(dict(target_states=1,states=[dict(key='one',namespace=str(namespace),namespace_sha256=sha256_file(namespace),source_state=source,nao=1)])))
    output=tmp_path/'reviews';output.mkdir()
    row=dict(status='benchmark_screening_energy_control_independently_checked',key='one',namespace=str(tmp_path/'other.json'),source_state=source,source_sha256={})
    (output/'one.json').write_text(json.dumps(row))
    with pytest.raises(ValueError,match='Wrong state'):module.collect(manifest,output)
    row['namespace']=str(namespace);row['source_state']=dict(source,charge=1)
    (output/'one.json').write_text(json.dumps(row))
    with pytest.raises(ValueError,match='Wrong state'):module.collect(manifest,output)
