"""Checkpoint reuse must preserve exact bytes and reject unsafe storage layouts."""
import h5py
import numpy as np
import pytest

from redox.backends.checkpoint_arrays import contiguous_payload


@pytest.mark.parametrize('shape', [(), (2, 3)])
def test_exact_readonly_payload(tmp_path, shape):
    path = tmp_path / 'checkpoint.h5'
    expected = np.arange(np.prod(shape) if shape else 1, dtype=np.float64).reshape(shape)
    with h5py.File(path, 'w') as handle:
        handle['value'] = expected
    content = path.read_bytes()
    with h5py.File(path, 'r') as handle:
        actual = contiguous_payload(handle, 'value', content)
        assert actual.tobytes() == expected.tobytes()
        assert not actual.flags.writeable
        with pytest.raises(ValueError):
            contiguous_payload(handle, 'value', content[:8])


@pytest.mark.parametrize('kwargs', [dict(dtype='f4'), dict(chunks=True), dict(compression='gzip')])
def test_reject_unsupported_layout(tmp_path, kwargs):
    path = tmp_path / 'checkpoint.h5'
    with h5py.File(path, 'w') as handle:
        handle.create_dataset('value', data=[1., 2.], **kwargs)
    with h5py.File(path, 'r') as handle:
        with pytest.raises(ValueError):
            contiguous_payload(handle, 'value', path.read_bytes())
