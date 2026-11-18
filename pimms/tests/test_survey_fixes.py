"""Fixes from the October 2026 survey of the engine's Python layer."""

import importlib
import sys

import pytest

import pimms


def test_a_cluster_kernel_that_will_not_import_is_announced(monkeypatch):
    """The pure-Python fallback is identical but many times slower, so a kernel
    that exists and fails to import (a partial rebuild) must not be hidden."""
    pytest.importorskip("pimms.cluster_kernels")
    from pimms import cluster_utils

    assert cluster_utils._HAVE_CLUSTER_KERNELS is True
    # `from . import cluster_kernels` returns the package attribute when there is
    # one, so remove it; None in sys.modules then makes the import raise ImportError
    monkeypatch.delattr(pimms, "cluster_kernels")
    monkeypatch.setitem(sys.modules, "pimms.cluster_kernels", None)
    try:
        with pytest.warns(UserWarning, match="cluster_kernels could not be imported"):
            importlib.reload(cluster_utils)
        assert cluster_utils._HAVE_CLUSTER_KERNELS is False
    finally:
        monkeypatch.undo()
        importlib.reload(cluster_utils)
    assert cluster_utils._HAVE_CLUSTER_KERNELS is True
