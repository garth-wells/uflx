"""Check HIP process isolation without requiring a GPU or MLIR bindings."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from uflx_mlir.gpu_driver import launch_isolated


def test_driver_import_does_not_load_mlir():
    """The worker entry point must not import the compiler, even indirectly."""
    driver = Path(__file__).parents[1] / "uflx_mlir" / "gpu_driver.py"
    code = """
import json, runpy, sys
runpy.run_path(sys.argv[1], run_name='driver_import_test')
print(json.dumps([name for name in sys.modules
                  if name == 'mlir' or name.startswith(('mlir.', 'uflx_mlir'))]))
"""
    result = subprocess.run(
        [sys.executable, "-I", "-c", code, str(driver)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert json.loads(result.stdout) == []


def test_worker_failure_preserves_output(tmp_path):
    """A failed child reports its error without changing the parent's output."""
    output = np.array([1.25, -2.5])
    arrays = [output, np.ones(6), np.ones(4), np.zeros(4, dtype=np.int32)]
    with pytest.raises(RuntimeError, match="HIP worker exited") as error:
        launch_isolated(
            str(tmp_path / "missing-hip-library.so"),
            b"invalid code",
            "test_kernel",
            arrays,
            (1, 1, 1),
            (4, 1, 1),
            "amd",
            0,
        )
    assert "missing-hip-library.so" in str(error.value)
    np.testing.assert_array_equal(output, [1.25, -2.5])
