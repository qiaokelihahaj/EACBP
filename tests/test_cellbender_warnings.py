from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from eacbp.auditor.advanced_qc import _warning_strings as auditor_warning_strings
from eacbp.capabilities.advanced_qc import _cellbender_quality_warnings
from scripts.verify_cellbender_runtime import _warning_strings as runtime_warning_strings


def _warnings(report: str, *, stdout: str = "", stderr: str = "", log: str = "") -> list[str]:
    return _cellbender_quality_warnings(
        stdout=stdout,
        stderr=stderr,
        log_text=log,
        report_text=report,
        sidecar_metrics={"convergence_indicator": 1.233},
    )


def test_cellbender_report_parser_ignores_css_and_keeps_visible_quality_items():
    report = """
    <html><head>
      <style>.warn { color: var(--jp-warn-color0); } /* WARNING in CSS */</style>
      <script>const warning = 'convergence warning from script';</script>
      <title>warning should not be read from head</title>
    </head><body>
      <p>Generally it is desirable for the ELBO to converge at a stable plateau.</p>
      <p><strong>WARNING</strong>: The training ELBO deviates from the max value.</p>
      <ul>
        <li>The test ELBO ends low and the output could be suboptimal.</li>
      </ul>
    </body></html>
    """

    warnings = _warnings(
        report,
        stderr="wrapper WARNING: stderr warning retained",
        log="CellBender WARNING: log warning retained",
    )

    assert any(item.startswith("stderr:") for item in warnings)
    assert any(item.startswith("log:") for item in warnings)
    assert any("report: WARNING: The training ELBO deviates" in item for item in warnings)
    assert any("report:" in item and "suboptimal" in item for item in warnings)
    assert not any("warn-color" in item or "script" in item or "head" in item for item in warnings)
    assert not any("desirable for the ELBO to converge" in item for item in warnings)
    assert any(item.startswith("metrics: convergence_indicator=") for item in warnings)


def test_cellbender_real_report_keeps_elbo_warnings_beyond_jupyter_css():
    report_path = Path("outputs/cellbender_smoke/tiny_output_report.html")
    if not report_path.is_file():
        pytest.skip("CellBender smoke report is absent; run the CellBender smoke workflow first")
    warnings = _warnings(report_path.read_text(encoding="utf-8"))

    assert len(warnings) <= 32
    assert any("report: WARNING: The training ELBO deviates" in item for item in warnings)
    assert any("report:" in item and "wrong direction" in item for item in warnings)
    assert any("report:" in item and "suboptimal" in item for item in warnings)
    assert any(item.startswith("metrics: convergence_indicator=") for item in warnings)
    assert not any("jp-warn-color" in item or "warn icon colors" in item for item in warnings)
    assert not any("Assessing convergence" in item or "lack of convergence" in item for item in warnings)


def test_cellbender_warning_arrays_are_flattened_item_by_item():
    value = np.asarray([["stderr: first", "report: second"]], dtype=object)

    assert runtime_warning_strings(value) == ["stderr: first", "report: second"]
    assert auditor_warning_strings(value) == ["stderr: first", "report: second"]
    assert runtime_warning_strings(np.asarray("single", dtype=object)) == ["single"]
