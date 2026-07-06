# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Machine-readable summary emitted by the ``compare-*.py`` validation scripts.

The scripts print a human-readable per-quantity error report to stdout as before.
When called with ``--report-json PATH`` they *additionally* write the same
information as a small JSON blob to ``PATH``. ``verify.py`` reads that blob so it
can show the achieved error for every comparison in its summary table -- even for
comparisons that pass -- without scraping stdout.

Writing the report never changes the exit code: pass/fail is still signalled the
old way (exit 0 / exit 1), so the plain ``compare-*.py`` CLI stays backwards
compatible.

Schema (version 1)::

    {
      "schema":     1,
      "category":   "energy",          # volume|fault|surface|energy|receiver
      "epsilon":    0.01,              # the threshold that was applied
      "passed":     true,             # error <= epsilon everywhere
      "max_error":  1.23e-04,         # worst error over all quantities
      "quantities": {                 # per-quantity achieved error
        "elastic_energy": 1.23e-04,
        "seismic_moment": 4.56e-06
      }
    }
"""

from __future__ import annotations

import json
import math
from typing import Mapping

SCHEMA_VERSION = 1


def _json_safe(x: float) -> float | str:
    """JSON has no inf/nan; encode them as strings so the file stays valid."""
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    if math.isnan(x):
        return "nan"
    return float(x)


def write_report_json(
    path: str,
    category: str,
    epsilon: float,
    passed: bool,
    quantities: Mapping[str, float],
) -> None:
    """Write the achieved-error summary for one comparison to ``path``."""
    numeric = [float(v) for v in quantities.values()]
    max_error = max(numeric) if numeric else 0.0
    payload = {
        "schema": SCHEMA_VERSION,
        "category": category,
        "epsilon": float(epsilon),
        "passed": bool(passed),
        "max_error": _json_safe(max_error),
        "quantities": {str(k): _json_safe(float(v)) for k, v in quantities.items()},
    }
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")


def read_report_json(path: str) -> dict:
    """Read a summary written by :func:`write_report_json`, decoding inf/nan."""
    with open(path) as fh:
        payload = json.load(fh)

    def _decode(v):
        if isinstance(v, str):
            return {"inf": math.inf, "-inf": -math.inf, "nan": math.nan}.get(v, v)
        return v

    payload["max_error"] = _decode(payload.get("max_error", 0.0))
    payload["quantities"] = {
        k: _decode(v) for k, v in payload.get("quantities", {}).items()
    }
    return payload
