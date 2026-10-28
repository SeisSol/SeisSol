#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Strong scaling series of SeisSol, from a campaign file.

    scaling.py plan    campaign.yaml   # the runs and jobs it would submit
    scaling.py check   campaign.yaml   # layouts, tools, SeisSol --check
    scaling.py submit  campaign.yaml   # build, check, submit what is open
    scaling.py status  campaign.yaml   # the state of every run
    scaling.py report  campaign.yaml   # tables and plots of what ran

Every command can be repeated: submit attempts again only the runs that did
not complete. See docs/strong-scaling.rst.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from strongscaling.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
