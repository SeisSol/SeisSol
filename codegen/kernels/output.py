# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Writing the files of the code generator."""

import os


def write_if_changed(path, content):
    """Writes `content` into the file at `path`, unless the file holds it already.

    The build compares the time of each file of the code generator with the
    times of what it built from it. CMake runs the code generator at every
    configure, to collect the files it generates; that run writes some headers
    and lists, and rewriting them with the same content would rebuild all that
    depends on them.
    """
    if os.path.exists(path):
        with open(path) as file:
            if file.read() == content:
                return
    with open(path, "w") as file:
        file.write(content)
