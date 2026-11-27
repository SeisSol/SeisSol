#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The include graph of the C++ sources of SeisSol, as layering.py and restructure.py read it."""

import os
import re
import subprocess

INCLUDE_RE = re.compile(r'^(\s*#\s*include\s*)(["<])([^">]+)([">])')
CPP_SUFFIXES = (".h", ".hpp", ".cuh", ".cpp", ".cu", ".cc", ".h.in", ".cpp.in")


def repo_root():
    result = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def tracked_files(root, *dirs):
    """The files git tracks below `dirs`, relative to the repository."""
    result = subprocess.run(
        ["git", "-C", root, "ls-files", "-z", "--", *dirs],
        check=True,
        capture_output=True,
        text=True,
    )
    return sorted(path for path in result.stdout.split("\0") if path)


def is_cpp(path):
    return path.endswith(CPP_SUFFIXES)


class Sources:
    """The tracked C++ files of src/, app/ and tests/."""

    def __init__(self, root):
        self.root = root
        self.files = [
            path for path in tracked_files(root, "src", "app", "tests") if is_cpp(path)
        ]
        self.existing = set(self.files)
        # relative to src/
        self.src = {path[4:] for path in self.files if path.startswith("src/")}

    def resolve(self, includer, spelled):
        """The file below src/ that `#include "spelled"` in `includer` (relative to the
        repository) names, relative to src/, or None; and whether the include names it relative
        to the directory of the includer.

        The compiler looks in the directory of the includer first, then in tests/ for the tests,
        then in src/."""
        if includer.startswith("src/"):
            here = os.path.dirname(includer[4:])
            candidate = os.path.normpath(os.path.join(here, spelled))
            if candidate in self.src:
                return candidate, candidate != os.path.normpath(spelled)
        else:
            candidate = os.path.normpath(
                os.path.join(os.path.dirname(includer), spelled)
            )
            if candidate in self.existing:
                return None, False
            if includer.startswith("tests/"):
                if os.path.normpath(os.path.join("tests", spelled)) in self.existing:
                    return None, False
        candidate = os.path.normpath(spelled)
        if candidate in self.src:
            return candidate, False
        return None, False


def read_includes(root, path):
    """[(line index, bracket, spelled name)] of the includes of a file (relative to the
    repository)."""
    result = []
    with open(os.path.join(root, path), encoding="utf-8", errors="replace") as file:
        for index, line in enumerate(file):
            match = INCLUDE_RE.match(line)
            if match:
                result.append((index, match.group(2), match.group(3)))
    return result
