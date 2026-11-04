# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The inputs of a step of the code generator, so that the build runs it again
whenever one of them changes: the modules it has imported, including those of
yateto and of the code generators behind it, and the data files it has read,
such as the matrices. A step writes them into a depfile in the format of make,
which CMake hands to the build tool (DEPFILE of add_custom_command)."""

import os
import sys
import sysconfig


class InputRecorder:
    """Records the files under the directories `roots` that the process opens
    for reading, from its construction on."""

    def __init__(self, roots):
        self._roots = [os.path.join(os.path.abspath(root), "") for root in roots]
        self.files = set()
        sys.addaudithook(self._audit)

    def _audit(self, event, args):
        if event != "open" or len(args) < 2:
            return
        path, mode = args[0], args[1]
        if not isinstance(path, (str, bytes, os.PathLike)) or not isinstance(mode, str):
            return
        if any(flag in mode for flag in "wax+"):
            return
        path = os.path.abspath(os.fsdecode(path))
        # the compiled modules that the import reads are written by the imports
        # as well; the modules themselves are what imported_files gives
        if path.endswith(".pyc"):
            return
        if any(path.startswith(root) for root in self._roots):
            self.files.add(path)


def imported_files():
    """The files of the modules that the process has imported, but those of the
    standard library (the installed packages may be in a directory of it)."""
    paths = sysconfig.get_paths()

    def directories(*keys):
        return [os.path.join(os.path.abspath(paths[key]), "") for key in keys]

    standard = directories("stdlib", "platstdlib")
    installed = directories("purelib", "platlib")
    files = set()
    for module in list(sys.modules.values()):
        file = getattr(module, "__file__", None)
        if file is None:
            continue
        file = os.path.abspath(file)
        if not os.path.isfile(file):
            continue
        if any(file.startswith(prefix) for prefix in installed) or not any(
            file.startswith(prefix) for prefix in standard
        ):
            files.add(file)
    return files


def _escape(path):
    return (
        path.replace("\\", "\\\\")
        .replace(" ", "\\ ")
        .replace("#", "\\#")
        .replace("$", "$$")
    )


def write_depfile(path, output, inputs):
    """Writes the depfile `path`: the file `output`, the first output of the
    step, by which the build tool knows the step, depends on the files `inputs`.
    It names no other output: Ninja stops at a depfile that names an output the
    step no longer has, as once a configuration has left the build."""
    dependencies = "".join(f" \\\n  {_escape(file)}" for file in sorted(inputs))
    temporary = f"{path}.tmp"
    with open(temporary, "w") as file:
        file.write(f"{_escape(os.path.abspath(output))}:{dependencies}\n")
    os.replace(temporary, path)
