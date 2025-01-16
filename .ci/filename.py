# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import argparse
import pathlib
import sys

# a script to check if headers use the .h ending
# and if their names are PascalCase.


def main():
    print("SeisSol Filename Conformity Checker")

    argParser = argparse.ArgumentParser()
    argParser.add_argument("--fix", action="store_true")
    argParser.add_argument("--dirs", action="store_true")
    argParser.add_argument("path")
    argParser.set_defaults(fix=False, dirs=False)
    args = argParser.parse_args()

    class Found:
        found = False

    foundobj = Found()

    def fixName(name):
        parts = name.split("_")
        return "".join(
            part[0].upper() + part[1:] for part in parts if len(part) > 0
        )

    def sanitizeFile(file, name):
        if file.name != name:
            print(f"{file} should be {name}")
            foundobj.found = True
            if args.fix:
                file.rename(file.parent / name)

    def checkFile(file, suffixmask, suffix):
        if suffixmask is None or file.suffix == suffixmask:
            sanitized = fixName(file.stem) + suffix
            sanitizeFile(file, sanitized)

    for root, dirs, files in pathlib.Path(args.path).walk(top_down=False):
        for prefile in files:
            checkFile(root / prefile, ".hpp", ".h")
            checkFile(root / prefile, ".h", ".h")
            checkFile(root / prefile, ".cpp", ".cpp")
        if args.dirs:
            for prefile in dirs:
                checkFile(root / prefile, None, "")

    if foundobj.found:
        print("Unconformant files found.")
        sys.exit(1)
    else:
        print("All files conformant.")


if __name__ == "__main__":
    main()
