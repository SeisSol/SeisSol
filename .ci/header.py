# SPDX-FileCopyrightText: 2024 Technical University of Munich
#
# SPDX-License-Identifier: BSD-3-Clause

import argparse
import datetime
import os
import pathlib
import re
import stat
import subprocess
import sys
from typing import List


class Settings:
    MAIN_DEVELOPER = "SeisSol Group"
    PROJECT_NAME = "SeisSol"
    PROJECT_WEBPAGE = "www.seissol.org"
    LICENSE = "BSD-3-Clause"

    CPP_HEADER_ENDINGS = [".h", ".hpp", ".cuh", ".hh"]
    CPP_SOURCE_ENDINGS = [".cpp", ".cu", ".cc"]
    PY_SOURCE_ENDINGS = [".py"]
    SHELL_SOURCE_ENDINGS = [".sh"]
    CMAKE_SOURCE_ENDINGS = [".cmake"]
    CMAKE_FILE_NAMES = ["CMakeLists.txt"]


class Span:
    def __init__(self, start, end=None):
        self.start = start
        if end is None:
            self.end = start
        else:
            self.end = end

    def unite(self, other):
        return Span(min(self.start, other.start), max(self.end, other.end))

    def extend(self, year):
        return Span(min(self.start, year), max(self.end, year))

    def __str__(self):
        if self.start == self.end:
            return f"{self.start}"
        else:
            return f"{self.start}-{self.end}"


def rootpath():
    path = subprocess.getoutput("git rev-parse --show-toplevel")
    return pathlib.Path(path)


def textCopyrights(copyrights, authorspans):
    crregex = re.compile(
        r"(?:(?:[Cc][Oo][Pp][Yy][Rr][Ii][Gg][Hh][Tt]\s+)|(?:SPDX\-FileCopyrightText\:\s*)|[\(\[]\s*[Cc]\s*[\)\]])+\s*(?:(\d+)(?:\s*-\s*(\d+))?)?[\s,;]*([\w\d\s\-_]*)"
    )
    for cr in copyrights:
        result = crregex.search(cr)
        holder = result.group(3).strip()
        year = result.group(1)
        year2 = result.group(2)
        if year2 is None:
            timespan = Span(int(year))
        else:
            timespan = Span(int(year), int(year2))

        if holder not in authorspans:
            authorspans[holder] = timespan
        else:
            authorspans[holder] = authorspans[holder].unite(timespan)


def textAuthors(authors):
    aregex = re.compile(
        r"(?:(?:[Aa][Uu][Tt][Hh][Oo][Rr]\s+)|(?:SPDX\-FileContributor\:\s*))(.*)"
    )
    alines = []
    for author in authors:
        result = aregex.search(author)
        holder = result.group(1).strip()
        if holder != "Author lists in /AUTHORS and /CITATION.cff":
            alines += [holder]
    return alines


def licenseHeader(authorspans, authors, commentstyle="//"):
    # the main developer always comes first
    if Settings.MAIN_DEVELOPER is not None:
        authorlist = [Settings.MAIN_DEVELOPER] + [
            author for author in authorspans if author != Settings.MAIN_DEVELOPER
        ]
    else:
        authorlist = [author for author in authorspans]

    # no titleline for now
    # titleline = f'{commentstyle} This file is part of {Settings.PROJECT_NAME} ( {Settings.PROJECT_WEBPAGE} ).'

    # always add:
    # SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
    #
    # SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

    crlines = []
    for author in authorlist:
        # TODO: unite year spans
        timespan = authorspans[author]
        crlines += [f"{commentstyle} SPDX-FileCopyrightText: {timespan} {author}"]
    alines = []
    for author in authors:
        alines += [f"{commentstyle} SPDX-FileContributor: {author}"]
    return (
        crlines
        + [
            f"{commentstyle}",
            f"{commentstyle} SPDX-License-Identifier: {Settings.LICENSE}",
            f"{commentstyle} SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/",
            f"{commentstyle}",
            f"{commentstyle} SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff",
        ]
        + alines
    )


def makeLicense(path, copyrights, authors, commentstyle="//"):
    authorspans = {}
    textCopyrights(copyrights, authorspans)
    authors = textAuthors(authors)

    # let the main developer always hold the copyright for the current year
    if Settings.MAIN_DEVELOPER is not None:
        currentYear = datetime.date.today().year
        if Settings.MAIN_DEVELOPER not in authorspans:
            authorspans[Settings.MAIN_DEVELOPER] = Span(currentYear)
        else:
            authorspans[Settings.MAIN_DEVELOPER] = authorspans[
                Settings.MAIN_DEVELOPER
            ].extend(currentYear)

    # ...or just set it to the first year
    if Settings.MAIN_DEVELOPER is not None:
        authorspans[Settings.MAIN_DEVELOPER] = Span(
            authorspans[Settings.MAIN_DEVELOPER].start
        )

    return licenseHeader(authorspans, authors, commentstyle)


def licenseCleaner(path, lines, commentstyle="//"):
    inComment = False
    copyrights = []
    authors = []

    rest = []

    for i, line in enumerate(lines):
        if (line.strip().startswith("#") and commentstyle != "#") or (
            not inComment
            and not line.strip().startswith(commentstyle)
            and not line.strip().startswith("/*")
        ):
            rest = lines[i:]
            break

        commentStart = line.find("/*")
        if commentStart != -1:
            inComment = True

        commentEnd = line.find("*/")
        if inComment or line.strip().startswith(commentstyle):
            copyrightpos = line.lower().find("copyright")
            if (
                copyrightpos != -1
                and (commentEnd == -1 or copyrightpos < commentEnd)
                and (commentStart == -1 or copyrightpos > commentStart)
            ):
                copyrights += [line]
            authorpos = line.lower().find("author")
            if authorpos == -1:
                authorpos = line.lower().find("contributor")
            if (
                authorpos != -1
                and (commentEnd == -1 or authorpos < commentEnd)
                and (commentStart == -1 or authorpos > commentStart)
            ):
                authors += [line]

        if commentEnd != -1:
            inComment = False

    return makeLicense(path, copyrights, authors, commentstyle) + rest


def cppCanonicalHeaders(path: pathlib.Path, lines: List[str]):
    # remove pragma once
    for i, line in enumerate(lines):
        if line.strip() == "#pragma once":
            preheader = lines[:i]
            innerlines = lines[i + 1 :]
            break

    # TODO: check if there's an include guard at all
    # remove default include guard
    ifndefline = None
    innerlinesPre = None
    for i, line in enumerate(lines):
        if line.strip().startswith("#ifndef"):
            ifndefline = i
        elif line.strip().startswith("#define") and ifndefline is not None:
            preheader = lines[:ifndefline]
            innerlinesPre = lines[i + 1 :]
            break
        elif not (line.strip().startswith("//") or line.strip() == ""):
            preheader = lines[:i]
            innerlinesPre = lines[i + 1 :]
            break
    innerlines = None
    for i, line in reversed(list(enumerate(innerlinesPre))):
        if line.strip().startswith("#endif"):
            innerlines = innerlinesPre[:i]
            break
        elif not (line.strip().startswith("//") or line.strip() == ""):
            innerlines = innerlinesPre[: (i + 1)]
            break

    # re-insert new include guard
    strpath = str(path.relative_to(rootpath()))
    sanitizedpath = strpath.replace("/", "_").replace(".", "_").upper()
    projectprefix = Settings.PROJECT_NAME.upper()
    guard = f"{projectprefix}_{sanitizedpath}_"

    return (
        preheader
        + [f"#ifndef {guard}", f"#define {guard}"]
        + innerlines
        + [f"#endif // {guard}"]
    )


def addShebang(path, lines, program):
    # only add shebang for files that are executable

    stats = path.stat().st_mode
    if (
        (stats & stat.S_IXUSR) != 0
        or (stats & stat.S_IXGRP) != 0
        or (stats & stat.S_IXOTH) != 0
    ):
        shebang = f"#!/usr/bin/env {program}"
        if lines[0].strip().startswith("#!"):
            return [shebang, ""] + lines[1:]
        else:
            return [shebang, ""] + lines
    else:
        return lines


def sanitizeLineEndings(lines):
    return [line if line.endswith("\n") else line + "\n" for line in lines]


def processFile(path, dryrun):
    pathobj = pathlib.Path(path).resolve()
    with open(path) as file:
        lines = file.readlines()

    if pathobj.suffix in Settings.PY_SOURCE_ENDINGS:
        lines = licenseCleaner(pathobj, lines, "#")
        lines = addShebang(pathobj, lines, "python3")
    if pathobj.suffix in Settings.SHELL_SOURCE_ENDINGS:
        lines = licenseCleaner(pathobj, lines, "#")
        lines = addShebang(pathobj, lines, "sh")
    if (
        pathobj.suffix in Settings.CMAKE_SOURCE_ENDINGS
        or pathobj.name in Settings.CMAKE_FILE_NAMES
    ):
        lines = licenseCleaner(pathobj, lines, "#")
    if pathobj.suffix in Settings.CPP_HEADER_ENDINGS + Settings.CPP_SOURCE_ENDINGS:
        lines = licenseCleaner(pathobj, lines, "//")
    if pathobj.suffix in Settings.CPP_HEADER_ENDINGS and not pathobj.name.endswith(
        ".t.h"
    ):
        lines = cppCanonicalHeaders(pathobj, lines)
    lines = sanitizeLineEndings(lines)

    with open(path) as file:
        linescomp = file.readlines()
    if not dryrun:
        with open(path, "w") as file:
            file.writelines(lines)
    return lines != linescomp


def main():
    print("SeisSol File Sanitizer")

    argParser = argparse.ArgumentParser()
    argParser.add_argument("--fix", action="store_true")
    argParser.add_argument("path")
    argParser.set_defaults(fix=False)
    args = argParser.parse_args()

    found = False

    for root, dirs, files in pathlib.Path(args.path).walk(top_down=False):
        for prefile in files:
            file = root / prefile
            if (
                not os.path.islink(file)
                and file.suffix
                in Settings.CPP_HEADER_ENDINGS + Settings.CPP_SOURCE_ENDINGS
            ):
                result = processFile(file, not args.fix)
                if result:
                    print(f"Reformatted: {file}")
                found |= result

    if found:
        print("Unsanitized files found.")
        sys.exit(1)
    else:
        print("All files conformant.")


if __name__ == "__main__":
    main()
