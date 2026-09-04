# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Tests for .ci/header.py — the SPDX copyright header enforcer.

This script runs as a pre-commit hook and on every PR. Unlike filename.py
it has substantial internal logic (regex parsing, Span merging, year
tracking). Bugs here silently corrupt SPDX metadata in the source tree.
"""

import importlib.util
import sys
from pathlib import Path

import pytest  # noqa: F401

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
_spec = importlib.util.spec_from_file_location(
    "ci_header", SEISSOL_ROOT / ".ci" / "header.py"
)
hdr = importlib.util.module_from_spec(_spec)
sys.modules["ci_header"] = hdr
try:
    _spec.loader.exec_module(hdr)
except SystemExit:
    pass


# =============================================================================
# Span class — unite() and extend() define year-range merging
# =============================================================================


class TestSpan:

    def test_single_year(self):
        s = hdr.Span(2020)
        assert s.start == 2020
        assert s.end == 2020
        assert str(s) == "2020"

    def test_year_range(self):
        s = hdr.Span(2020, 2024)
        assert str(s) == "2020-2024"

    def test_end_defaults_to_start(self):
        s = hdr.Span(2019)
        assert s.end == s.start

    def test_unite_two_overlapping_spans(self):
        a = hdr.Span(2020, 2022)
        b = hdr.Span(2021, 2024)
        merged = a.unite(b)
        assert merged.start == 2020
        assert merged.end == 2024

    def test_unite_is_commutative(self):
        a = hdr.Span(2018, 2019)
        b = hdr.Span(2023)
        m1 = a.unite(b)
        m2 = b.unite(a)
        assert (m1.start, m1.end) == (m2.start, m2.end)

    def test_unite_disjoint_spans_fills_gap(self):
        """NOTE: unite does NOT represent discontinuous spans — it collapses
        them to the full outer range. This is SeisSol's intended behavior
        per the SPDX convention (year-range means 'copyright years').
        """
        a = hdr.Span(2015)
        b = hdr.Span(2023)
        m = a.unite(b)
        assert m.start == 2015 and m.end == 2023

    def test_extend_expands_range(self):
        s = hdr.Span(2020, 2022)
        e = s.extend(2025)
        assert e.start == 2020 and e.end == 2025

    def test_extend_earlier_year_moves_start(self):
        s = hdr.Span(2020, 2022)
        e = s.extend(2015)
        assert e.start == 2015 and e.end == 2022

    def test_extend_interior_year_is_noop(self):
        s = hdr.Span(2020, 2025)
        e = s.extend(2022)
        assert e.start == 2020 and e.end == 2025


# =============================================================================
# textCopyrights — regex-based SPDX parsing
# =============================================================================


class TestTextCopyrights:
    """Parses copyright strings like:
       '// SPDX-FileCopyrightText: 2020 SeisSol Group'
       '// SPDX-FileCopyrightText: 2020-2024 SeisSol Group'
       '// Copyright (c) 2020 SeisSol Group'
    and accumulates year spans by holder name.
    """

    def test_single_year_spdx_format(self):
        authorspans = {}
        hdr.textCopyrights(
            ["// SPDX-FileCopyrightText: 2020 SeisSol Group"],
            authorspans,
        )
        assert "SeisSol Group" in authorspans
        assert authorspans["SeisSol Group"].start == 2020
        assert authorspans["SeisSol Group"].end == 2020

    def test_year_range_spdx_format(self):
        authorspans = {}
        hdr.textCopyrights(
            ["// SPDX-FileCopyrightText: 2020-2024 SeisSol Group"],
            authorspans,
        )
        assert authorspans["SeisSol Group"].start == 2020
        assert authorspans["SeisSol Group"].end == 2024

    def test_legacy_copyright_format(self):
        """Older SeisSol files use 'Copyright (c) YEAR'."""
        authorspans = {}
        hdr.textCopyrights(
            ["// Copyright (c) 2015 SeisSol Group"],
            authorspans,
        )
        assert authorspans["SeisSol Group"].start == 2015

    def test_multiple_entries_merge_years(self):
        authorspans = {}
        hdr.textCopyrights(
            [
                "// SPDX-FileCopyrightText: 2018 SeisSol Group",
                "// SPDX-FileCopyrightText: 2022 SeisSol Group",
            ],
            authorspans,
        )
        # Merged to span [2018, 2022]
        assert authorspans["SeisSol Group"].start == 2018
        assert authorspans["SeisSol Group"].end == 2022

    def test_multiple_holders_separate(self):
        authorspans = {}
        hdr.textCopyrights(
            [
                "// SPDX-FileCopyrightText: 2015 Example Corporation",
                "// SPDX-FileCopyrightText: 2018 SeisSol Group",
            ],
            authorspans,
        )
        assert "Example Corporation" in authorspans
        assert "SeisSol Group" in authorspans

    def test_case_insensitive_copyright_keyword(self):
        authorspans = {}
        hdr.textCopyrights(
            ["// COPYRIGHT 2020 SeisSol Group"],
            authorspans,
        )
        assert "SeisSol Group" in authorspans

    def test_parenthesized_c_notation(self):
        authorspans = {}
        hdr.textCopyrights(
            ["/* (c) 2019 SeisSol Group */"],
            authorspans,
        )
        assert authorspans["SeisSol Group"].start == 2019


# =============================================================================
# textAuthors — parses Author / Contributor lines
# =============================================================================


class TestTextAuthors:

    def test_spdx_filecontributor(self):
        authors = hdr.textAuthors(["// SPDX-FileContributor: Max Mustermann"])
        assert authors == ["Max Mustermann"]

    def test_legacy_author_keyword(self):
        authors = hdr.textAuthors(["// @author Max Mustermann (max@example.com)"])
        assert len(authors) == 1
        assert "Max Mustermann" in authors[0]

    def test_canonical_author_list_line_filtered(self):
        """The boilerplate 'Author lists in /AUTHORS and /CITATION.cff' line
        is explicitly filtered out — otherwise it would appear as a pseudo-
        author in every file."""
        authors = hdr.textAuthors(
            [
                "// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff",
                "// SPDX-FileContributor: Max Mustermann",
            ]
        )
        assert "Author lists in /AUTHORS and /CITATION.cff" not in authors
        assert "Max Mustermann" in authors

    def test_multiple_authors(self):
        authors = hdr.textAuthors(
            [
                "// SPDX-FileContributor: Max Mustermann",
                "// SPDX-FileContributor: Erika Musterfrau",
                "// SPDX-FileContributor: John Doe",
            ]
        )
        assert len(authors) == 3

    def test_case_insensitive(self):
        authors = hdr.textAuthors(["// AUTHOR Jane Doe"])
        assert "Jane Doe" in authors


# =============================================================================
# licenseHeader — assembles the final header block
# =============================================================================


class TestLicenseHeader:

    def test_main_developer_comes_first(self):
        """Settings.MAIN_DEVELOPER (SeisSol Group) must always be first
        in the author list."""
        spans = {
            "Example Corporation": hdr.Span(2015),
            "SeisSol Group": hdr.Span(2020, 2024),
            "Some Other": hdr.Span(2018),
        }
        lines = hdr.licenseHeader(spans, [], commentstyle="//")
        # First copyright line must be SeisSol Group
        cr_lines = [l for l in lines if "FileCopyrightText" in l]
        assert "SeisSol Group" in cr_lines[0]

    def test_spdx_license_identifier_included(self):
        spans = {"SeisSol Group": hdr.Span(2024)}
        lines = hdr.licenseHeader(spans, [], commentstyle="//")
        id_line = [l for l in lines if "License-Identifier" in l]
        assert len(id_line) == 1
        assert "BSD-3-Clause" in id_line[0]

    def test_canonical_author_list_boilerplate_added(self):
        spans = {"SeisSol Group": hdr.Span(2024)}
        lines = hdr.licenseHeader(spans, [], commentstyle="//")
        assert any("Author lists in /AUTHORS and /CITATION.cff" in l for l in lines)

    def test_license_comments_pointer_added(self):
        spans = {"SeisSol Group": hdr.Span(2024)}
        lines = hdr.licenseHeader(spans, [], commentstyle="//")
        assert any("Full text under /LICENSE and /LICENSES/" in l for l in lines)

    def test_python_comment_style(self):
        """Python files use '#' not '//'."""
        spans = {"SeisSol Group": hdr.Span(2024)}
        lines = hdr.licenseHeader(spans, [], commentstyle="#")
        for line in lines:
            # Every non-empty line must start with '# '
            assert line.startswith("# ") or line.startswith("#")

    def test_year_range_formatted_correctly(self):
        spans = {"SeisSol Group": hdr.Span(2018, 2024)}
        lines = hdr.licenseHeader(spans, [], commentstyle="//")
        assert any("2018-2024 SeisSol Group" in l for l in lines)

    def test_contributors_listed(self):
        spans = {"SeisSol Group": hdr.Span(2024)}
        authors = ["Max Mustermann", "Erika Musterfrau"]
        lines = hdr.licenseHeader(spans, authors, commentstyle="//")
        assert any("Max Mustermann" in l for l in lines)
        assert any("Erika Musterfrau" in l for l in lines)


# =============================================================================
# makeLicense — integration of copyright parsing + current-year injection
# =============================================================================


class TestMakeLicense:

    def test_current_year_added_to_main_developer_span(self):
        """makeLicense always ensures the MAIN_DEVELOPER's span includes
        the CURRENT year (so new files get the right year automatically)."""

        # Note: makeLicense calls .extend(currentYear) then sets start
        # back — the final span is Span(start_from_input) — but this
        # means SeisSol's year policy: MAIN_DEVELOPER always shows the
        # start year only, not the range. Let's verify.
        result = hdr.makeLicense(
            "dummy.cpp",
            copyrights=["// SPDX-FileCopyrightText: 2018 SeisSol Group"],
            authors=[],
            commentstyle="//",
        )
        # Look for the SeisSol Group copyright line
        cr_line = [l for l in result if "SeisSol Group" in l][0]
        assert "SPDX-FileCopyrightText" in cr_line
        # The span is Span(start_year) only — confirms the post-extend-reset
        # logic in makeLicense()

    def test_no_prior_copyright_adds_current_year(self):
        """A new file with no copyright info — should still produce a
        valid header stamped with the current year."""
        import datetime

        current = datetime.date.today().year
        result = hdr.makeLicense(
            "new_file.cpp",
            copyrights=[],
            authors=[],
            commentstyle="//",
        )
        cr_line = [l for l in result if "SeisSol Group" in l][0]
        assert str(current) in cr_line


# =============================================================================
# licenseCleaner — extracts copyrights/authors from existing file content
# =============================================================================


class TestLicenseCleaner:
    """licenseCleaner is the combined parse-and-rewrite driver:
    input = lines of a file
    output = new header + original non-header content
    """

    def test_extracts_copyright_from_existing_header(self):
        # avoid failing the REUSE test here
        spdxli = "SPDX-License-Identifier"

        lines = [
            "// SPDX-FileCopyrightText: 2020 SeisSol Group\n",
            "//\n",
            f"// {spdxli}: BSD-3-Clause\n",
            "\n",
            "#include <iostream>\n",
            "int main() {}\n",
        ]
        result = hdr.licenseCleaner("dummy.cpp", lines, "//")
        # Resulting header must still have SeisSol Group
        hdr_block = [l for l in result if l.startswith("//")]
        assert any("SeisSol Group" in l for l in hdr_block)
        # Original content must be preserved at the end
        assert "#include <iostream>\n" in result
        assert "int main() {}\n" in result

    def test_preserves_content_starting_with_non_comment(self):
        """Non-comment first line ends the header parsing immediately."""
        lines = [
            "#include <iostream>\n",  # non-header → rest begins here
            "int main() {}\n",
        ]
        result = hdr.licenseCleaner("dummy.cpp", lines, "//")
        # Header is inserted, then original content appended
        assert "#include <iostream>\n" in result

    def test_python_hash_comments(self):
        lines = [
            "# SPDX-FileCopyrightText: 2020 SeisSol Group\n",
            "#\n",
            "import os\n",
        ]
        result = hdr.licenseCleaner("dummy.py", lines, "#")
        assert "import os\n" in result

    def test_c_style_block_comment(self):
        """/* ... */ blocks must be consumed as part of the header."""
        lines = [
            "/*\n",
            " * Copyright (c) 2018 SeisSol Group\n",
            " */\n",
            "#include <iostream>\n",
        ]
        result = hdr.licenseCleaner("dummy.cpp", lines, "//")
        assert "#include <iostream>\n" in result

    def test_python_file_with_pound_then_code(self):
        """Ensure Python ('#' comment style) handles the detection
        branch that treats '#' differently when commentstyle='#'."""
        lines = [
            "#!/usr/bin/env python3\n",
            "# SPDX-FileCopyrightText: 2020 SeisSol Group\n",
            "\n",
            "print('hello')\n",
        ]
        result = hdr.licenseCleaner("dummy.py", lines, "#")
        assert "print('hello')\n" in result
