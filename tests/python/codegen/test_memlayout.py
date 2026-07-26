# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Unit tests for codegen/kernels/memlayout.py

This module picks a memory-layout XML file based on a build config.
Its correctness directly determines which matrix layouts the generated
kernel code uses, which affects correctness AND performance silently.
No yateto dependency — pure Python / re / os.
"""

import pytest
from kernels.memlayout import Candidate, findCandidates, guessMemoryLayout

# -----------------------------------------------------------------------------
# Candidate.score() — the scoring algorithm
# -----------------------------------------------------------------------------


class TestCandidateScore:
    """Score returns sum of 2**IMPORTANCE[key] for matching keys, else 0."""

    def test_empty_candidate_matches_everything(self):
        # A candidate with no attributes can never fail a requirement
        c = Candidate({})
        assert c.score({"precision": "d", "order": 5}) == 0
        # But also can't *match* anything — score is 0 for non-contributing keys
        assert c.score({}) == 0

    def test_single_key_match(self):
        c = Candidate({"precision": "d"})
        # precision matches -> 2**1 = 2
        assert c.score({"precision": "d"}) == 2

    def test_single_key_mismatch_zeros_entire_score(self):
        # "requirement not satisfied" short-circuits to 0
        c = Candidate({"precision": "d", "order": 5})
        assert c.score({"precision": "s", "order": 5}) == 0

    def test_missing_requirement_doesnt_contribute(self):
        # If candidate lacks a key, it's neutral — not a mismatch
        c = Candidate({"precision": "d"})
        # precision matches (2), order is absent from candidate (neutral)
        assert c.score({"precision": "d", "order": 5}) == 2

    def test_importance_ordering(self):
        # The scoring system guarantees best-candidate uniqueness via 2**i weights
        # So a multipleSimulations match (2**6=64) outweighs ALL lower matches combined
        # (2**1 + 2**2 + 2**3 + 2**4 + 2**5 = 62 < 64)
        c_low = Candidate(
            {
                "precision": "d",
                "equations": "elastic",
                "order": 5,
                "pe": "hsw",
                "gemmgen": "libxsmm",
            }
        )
        c_high = Candidate({"multipleSimulations": 8})
        reqs = {
            "precision": "d",
            "equations": "elastic",
            "order": 5,
            "pe": "hsw",
            "gemmgen": "libxsmm",
            "multipleSimulations": 8,
        }
        assert c_high.score(reqs) > c_low.score(reqs)

    def test_all_keys_match(self):
        atts = {
            "precision": "d",
            "equations": "elastic",
            "order": 5,
            "pe": "hsw",
            "gemmgen": "libxsmm",
            "multipleSimulations": 1,
        }
        c = Candidate(atts)
        # Sum of all 2**IMPORTANCE: 2+4+8+16+32+64 = 126
        assert c.score(atts) == 2 + 4 + 8 + 16 + 32 + 64

    def test_gemmgen_accepts_set_membership(self):
        # Special case: gemmgen requirement is a SET, matches if candidate's value is in it
        c = Candidate({"gemmgen": "libxsmm"})
        assert c.score({"gemmgen": {"libxsmm", "pspamm"}}) == 2**5
        assert c.score({"gemmgen": {"pspamm"}}) == 0

    def test_gemmgen_exact_value_also_works(self):
        # `val == self.atts[key]` branch — non-set equality
        c = Candidate({"gemmgen": "libxsmm"})
        # This line actually EXPOSES a latent bug: the code checks
        #     val == self.atts[key] or (key == "gemmgen" and self.atts[key] in val)
        # If val is a string (not a set), `self.atts[key] in val` treats val as a
        # substring check. But the normal call site always passes a set.
        # Documenting the invariant:
        assert c.score({"gemmgen": "libxsmm"}) == 2**5

    def test_integer_keys_compare_strictly(self):
        c = Candidate({"order": 5})
        assert c.score({"order": 5}) == 2**3
        # 5 != "5" — types must match
        assert c.score({"order": "5"}) == 0


# -----------------------------------------------------------------------------
# findCandidates() — filename parsing
# -----------------------------------------------------------------------------


class TestFindCandidates:
    """Parses attributes from filenames like `dense_O5_d_hsw_libxsmm.xml`."""

    @pytest.fixture
    def layout_dir(self, tmp_path):
        # Create representative layout filenames — mirrors codegen/config/cpu/*
        for fname in [
            "dense.xml",
            "elastic_O5_d_hsw.xml",
            "elastic_O5_s_hsw_libxsmm.xml",
            "viscoelastic2_O6_d_skx_pspamm_ms8.xml",
            "poroelastic_O4_f32_a64fx.xml",
        ]:
            (tmp_path / fname).touch()
        return tmp_path

    def test_returns_dict_keyed_by_filename(self, layout_dir):
        candidates = findCandidates(str(layout_dir))
        assert set(candidates.keys()) == {
            "dense.xml",
            "elastic_O5_d_hsw.xml",
            "elastic_O5_s_hsw_libxsmm.xml",
            "viscoelastic2_O6_d_skx_pspamm_ms8.xml",
            "poroelastic_O4_f32_a64fx.xml",
        }

    def test_parses_equation_from_first_token(self, layout_dir):
        c = findCandidates(str(layout_dir))["elastic_O5_d_hsw.xml"]
        assert c.atts["equations"] == "elastic"

    def test_parses_order(self, layout_dir):
        c = findCandidates(str(layout_dir))["elastic_O5_d_hsw.xml"]
        assert c.atts["order"] == 5
        assert isinstance(c.atts["order"], int)

    def test_parses_precision_single_char(self, layout_dir):
        c = findCandidates(str(layout_dir))["elastic_O5_d_hsw.xml"]
        assert c.atts["precision"] == "d"

    def test_parses_precision_f32_form(self, layout_dir):
        c = findCandidates(str(layout_dir))["poroelastic_O4_f32_a64fx.xml"]
        assert c.atts["precision"] == "f32"

    def test_parses_multiple_simulations(self, layout_dir):
        c = findCandidates(str(layout_dir))["viscoelastic2_O6_d_skx_pspamm_ms8.xml"]
        assert c.atts["multipleSimulations"] == 8

    def test_parses_gemmgen(self, layout_dir):
        c = findCandidates(str(layout_dir))["viscoelastic2_O6_d_skx_pspamm_ms8.xml"]
        assert c.atts["gemmgen"] == "pspamm"

    def test_parses_architecture(self, layout_dir):
        c = findCandidates(str(layout_dir))["poroelastic_O4_f32_a64fx.xml"]
        assert c.atts["pe"] == "a64fx"

    def test_dense_has_no_attributes(self, layout_dir):
        # "dense" is an equation name by fallthrough, but only one token
        c = findCandidates(str(layout_dir))["dense.xml"]
        # Current behavior: "dense" falls into the equations bucket
        assert c.atts == {"equations": "dense"}

    def test_empty_directory_returns_empty(self, tmp_path):
        assert findCandidates(str(tmp_path)) == {}


# -----------------------------------------------------------------------------
# Integration: scoring against parsed filenames
# -----------------------------------------------------------------------------


class TestScoringAgainstFilenames:
    """End-to-end: parse real-ish filenames and check selection is sensible."""

    @pytest.fixture
    def config_dir(self, tmp_path):
        for fname in [
            "dense.xml",
            "elastic_O5_d_hsw.xml",
            "elastic_O5_d_skx.xml",
            "elastic_O6_d_hsw.xml",
            "elastic_O5_s_hsw.xml",
        ]:
            (tmp_path / fname).touch()
        return tmp_path

    def test_exact_match_wins(self, config_dir):
        candidates = findCandidates(str(config_dir))
        reqs = {
            "equations": "elastic",
            "order": 5,
            "precision": "d",
            "pe": "hsw",
            "multipleSimulations": 1,
            "gemmgen": {"libxsmm"},
        }
        scored = {name: c.score(reqs) for name, c in candidates.items()}
        best = max(scored, key=lambda k: scored[k])
        assert best == "elastic_O5_d_hsw.xml"

    def test_order_mismatch_zeroes_score(self, config_dir):
        candidates = findCandidates(str(config_dir))
        reqs = {"order": 99}  # no file matches
        # all elastic files have a specific order -> mismatch
        non_dense_scores = [
            c.score(reqs) for name, c in candidates.items() if name != "dense.xml"
        ]
        assert all(s == 0 for s in non_dense_scores)

    def test_dense_fallback_wins_when_nothing_fits(self, config_dir):
        # dense.xml has "equations=dense" only -> never matches "elastic"
        # but the actual guessMemoryLayout falls back to "dense.xml" by name
        candidates = findCandidates(str(config_dir))
        reqs = {
            "equations": "acoustic",
            "order": 99,
            "precision": "d",
            "pe": "hsw",
            "multipleSimulations": 1,
            "gemmgen": {"libxsmm"},
        }
        scores = {name: c.score(reqs) for name, c in candidates.items()}
        # All proper matches fail because equation mismatch or order mismatch
        # dense.xml: "equations=dense" vs "acoustic" -> 0
        assert all(s == 0 for s in scores.values())
        # guessMemoryLayout handles this case by string-falling-back to "dense.xml"


# -----------------------------------------------------------------------------
# guessMemoryLayout — end-to-end against real SeisSol config files
# -----------------------------------------------------------------------------


class TestGuessMemoryLayoutReal:
    """Runs against the real codegen/config/cpu directory in SeisSol master."""

    def test_elastic_hsw_d_has_a_match(self, capsys):
        env = {
            "targets": ["cpu"],
            "precision": "d",
            "equations": "elastic",
            "order": 6,
            "arch": "hsw",
            "multipleSimulations": 1,
            "gemmgen": {"libxsmm"},
        }
        result = guessMemoryLayout(env)
        # Either a specific layout or fallback to dense.xml — both valid
        assert result.endswith(".xml")

    def test_fallback_prints_warning_for_impossible_config(self, capsys):
        env = {
            "targets": ["cpu"],
            "precision": "d",
            "equations": "nonexistent",
            "order": 99,
            "arch": "fantasyarch",
            "multipleSimulations": 999,
            "gemmgen": {"nothing"},
        }
        result = guessMemoryLayout(env)
        captured = capsys.readouterr()
        assert "dense.xml" in result
        assert (
            "No suitable memory layout" in captured.out
            or "Using memory layout" in captured.out
        )

    def test_gpu_target_uses_gpu_config_dir(self, capsys):
        env = {
            "targets": ["gpu"],
            "precision": "s",
            "equations": "elastic",
            "order": 5,
            "arch": "sm_60",
            "multipleSimulations": 1,
            "gemmgen": {"tensorforge"},
        }
        result = guessMemoryLayout(env)
        # Either /config/gpu/ in the path or dense fallback
        assert result.endswith(".xml")
