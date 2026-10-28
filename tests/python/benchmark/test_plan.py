# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The planner picks the node counts and orders the runs of a job."""

import pytest
from strongscaling import config, plan


def test_a_job_per_node_count(campaign_file):
    campaign = config.load_campaign(campaign_file())
    result = plan.make_plan(campaign)
    assert [job.name for job in result.jobs] == ["n0001", "n0002", "n0004"]
    assert len(result.runs) == 3 * 2 * 2
    assert {run.job for run in result.runs} == {"n0001", "n0002", "n0004"}


def test_repetitions_take_turns(campaign_file):
    campaign = config.load_campaign(campaign_file())
    job = plan.make_plan(campaign).jobs[0]
    # the first repetition of every variant, then the second, in a rotated order
    assert job.runs == [
        "n0001-fake-a-r1",
        "n0001-fake-b-r1",
        "n0001-fake-b-r2",
        "n0001-fake-a-r2",
    ]


def test_blocked_order(campaign_file):
    campaign = config.load_campaign(
        campaign_file(scaling={"nodes": [1], "repetitions": 2, "order": "blocked"})
    )
    assert plan.make_plan(campaign).jobs[0].runs == [
        "n0001-fake-a-r1",
        "n0001-fake-a-r2",
        "n0001-fake-b-r1",
        "n0001-fake-b-r2",
    ]


def test_guards(campaign_file):
    campaign = config.load_campaign(
        campaign_file(
            scaling={
                "nodes": {"min": 1, "max": 64},
                "max_cells_per_node": 600,
                "min_cells_per_rank": 100,
            }
        )
    )
    result = plan.make_plan(campaign, cells=1000)
    # 1000 cells: 1 node holds too many, 16 and more have too few per rank
    assert result.node_counts == [2, 4, 8]
    assert any(
        "1 nodes" in note and "max_cells_per_node" in note for note in result.notes
    )
    assert any(
        "16 nodes" in note and "min_cells_per_rank" in note for note in result.notes
    )


def test_guards_need_the_mesh(campaign_file):
    campaign = config.load_campaign(
        campaign_file(scaling={"nodes": [1, 2], "min_cells_per_rank": 100})
    )
    result = plan.make_plan(campaign)
    assert result.node_counts == [1, 2]
    assert any("unknown" in note for note in result.notes)


def test_smallest_node_count_from_memory(campaign_file):
    campaign = config.load_campaign(
        campaign_file(
            scaling={"nodes": {"min": "auto", "max": 16}, "max_cells_per_node": 300}
        )
    )
    assert plan.make_plan(campaign, cells=1000).node_counts == [4, 8, 16]
    with pytest.raises(config.ConfigError, match="auto"):
        plan.make_plan(campaign)


def test_definitions_follow_the_settings(campaign_file):
    first = plan.make_plan(config.load_campaign(campaign_file()))
    changed = plan.make_plan(
        config.load_campaign(
            campaign_file(
                variants=[
                    {"name": "a", "env": {"X": 1}},
                    {"name": "b", "env": {"SEISSOL_COMMTHREAD": 0}},
                ]
            )
        )
    )
    definitions = {run.name: run.definition for run in first.runs}
    for run in changed.runs:
        if run.variant == "a":
            assert definitions[run.name] != run.definition
        else:
            assert definitions[run.name] == run.definition


def test_plan_description(campaign_file):
    lines = plan.describe(
        plan.make_plan(config.load_campaign(campaign_file()), cells=1000)
    )
    assert lines[0].startswith("campaign fake: 12 runs in 3 jobs")
    assert any(line.split()[:3] == ["4", "4", "250"] for line in lines)
