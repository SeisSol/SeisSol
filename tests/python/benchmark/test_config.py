# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The campaign file and the machine profiles are read strictly."""

import json
from pathlib import Path

import pytest
from strongscaling import config

MACHINES = Path(config.MACHINE_DIRECTORY)


def test_durations():
    assert config.parse_duration(90, "t") == 5400
    assert config.parse_duration("01:30:00", "t") == 5400
    assert config.parse_duration("1-00:00", "t") == 86400
    assert config.format_duration(5400) == "01:30:00"
    assert config.format_duration(90000) == "1-01:00:00"
    with pytest.raises(config.ConfigError):
        config.parse_duration("ninety", "t")
    with pytest.raises(config.ConfigError):
        config.parse_duration(0, "t")


def test_sizes():
    assert config.parse_size("96G", "m") == 96 * 2**30
    assert config.parse_size("512 GiB", "m") == 512 * 2**30
    assert config.parse_size(1000, "m") == 1000
    with pytest.raises(config.ConfigError):
        config.parse_size("lots", "m")


@pytest.mark.parametrize("path", sorted(MACHINES.glob("*.yaml")), ids=lambda p: p.stem)
def test_shipped_profiles_parse(path):
    pytest.importorskip("yaml")
    machine = config.load_machine(path.stem, Path("."))
    assert machine.name == path.stem
    layout = machine.layout
    assert layout.threads_per_rank <= layout.cpus_per_rank
    assert layout.ranks_per_node * layout.cpus_per_rank <= (
        machine.cores_per_node * machine.threads_per_core
    )


def test_a_campaign(campaign_file):
    campaign = config.load_campaign(campaign_file())
    assert campaign.name == "fake"
    assert campaign.case.overrides == {"abortcriteria.endtime": "1.0"}
    assert campaign.case.parameter_file.name == "parameters.par"
    assert campaign.case.parameter_file.is_absolute()
    assert [variant.name for variant in campaign.variants] == ["a", "b"]
    assert campaign.variants[1].env == {"SEISSOL_COMMTHREAD": "0"}
    assert campaign.scaling.nodes == [1, 2, 4]
    assert campaign.job.time_limit == 600
    assert campaign.directory == campaign_file().parent / "fake"


def test_unknown_keys_are_errors(campaign_file):
    path = campaign_file(scaling={"nodes": [1, 2], "repetition": 2})
    with pytest.raises(config.ConfigError, match="repetition"):
        config.load_campaign(path)


def test_the_output_is_the_campaigns(campaign_file):
    path = campaign_file(
        case={
            "parameter_file": "case/parameters.par",
            "set": {"output.outputfile": "x"},
        }
    )
    with pytest.raises(config.ConfigError, match="outputfile"):
        config.load_campaign(path)


def test_a_build_is_built_or_given(campaign_file):
    with pytest.raises(config.ConfigError, match="either binary"):
        config.load_campaign(campaign_file(builds=[{"name": "x"}]))


def test_node_ranges(campaign_file):
    campaign = config.load_campaign(
        campaign_file(scaling={"nodes": {"min": 2, "max": 32}})
    )
    assert campaign.scaling.node_range == (2, 32, 2)
    with pytest.raises(config.ConfigError):
        config.load_campaign(campaign_file(scaling={"nodes": {"min": 8, "max": 2}}))


def test_variant_layouts(campaign_file):
    machine = config.load_campaign(campaign_file()).machine
    # half the ranks of a node get twice the CPUs, and the workers follow
    variant = config.Variant(name="v", ranks_per_node=2)
    base = config.Variant(name="w")
    assert base.layout(machine).to_dict() == {
        "ranks_per_node": 1,
        "threads_per_rank": 3,
        "cpus_per_rank": 4,
        "gpus_per_rank": 0,
    }
    layout = variant.layout(machine)
    assert (layout.ranks_per_node, layout.cpus_per_rank, layout.threads_per_rank) == (
        2,
        2,
        1,
    )


def test_profiles_merge(tmp_path):
    profile = tmp_path / "machine.json"
    profile.write_text(
        json.dumps(
            {
                "scheduler": "slurm",
                "node": {"cores": 8},
                "slurm": {"partition": [{"max_nodes": 4, "partition": "small"}]},
            }
        )
    )
    machine = config.load_machine(
        {"profile": str(profile), "slurm": {"account": "project"}}, tmp_path
    )
    assert machine.slurm["account"] == "project"
    assert config.partition_for(machine.slurm, 4) == "small"
    with pytest.raises(config.ConfigError, match="no partition"):
        config.partition_for(machine.slurm, 5)
