# SPDX-FileCopyrightText: 2016-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause


# @file
# This file is part of SeisSol.
#
# @author Carsten Uphoff (c.uphoff AT tum.de)
#

import os
import re

import kernels.arch as arch


class Candidate(object):
    """A Candidate measures if a memory layout is suitable
    for a build configuration. If a build configuration
    shares an attribute with a Candidate, the Candidate
    gets a higher score.
    The scoring system is chosen such that the best
    Candidate is unique (i.e. 2**Importance).
    """

    IMPORTANCE = {
        "precision": 1,
        "equations": 2,
        "order": 3,
        "pe": 4,
        "multipleSimulations": 5,
    }

    def __init__(self, atts):
        self.atts = atts

    def score(self, reqs):
        sc = 0
        for key, val in reqs.items():
            if key in self.atts:
                if val == self.atts[key]:
                    sc += 2 ** self.IMPORTANCE[key]
                else:  # requirement not satisifed
                    return 0
        return sc

    def __repr__(self):
        return repr(self.atts)


def findCandidates(search_path):
    """Determine Candidate attributes from file name."""

    archs = arch.getArchitectures()
    pes = [arch.getCpu(a) for a in archs]

    candidates = dict()
    for c in os.listdir(search_path):
        name, ext = os.path.splitext(c)
        atts = dict()
        for att in name.split("_"):
            multipleSimulations = re.match("ms([0-9]+)", att)
            order = re.match("O([0-9]+)", att)
            if multipleSimulations:
                atts["multipleSimulations"] = int(multipleSimulations.group(1))
            elif order:
                atts["order"] = int(order.group(1))
            elif att.lower() in ["s", "d"]:
                atts["precision"] = att.lower()
            elif att.lower() in pes:
                atts["pe"] = att.lower()
            else:
                atts["equations"] = att
        candidates[c] = Candidate(atts)
    return candidates


def guessMemoryLayout(env):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(script_dir, "..", "config")

    if "gpu" in env["targets"]:
        print(
            "INFO: Found gpu as a target. "
            + "Memory layout will fall back to all dense"
        )
        return os.path.join(path, "dense.xml")

    values = {
        "precision": env["arch"][0].lower(),
        "equations": env["equations"].lower(),
        "order": int(env["order"]),
        "pe": arch.getCpu(env["arch"]),
        "multipleSimulations": int(env["multipleSimulations"]),
    }

    candidates = findCandidates(search_path=path)
    bestFit = max(candidates.keys(), key=lambda key: candidates[key].score(values))
    bestScore = candidates[bestFit].score(values)

    if bestScore == 0:
        print(
            "WARNING: No suitable memory layout found."
            + "(Will fall back to all dense.)"
        )
        bestFit = "dense.xml"
    print("Using memory layout {}".format(bestFit))
    return os.path.join(path, bestFit)
