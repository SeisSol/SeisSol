# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

import os
import re


class Candidate(object):
    """Measures if a memory layout is suitable for a build configuration.

    If a build configuration
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
        "gemmgen": 5,
        "multipleSimulations": 6,
    }

    def __init__(self, atts):
        self.atts = atts

    def score(self, reqs):
        sc = 0
        for key, val in reqs.items():
            if key in self.atts:
                if val == self.atts[key] or (
                    key == "gemmgen" and self.atts[key] in val
                ):
                    sc += 2 ** self.IMPORTANCE[key]
                else:  # requirement not satisifed
                    return 0
        return sc

    def __repr__(self):
        return repr(self.atts)


def findCandidates(search_path):
    """Determine Candidate attributes from file name."""

    # for now, keep a small subset of possible (old) archs here
    pes = [
        "noarch",
        "wsm",
        "snb",
        "knc",
        "hsw",
        "knl",
        "skx",
        "thunderx2t99",
        "a64fx",
    ]
    gemmgen = ["pspamm", "libxsmm", "tensorforge"]

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
            elif att.lower() in ["s", "d", "f32", "f64"]:
                atts["precision"] = att.lower()
            elif att.lower() in pes:
                atts["pe"] = att.lower()
            elif att.lower() in gemmgen:
                atts["gemmgen"] = att.lower()
            else:
                atts["equations"] = att
        candidates[c] = Candidate(atts)
    return candidates


def guessMemoryLayout(env):
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # look at the CPU or GPU configs; dependent on the considered target
    subfolder = "gpu" if "gpu" in env["targets"] else "cpu"
    path = os.path.join(script_dir, "..", "config", subfolder)

    values = {
        "precision": env["precision"].lower(),
        "equations": env["equations"].lower(),
        "order": int(env["order"]),
        "pe": env["arch"].lower(),
        "multipleSimulations": int(env["multipleSimulations"]),
        "gemmgen": set(gg.lower() for gg in env["gemmgen"]),
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
