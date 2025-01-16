# SPDX-FileCopyrightText: 2016-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff


def getArchitectures():
    # wsm = Westmere
    # snb = Sandy Bridge
    # knc = Knights Corner (Xeon Phi)
    # hsw = Haswell
    # knl = Knight Landing (Xeon Phi)
    cpus = [
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
    precisions = ["s", "d"]
    return [p + c for c in cpus for p in precisions]


def getRealSize(architecture):
    realSize = {"s": 4, "d": 8}
    return realSize[architecture[0].lower()]


def getCpu(architecture):
    return architecture[1:]
