# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Strong scaling series of SeisSol on a cluster, from one campaign file.

A campaign names a case (a parameter file), the builds and run-time variants to
compare, the node counts, and the machine; the modules plan its runs, build and
check SeisSol, submit one job per node count, and follow the jobs. See
docs/strong-scaling.rst.
"""
