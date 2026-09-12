# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff
# SPDX-FileContributor: Sebastian Wolf

from kernels.equations.elastic import ElasticADERDG


class AnisotropicADERDG(ElasticADERDG):
    def name(self):
        return "anisotropic"


def kernel_class(**kwargs):
    return AnisotropicADERDG
