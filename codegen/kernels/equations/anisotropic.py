# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff
# SPDX-FileContributor: Sebastian Wolf

import numpy as np
from kernels.equations.elastic import ElasticADERDG


class AnisotropicADERDG(ElasticADERDG):
    def name(self):
        return "anisotropic"

    def extractTractions(self):
        # The pattern is dense across the traction triple here, because an
        # anisotropic face couples all three of them; only the width follows
        # the layout.
        extractTractionsSPP = np.zeros((3, self.numTransportQuantities()))
        extractTractionsSPP[0, 0] = 1
        extractTractionsSPP[1, 0] = 1
        extractTractionsSPP[2, 0] = 1
        extractTractionsSPP[0, 3] = 1
        extractTractionsSPP[1, 3] = 1
        extractTractionsSPP[2, 3] = 1
        extractTractionsSPP[0, 5] = 1
        extractTractionsSPP[1, 5] = 1
        extractTractionsSPP[2, 5] = 1
        return extractTractionsSPP


def kernel_class(**kwargs):
    return AnisotropicADERDG
