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

    def tractionMatrixSpp(self):
        # b = eta * Y is dense for an anisotropic impedance, so every traction row carries all
        # three columns
        tractionMatrixSpp = np.zeros((self.numQuantities(), 3))
        for row in (0, 3, 5):
            tractionMatrixSpp[row, :] = 1
        return tractionMatrixSpp


def kernel_class(**kwargs):
    return AnisotropicADERDG
