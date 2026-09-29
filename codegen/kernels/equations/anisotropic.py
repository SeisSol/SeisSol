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
        # three columns. The pattern's shape is the layout's; only the rows it fills are this
        # material's.
        pattern = np.zeros_like(self.extractTractions().T[:, :3])
        for row in (0, 3, 5):
            pattern[row, :] = 1
        return pattern


def kernel_class(**kwargs):
    return AnisotropicADERDG
