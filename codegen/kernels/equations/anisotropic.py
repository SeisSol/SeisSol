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

    def fluxDecomposes(self):
        """The scalars a face would carry are what an isotropic Riemann problem
        leaves. An anisotropic one has no such split: over a hundred material
        pairs with a five percent anisotropic perturbation its operator occupies
        twenty-five of eighty-one entries instead of thirteen and spans sixteen
        dimensions instead of ten, and reading the ten off it leaves six percent
        of the operator behind. A decomposition for it would be written from the
        stiffness tensor, not from these scalars.
        """
        return False

    def tractionMatrixSpp(self):
        # b = eta * Y is dense for an anisotropic impedance, so every traction row carries all
        # three columns
        tractionMatrixSpp = np.zeros((self.numQuantities(), 3))
        for row in (0, 3, 5):
            tractionMatrixSpp[row, :] = 1
        return tractionMatrixSpp


def kernel_class(**kwargs):
    return AnisotropicADERDG
