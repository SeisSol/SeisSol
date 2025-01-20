# SPDX-FileCopyrightText: 2016-2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

import numpy as np
from kernels.aderdg import LinearADERDG
from yateto import Tensor
from yateto.input import memoryLayoutFromFile, parseXMLMatrixFile


class ViscoelasticADERDG(LinearADERDG):
    def __init__(
        self,
        order,
        multipleSimulations,
        matricesDir,
        memLayout,
        numberOfMechanisms,
        **kwargs
    ):
        self.numberOfMechanisms = numberOfMechanisms
        self.numberOfElasticQuantities = 9

        super().__init__(order, multipleSimulations, matricesDir)
        clones = {
            "star": ["star(0)", "star(1)", "star(2)"],
        }
        self.db.update(
            parseXMLMatrixFile(
                "{}/matrices_viscoelastic.xml".format(matricesDir), clones
            )
        )

        star_spp = self.db.star[0].spp().as_ndarray()
        star_rows, star_cols = star_spp.shape
        aniso_cols = star_cols - self.numberOfElasticQuantities
        star_spp_new = np.zeros(
            (self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool
        )
        star_spp_new[0:star_rows, 0:star_cols] = star_spp
        """ The last 6 columns of star_spp contain the prototype
        sparsity pattern for a mechanism. Therefore, the spp is repeated
        for every mechanism. """
        for mech in range(1, numberOfMechanisms):
            offset0 = self.numberOfElasticQuantities
            offsetm = self.numberOfElasticQuantities + mech * aniso_cols
            star_spp_new[0:star_rows, offsetm : (offsetm + aniso_cols)] = star_spp[
                0:star_rows, offset0 : (offset0 + aniso_cols)
            ]
        for dim in range(3):
            self.db.star[dim] = Tensor(
                self.db.star[dim].name(), star_spp_new.shape, spp=star_spp_new
            )

        source_spp = np.zeros(
            (self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool
        )
        ET_spp = self.db["ET"].spp().as_ndarray()
        """ ET is a prototype sparsity pattern for a mechanism.
        Therefore, repeated for every mechanism. See Kaeser
        and Dumbser 2006, III. Viscoelastic attenuation."""
        for mech in range(numberOfMechanisms):
            offset = self.numberOfElasticQuantities + mech * aniso_cols
            r = slice(offset, offset + aniso_cols)
            source_spp[r, 0:aniso_cols] = ET_spp
            source_spp[r, r] = np.identity(aniso_cols, dtype=bool)
        self.db.ET = Tensor("ET", source_spp.shape, spp=source_spp)

        memoryLayoutFromFile(memLayout, self.db, clones)

        self.kwargs = kwargs

    def numberOfQuantities(self):
        return 9 + 6 * self.numberOfMechanisms

    def starMatrix(self, dim):
        return self.db.star[dim]

    def sourceMatrix(self):
        return self.db.ET

    def name(self):
        return "viscoelastic2"

    def godunov_spp(self):
        spp = np.zeros(
            (self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool
        )
        spp[0 : self.numberOfElasticQuantities, :] = True
        return spp

    def flux_solver_spp(self):
        return self.godunov_spp()

    def transformation_spp(self):
        spp = np.zeros(
            (self.numberOfQuantities(), self.numberOfQuantities()), dtype=bool
        )
        spp[0:6, 0:6] = 1
        spp[6:9, 6:9] = 1
        for mechs in range(self.numberOfMechanisms):
            offset = 9 + mechs * 6
            spp[offset : (offset + 6), offset : (offset + 6)] = 1
        return spp

    def transformation_inv_spp(self):
        return self.transformation_spp()


EQUATION_CLASS = ViscoelasticADERDG
