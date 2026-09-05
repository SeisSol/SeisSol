# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

"""Equation classes for materials with relaxation.

Adding attenuation to a material does not depend on which material it is: the
mechanism block is repeated along the quantity axis, the star matrix grows one
prototype column block per mechanism, and the source picks up the relaxation
diagonal. The only things that differ are the base material's quantity groups
and what one mechanism carries, so both are arguments here rather than a
second copy of the file.
"""

import numpy as np
from kernels.aderdg.linearck import LinearCK
from kernels.aderdg.linearckanelastic import LinearCKAnelastic
from kernels.quantities import layout, total_extent
from yateto import Tensor


def visco_classes(equation_name, primary_groups, mechanism_groups):
    """Returns the fused and the split equation class for a relaxing material.

    The fused class carries the memory variables in Q and needs the star and
    source patterns widened to match. The split class keeps them in their own
    tensor, so its layout is the base material's and nothing has to be widened.
    """

    num_elastic = total_extent(layout(primary_groups))
    num_per_mechanism = total_extent(layout(mechanism_groups))

    class Fused(LinearCK):
        def __init__(
            self,
            order,
            multipleSimulations,
            matricesDir,
            memLayout,
            numMechanisms,
            **kwargs,
        ):
            self.numMechanisms = numMechanisms
            self.numElasticQuantities = num_elastic

            super().__init__(order, multipleSimulations, matricesDir)
            clones = dict(self.StarClones)
            self.db.update(self.readMatrices(matricesDir, clones))

            self._widenStar()
            self._widenSource()

            self.finishConfigure(memLayout, clones, kwargs)

        def _widenStar(self):
            """The trailing columns of the star hold one mechanism's prototype
            pattern, so repeat them for every further mechanism."""
            star_spp = self.db.star[0].spp().as_ndarray()
            star_rows, star_cols = star_spp.shape
            widened = np.zeros((self.numQuantities(), self.numQuantities()), dtype=bool)
            widened[0:star_rows, 0:star_cols] = star_spp
            for mech in range(1, self.numMechanisms):
                source = slice(num_elastic, num_elastic + num_per_mechanism)
                target = slice(
                    num_elastic + mech * num_per_mechanism,
                    num_elastic + (mech + 1) * num_per_mechanism,
                )
                widened[0:star_rows, target] = star_spp[0:star_rows, source]
            for dim in range(3):
                self.db.star[dim] = Tensor(
                    self.db.star[dim].name(), widened.shape, spp=widened
                )

        def _widenSource(self):
            """ET holds one mechanism's prototype pattern. Every mechanism gets
            a copy of it, plus its own relaxation on the diagonal. See Kaeser
            and Dumbser 2006, III. Viscoelastic attenuation."""
            widened = np.zeros((self.numQuantities(), self.numQuantities()), dtype=bool)
            prototype = self.db["ET"].spp().as_ndarray()
            for mech in range(self.numMechanisms):
                offset = num_elastic + mech * num_per_mechanism
                rows = slice(offset, offset + num_per_mechanism)
                widened[rows, 0:num_per_mechanism] = prototype
                widened[rows, rows] = np.identity(num_per_mechanism, dtype=bool)
            self.db.ET = Tensor("ET", widened.shape, spp=widened)

        def primaryGroups(self):
            return list(primary_groups)

        def mechanismGroups(self):
            return list(mechanism_groups)

        def mechanismRepetitions(self):
            return self.numMechanisms

        def sourceMatrix(self):
            return self.db.ET

        def name(self):
            return equation_name

        def godunov_spp(self):
            spp = np.zeros((self.numQuantities(), self.numQuantities()), dtype=bool)
            spp[0:num_elastic, :] = True
            return spp

        def flux_solver_spp(self):
            return self.godunov_spp()

    class Split(LinearCKAnelastic):
        def __init__(
            self,
            order,
            multipleSimulations,
            matricesDir,
            memLayout,
            numMechanisms,
            **kwargs,
        ):
            super().__init__(
                order,
                multipleSimulations,
                matricesDir,
                memLayout,
                numMechanisms,
                **kwargs,
            )
            self.configure(matricesDir, memLayout, kwargs)

        def primaryGroups(self):
            return list(primary_groups)

        def mechanismGroups(self):
            return list(mechanism_groups)

        def name(self):
            return equation_name

    stem = equation_name.capitalize()
    Fused.__name__ = Fused.__qualname__ = f"{stem}ADERDG"
    Split.__name__ = Split.__qualname__ = f"{stem}AnelasticADERDG"
    return Fused, Split


def select(fused, split, visco_mode):
    """Picks the class the requested VISCO_MODE asks for."""
    mode = visco_mode.lower()
    if mode == "extend":
        return fused
    if mode == "split":
        return split
    raise NotImplementedError(f"Unknown visco_mode {visco_mode}.")
