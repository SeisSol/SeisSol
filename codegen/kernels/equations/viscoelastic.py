# SPDX-FileCopyrightText: 2016 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff
from kernels.equations.visco import select, visco_classes
from kernels.quantities import FaceRole, QuantityGroup, QuantityKind

PrimaryGroups = [
    QuantityGroup("s", QuantityKind.SYM_TENSOR2, FaceRole.TRACTION),
    QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
]

MechanismGroups = [QuantityGroup("theta", QuantityKind.SYM_TENSOR2)]

ViscoelasticADERDG, ViscoelasticAnelasticADERDG = visco_classes(
    "viscoelastic", PrimaryGroups, MechanismGroups
)


def kernel_class(**kwargs):
    return select(ViscoelasticADERDG, ViscoelasticAnelasticADERDG, kwargs["visco_mode"])
