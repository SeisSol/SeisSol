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
    QuantityGroup("pprime", QuantityKind.SCALAR, FaceRole.TRACTION),
    QuantityGroup("v", QuantityKind.VECTOR, FaceRole.VELOCITY),
]

MechanismGroups = [QuantityGroup("theta", QuantityKind.SCALAR)]

ViscoacousticADERDG, ViscoacousticAnelasticADERDG = visco_classes(
    "viscoacoustic", PrimaryGroups, MechanismGroups
)


def kernel_class(**kwargs):
    return select(
        ViscoacousticADERDG, ViscoacousticAnelasticADERDG, kwargs["visco_mode"]
    )
