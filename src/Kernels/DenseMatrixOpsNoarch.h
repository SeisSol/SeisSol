// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_DENSEMATRIXOPSNOARCH_H_
#define SEISSOL_SRC_KERNELS_DENSEMATRIXOPSNOARCH_H_

#define DMO_INCREMENT64 1
#define DMO_STREAM64(IN, OUT) *(OUT) = *(IN);
#define DMO_INCREMENT32 1
#define DMO_STREAM64(IN, OUT) *(OUT) = *(IN);

#endif // SEISSOL_SRC_KERNELS_DENSEMATRIXOPSNOARCH_H_
