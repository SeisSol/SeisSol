// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITMODEL_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITMODEL_H_

namespace seissol {
class SeisSol;
} // namespace seissol

namespace seissol::initializer::initprocedure {
void initModel(seissol::SeisSol& seissolInstance);
} // namespace seissol::initializer::initprocedure

#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INITMODEL_H_
