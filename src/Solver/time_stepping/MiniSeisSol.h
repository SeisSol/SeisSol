// SPDX-FileCopyrightText: 2017-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_SOLVER_TIME_STEPPING_MINISEISSOL_H_
#define SEISSOL_SRC_SOLVER_TIME_STEPPING_MINISEISSOL_H_

#include "Initializer/MemoryManager.h"

namespace seissol {
  void localIntegration(GlobalData* globalData,
                        initializer::LTS& lts,
                        initializer::Layer& layer,
                        seissol::SeisSol& seissolInstance);

  void localIntegrationOnDevice(CompoundGlobalData& globalData,
                                initializer::LTS& lts,
                                initializer::Layer& layer,
                                seissol::SeisSol& seissolInstance,
                                seissol::parallel::runtime::StreamRuntime& runtime);
  
  void fillWithStuff(real* buffer,
                     unsigned nValues);

  void fakeData(initializer::LTS& lts,
                initializer::Layer& layer,
                FaceType faceTp = FaceType::Regular);
  
  double miniSeisSol(initializer::MemoryManager& memoryManager,
                     bool usePlasticity,
                     seissol::SeisSol& seissolInstance);
  constexpr real miniSeisSolTimeStep = 1.0;
} //namespace seissol


#endif // SEISSOL_SRC_SOLVER_TIME_STEPPING_MINISEISSOL_H_

