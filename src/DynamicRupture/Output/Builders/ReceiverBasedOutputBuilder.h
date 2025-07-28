// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_BUILDERS_RECEIVERBASEDOUTPUTBUILDER_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_BUILDERS_RECEIVERBASEDOUTPUTBUILDER_H_

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "DynamicRupture/Output/OutputAux.h"
#include "Geometry/MeshReader.h"
#include "Initializer/InputAux.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Descriptor/LTS.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"
#include "Parallel/MPI.h"
#include <Memory/Tree/Backmap.h>
#include <memory>
#include <vector>

namespace seissol::dr::output {
class ReceiverBasedOutputBuilder {
  public:
  virtual ~ReceiverBasedOutputBuilder() = default;
  virtual void build(std::shared_ptr<ReceiverOutputData> outputData) = 0;

  void setMeshReader(const seissol::geometry::MeshReader* reader);
  void setLtsData(LTS::Storage& userWpTree,
                  LTS::Backmap& userWpLut,
                  DynamicRupture::Storage& userDrTree);

  void setVariableList(const std::vector<std::size_t>& variables);
  void setFaceToLtsMap(std::vector<std::size_t>* faceToLtsMap);

  protected:
  virtual void initTimeCaching() = 0;

  void initBasisFunctions();
  void initFaultDirections();
  void initRotationMatrices();
  void initOutputVariables(std::array<bool, std::tuple_size<DrVarsT>::value>& outputMask);
  void initJacobian2dMatrices();
  void assignNearestInternalGaussianPoints();
  void assignFaultTags();
  void assignFusedIndices();

  const seissol::geometry::MeshReader* meshReader{};
  LTS::Storage* wpTree;
  LTS::Backmap* wpLut;
  DynamicRupture::Storage* drTree;
  std::shared_ptr<ReceiverOutputData> outputData;
  std::vector<std::size_t> variables;
  std::vector<std::size_t>* faceToLtsMap{nullptr};
  int localRank{-1};
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_BUILDERS_RECEIVERBASEDOUTPUTBUILDER_H_
