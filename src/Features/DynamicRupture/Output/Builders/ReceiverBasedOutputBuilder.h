// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_BUILDERS_RECEIVERBASEDOUTPUTBUILDER_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_BUILDERS_RECEIVERBASEDOUTPUTBUILDER_H_

#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/InputAux.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Backmap.h"
#include "Core/Numerical/Transformation.h"
#include "Core/Parallel/MPI.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Output/DataTypes.h"
#include "Features/DynamicRupture/Output/OutputAux.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "WavePropagation/Physics/Common.h"

#include <memory>
#include <vector>

namespace seissol::dr::output {
class ReceiverBasedOutputBuilder {
  public:
  virtual ~ReceiverBasedOutputBuilder() = default;

  void setMeshReader(const seissol::geometry::MeshReader* reader);
  void setLtsData(LTS::Storage& userWpStorage,
                  LTS::Backmap& userWpBackmap,
                  DynamicRupture::Storage& userDrStorage);

  void setVariableList(const std::vector<std::size_t>& variables);
  void setFaceToLtsMap(::seissol::initializer::StorageBackmap<1>* faceToLtsMap);
  /// The simulations of the run, of which the faces fuse `fusedSimulations` each; there is a
  /// receiver for each simulation.
  void setSimulationCount(std::size_t numSimulations, std::size_t fusedSimulations) {
    numSimulations_ = numSimulations;
    fusedSimulations_ = fusedSimulations;
  }

  protected:
  virtual void initTimeCaching() = 0;

  void initTopology();
  void initBasisFunctions();
  void initDeviceCollectors(bool elementwise);
  void initFaultDirections();
  void initRotationMatrices();
  void initOutputVariables(std::array<bool, std::tuple_size_v<DrVarsT>>& outputMask);
  void initJacobian2dMatrices();
  void assignNearestGaussianPoints();
  void assignNearestInternalGaussianPoints();
  void assignFaultTags();
  void assignFusedIndices();

  /// The configuration a fault face computes in.
  [[nodiscard]] ConfigId configOfFace(std::size_t faultFaceIndex) const;

  /// The copy of the fault faces which holds the simulation of a receiver.
  void assignCopies();

  const seissol::geometry::MeshReader* meshReader_{nullptr};
  LTS::Storage* wpStorage_{nullptr};
  LTS::Backmap* wpBackmap_{nullptr};
  DynamicRupture::Storage* drStorage_{nullptr};
  std::shared_ptr<ReceiverOutputData> outputData_;
  std::vector<std::size_t> variables_;
  ::seissol::initializer::StorageBackmap<1>* faceToLtsMap_{nullptr};
  std::size_t numSimulations_{1};
  std::size_t fusedSimulations_{1};
  int localRank_{-1};
};
} // namespace seissol::dr::output

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_OUTPUT_BUILDERS_RECEIVERBASEDOUTPUTBUILDER_H_
