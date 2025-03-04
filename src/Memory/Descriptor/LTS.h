// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_LTS_H_
#define SEISSOL_SRC_INITIALIZER_LTS_H_

#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Model/Plasticity.h"
#include "generated_code/tensor.h"
#include <Initializer/CellLocalInformation.h>

#ifdef ACL_DEVICE
#include "Parallel/Helper.h"
#endif

namespace seissol::tensor {
class Qane;
} // namespace seissol::tensor

namespace seissol::initializer {

enum class AllocationPreset {
  Global,
  Timedofs,
  Constant,
  Dofs,
  TimedofsConstant,
  ConstantShared,
  Timebucket,
  Plasticity,
  PlasticityData
};

inline auto allocationModeWP(AllocationPreset preset,
                             int convergenceOrder = seissol::ConvergenceOrder) {
#ifndef ACL_DEVICE
  switch (preset) {
  case seissol::initializer::AllocationPreset::Global:
    [[fallthrough]];
  case seissol::initializer::AllocationPreset::TimedofsConstant:
    return AllocationMode::HostOnlyHBM;
  case seissol::initializer::AllocationPreset::Plasticity:
    return AllocationMode::HostOnly;
  case seissol::initializer::AllocationPreset::PlasticityData:
    return AllocationMode::HostOnly;
  case seissol::initializer::AllocationPreset::Timebucket:
    [[fallthrough]];
  case seissol::initializer::AllocationPreset::Timedofs:
    return (convergenceOrder <= 7 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
  case seissol::initializer::AllocationPreset::Constant:
    [[fallthrough]];
  case seissol::initializer::AllocationPreset::ConstantShared:
    return (convergenceOrder <= 4 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
  case seissol::initializer::AllocationPreset::Dofs:
    return (convergenceOrder <= 3 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
  default:
    return AllocationMode::HostOnly;
  }
#else
  switch (preset) {
  case seissol::initializer::AllocationPreset::Global:
    [[fallthrough]];
  case seissol::initializer::AllocationPreset::Constant:
    [[fallthrough]];
  case seissol::initializer::AllocationPreset::TimedofsConstant:
    return AllocationMode::HostOnly;
  case seissol::initializer::AllocationPreset::Dofs:
    [[fallthrough]];
  case seissol::initializer::AllocationPreset::PlasticityData:
    return useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplitPinned;
  case seissol::initializer::AllocationPreset::Timebucket:
    return useMPIUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit;
  default:
    return useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit;
  }
#endif
}

struct LTS {
  Variable<real[tensor::Q::size()]> dofs;
  // size is zero if Qane is not defined
  Variable<real[zeroLengthArrayHandler(kernels::size<tensor::Qane>())]> dofsAne;
  Variable<real*> buffers;
  Variable<real*> derivatives;
  Variable<CellLocalInformation> cellInformation;
  Variable<SecondaryCellLocalInformation> secondaryInformation;
  Variable<real* [4]> faceNeighbors;
  Variable<LocalIntegrationData> localIntegration;
  Variable<NeighboringIntegrationData> neighboringIntegration;
  Variable<CellMaterialData> material;
  Variable<seissol::model::PlasticityData> plasticity;
  Variable<CellDRMapping[4]> drMapping;
  Variable<CellBoundaryMapping[4]> boundaryMapping;
  Variable<real[tensor::QStress::size() + tensor::QEtaModal::size()]> pstrain;
  Variable<real* [4]> faceDisplacements;
  Bucket buffersDerivatives;
  Bucket faceDisplacementsBuffer;

  Variable<real*> buffersDevice;
  Variable<real*> derivativesDevice;
  Variable<real* [4]> faceDisplacementsDevice;
  Variable<real* [4]> faceNeighborsDevice;
  Variable<CellDRMapping[4]> drMappingDevice;
  Variable<CellBoundaryMapping[4]> boundaryMappingDevice;

#ifdef ACL_DEVICE
  ScratchpadMemory integratedDofsScratch;
  ScratchpadMemory derivativesScratch;
  ScratchpadMemory nodalAvgDisplacements;
  ScratchpadMemory analyticScratch;
#endif

  /// \todo Memkind
  void addTo(LTSTree& tree, bool usePlasticity) {
    LayerMask plasticityMask;
    if (usePlasticity) {
      plasticityMask = LayerMask(Ghost);
    } else {
      plasticityMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    }

    tree.addVar(dofs, LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));
    if (kernels::size<tensor::Qane>() > 0) {
      tree.addVar(
          dofsAne, LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));
    }
    tree.addVar(
        buffers, LayerMask(), 1, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    tree.addVar(
        derivatives, LayerMask(), 1, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    tree.addVar(
        cellInformation, LayerMask(), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.addVar(secondaryInformation, LayerMask(), 1, AllocationMode::HostOnly, true);
    tree.addVar(faceNeighbors,
                LayerMask(Ghost),
                1,
                allocationModeWP(AllocationPreset::TimedofsConstant),
                true);
    tree.addVar(localIntegration,
                LayerMask(Ghost),
                1,
                allocationModeWP(AllocationPreset::ConstantShared),
                true);
    tree.addVar(neighboringIntegration,
                LayerMask(Ghost),
                1,
                allocationModeWP(AllocationPreset::ConstantShared),
                true);
    tree.addVar(material, LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.addVar(
        plasticity, plasticityMask, 1, allocationModeWP(AllocationPreset::Plasticity), true);
    tree.addVar(drMapping, LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.addVar(
        boundaryMapping, LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.addVar(
        pstrain, plasticityMask, PagesizeHeap, allocationModeWP(AllocationPreset::PlasticityData));
    tree.addVar(faceDisplacements, LayerMask(Ghost), PagesizeHeap, AllocationMode::HostOnly, true);

    // TODO(David): remove/rename "constant" flag (the data is temporary; and copying it for IO is
    // handled differently)
    tree.addBucket(
        buffersDerivatives, PagesizeHeap, allocationModeWP(AllocationPreset::Timebucket), true);
    tree.addBucket(
        faceDisplacementsBuffer, PagesizeHeap, allocationModeWP(AllocationPreset::Timedofs));

    tree.addVar(buffersDevice, LayerMask(), 1, AllocationMode::HostOnly, true);
    tree.addVar(derivativesDevice, LayerMask(), 1, AllocationMode::HostOnly, true);
    tree.addVar(faceDisplacementsDevice, LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.addVar(faceNeighborsDevice, LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.addVar(drMappingDevice, LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.addVar(boundaryMappingDevice, LayerMask(Ghost), 1, AllocationMode::HostOnly, true);

#ifdef ACL_DEVICE
    tree.addScratchpadMemory(integratedDofsScratch, 1, AllocationMode::HostDeviceSplit);
    tree.addScratchpadMemory(derivativesScratch, 1, AllocationMode::DeviceOnly);
    tree.addScratchpadMemory(nodalAvgDisplacements, 1, AllocationMode::DeviceOnly);
    tree.addScratchpadMemory(analyticScratch, 1, AllocationMode::HostDevicePinned);
#endif
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) const {
    manager.registerData("dofs", tree, dofs);
    if constexpr (kernels::size<tensor::Qane>() > 0) {
      manager.registerData("dofsAne", tree, dofsAne);
    }
    // check plasticity usage over the layer mask (for now)
    if (plasticity.mask == LayerMask(Ghost)) {
      manager.registerData("pstrain", tree, pstrain);
    }
  }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_LTS_H_
