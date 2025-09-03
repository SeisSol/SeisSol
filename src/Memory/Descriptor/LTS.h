// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_
#define SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_

#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Model/Plasticity.h"
#include "GeneratedCode/tensor.h"
#include <Initializer/CellLocalInformation.h>
#include <Parallel/Helper.h>

#ifdef ACL_DEVICE
#include "Parallel/Helper.h"
#endif

namespace seissol::tensor {
class Qane;
} // namespace seissol::tensor

namespace seissol::initializer {

struct LTS {
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

  static auto allocationModeWP(AllocationPreset preset,
                               int convergenceOrder = seissol::ConvergenceOrder) {
    if constexpr (!isDeviceOn()) {
      switch (preset) {
      case AllocationPreset::Global:
        [[fallthrough]];
      case AllocationPreset::TimedofsConstant:
        return AllocationMode::HostOnlyHBM;
      case AllocationPreset::Plasticity:
        return AllocationMode::HostOnly;
      case AllocationPreset::PlasticityData:
        return AllocationMode::HostOnly;
      case AllocationPreset::Timebucket:
        [[fallthrough]];
      case AllocationPreset::Timedofs:
        return (convergenceOrder <= 7 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
      case AllocationPreset::Constant:
        [[fallthrough]];
      case AllocationPreset::ConstantShared:
        return (convergenceOrder <= 4 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
      case AllocationPreset::Dofs:
        return (convergenceOrder <= 3 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
      default:
        return AllocationMode::HostOnly;
      }
    } else {
      switch (preset) {
      case AllocationPreset::Global:
        [[fallthrough]];
      case AllocationPreset::Constant:
        [[fallthrough]];
      case AllocationPreset::TimedofsConstant:
        return AllocationMode::HostOnly;
      case AllocationPreset::Dofs:
        [[fallthrough]];
      case AllocationPreset::PlasticityData:
        return useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplitPinned;
      case AllocationPreset::Timebucket:
        return useMPIUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit;
      default:
        return useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit;
      }
    }
  }
  Variable<real[tensor::Q::size()]> dofs;
  // size is zero if Qane is not defined
  Variable<real[zeroLengthArrayHandler(kernels::size<tensor::Qane>())]> dofsAne;
  Variable<real*> buffers;
  Variable<real*> derivatives;
  Variable<CellLocalInformation> cellInformation;
  Variable<SecondaryCellLocalInformation> secondaryInformation;
  Variable<real* [Cell::NumFaces]> faceNeighbors;
  Variable<LocalIntegrationData> localIntegration;
  Variable<NeighboringIntegrationData> neighboringIntegration;
  Variable<CellMaterialData> material;
  Variable<seissol::model::PlasticityData> plasticity;
  Variable<CellDRMapping[Cell::NumFaces]> drMapping;
  Variable<CellBoundaryMapping[Cell::NumFaces]> boundaryMapping;
  Variable<real[tensor::QStress::size() + tensor::QEtaModal::size()]> pstrain;
  Variable<real* [Cell::NumFaces]> faceDisplacements;
  Bucket<real> buffersDerivatives;

  Variable<real*> buffersDevice;
  Variable<real*> derivativesDevice;
  Variable<real* [Cell::NumFaces]> faceDisplacementsDevice;
  Variable<real* [Cell::NumFaces]> faceNeighborsDevice;
  Variable<CellDRMapping[Cell::NumFaces]> drMappingDevice;
  Variable<CellBoundaryMapping[Cell::NumFaces]> boundaryMappingDevice;

  Scratchpad<real> integratedDofsScratch;
  Scratchpad<real> derivativesScratch;
  Scratchpad<real> nodalAvgDisplacements;
  Scratchpad<real> analyticScratch;
  Scratchpad<real> derivativesExtScratch;
  Scratchpad<real> derivativesAneScratch;
  Scratchpad<real> idofsAneScratch;
  Scratchpad<real> dofsExtScratch;

  Scratchpad<unsigned> flagScratch;
  Scratchpad<real> prevDofsScratch;
  Scratchpad<real> qEtaNodalScratch;
  Scratchpad<real> qStressNodalScratch;

  Scratchpad<real> rotateDisplacementToFaceNormalScratch;
  Scratchpad<real> rotateDisplacementToGlobalScratch;
  Scratchpad<real> rotatedFaceDisplacementScratch;
  Scratchpad<real> dofsFaceNodalScratch;
  Scratchpad<real> prevCoefficientsScratch;
  Scratchpad<real> dofsFaceBoundaryNodalScratch;

  void addTo(LTSTree& tree, bool usePlasticity) {
    LayerMask plasticityMask;
    if (usePlasticity) {
      plasticityMask = LayerMask(Ghost);
    } else {
      plasticityMask = LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior);
    }

    tree.add(dofs, LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));
    if (kernels::size<tensor::Qane>() > 0) {
      tree.add(dofsAne, LayerMask(Ghost), PagesizeHeap, allocationModeWP(AllocationPreset::Dofs));
    } else {
      tree.add(dofsAne,
               LayerMask(Ghost) | LayerMask(Copy) | LayerMask(Interior),
               PagesizeHeap,
               allocationModeWP(AllocationPreset::Dofs));
    }
    tree.add(buffers, LayerMask(), 1, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    tree.add(
        derivatives, LayerMask(), 1, allocationModeWP(AllocationPreset::TimedofsConstant), true);
    tree.add(cellInformation, LayerMask(), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.add(secondaryInformation, LayerMask(), 1, AllocationMode::HostOnly, true);
    tree.add(faceNeighbors,
             LayerMask(Ghost),
             1,
             allocationModeWP(AllocationPreset::TimedofsConstant),
             true);
    tree.add(localIntegration,
             LayerMask(Ghost),
             1,
             allocationModeWP(AllocationPreset::ConstantShared),
             true);
    tree.add(neighboringIntegration,
             LayerMask(Ghost),
             1,
             allocationModeWP(AllocationPreset::ConstantShared),
             true);
    tree.add(material, LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.add(plasticity, plasticityMask, 1, allocationModeWP(AllocationPreset::Plasticity), true);
    tree.add(drMapping, LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.add(
        boundaryMapping, LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.add(
        pstrain, plasticityMask, PagesizeHeap, allocationModeWP(AllocationPreset::PlasticityData));
    tree.add(faceDisplacements, LayerMask(Ghost), PagesizeHeap, AllocationMode::HostOnly, true);

    // TODO(David): remove/rename "constant" flag (the data is temporary; and copying it for IO is
    // handled differently)
    tree.add(buffersDerivatives,
             LayerMask(),
             PagesizeHeap,
             allocationModeWP(AllocationPreset::Timebucket),
             true);

    tree.add(buffersDevice, LayerMask(), 1, AllocationMode::HostOnly, true);
    tree.add(derivativesDevice, LayerMask(), 1, AllocationMode::HostOnly, true);
    tree.add(faceDisplacementsDevice, LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.add(faceNeighborsDevice, LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.add(drMappingDevice, LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.add(boundaryMappingDevice, LayerMask(Ghost), 1, AllocationMode::HostOnly, true);

    if constexpr (isDeviceOn()) {
      tree.add(derivativesExtScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(derivativesAneScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(idofsAneScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(dofsExtScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(integratedDofsScratch, LayerMask(), 1, AllocationMode::HostDeviceSplit);
      tree.add(derivativesScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(nodalAvgDisplacements, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(analyticScratch, LayerMask(), 1, AllocationMode::HostDevicePinned);

      tree.add(flagScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(prevDofsScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(qEtaNodalScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(qStressNodalScratch, LayerMask(), 1, AllocationMode::DeviceOnly);

      tree.add(rotateDisplacementToFaceNormalScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(rotateDisplacementToGlobalScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(rotatedFaceDisplacementScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(dofsFaceNodalScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(prevCoefficientsScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
      tree.add(dofsFaceBoundaryNodalScratch, LayerMask(), 1, AllocationMode::DeviceOnly);
    }
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) const {
    manager.registerData("dofs", tree, dofs);
    if constexpr (kernels::size<tensor::Qane>() > 0) {
      manager.registerData("dofsAne", tree, dofsAne);
    }
    // check plasticity usage over the layer mask (for now)
    if (tree->info(plasticity).mask == LayerMask(Ghost)) {
      manager.registerData("pstrain", tree, pstrain);
    }
  }
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_LTS_H_
