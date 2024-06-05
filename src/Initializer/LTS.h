/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/

#ifndef INITIALIZER_LTS_H_
#define INITIALIZER_LTS_H_

#include "Initializer/tree/LTSTree.hpp"
#include "Initializer/tree/Layer.hpp"
#include "Initializer/typedefs.hpp"
#include "Kernels/common.hpp"
#include "Model/plasticity.hpp"
#include "generated_code/tensor.h"
#include "tree/Layer.hpp"

namespace seissol {
namespace tensor {
class Qane;
}
} // namespace seissol

namespace seissol::initializer {

enum class AllocationPreset {
  Global,
  Timedofs,
  Constant,
  Dofs,
  TimedofsConstant,
  Timebucket,
  Plasticity
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
  case seissol::initializer::AllocationPreset::Timebucket:
    [[fallthrough]];
  case seissol::initializer::AllocationPreset::Timedofs:
    return (convergenceOrder <= 7 ? AllocationMode::HostOnlyHBM : AllocationMode::HostOnly);
  case seissol::initializer::AllocationPreset::Constant:
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
  case seissol::initializer::AllocationPreset::Timebucket:
    return useMPIUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplitPinned;
  default:
    return useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit;
  }
#endif
}

struct LTS {
  Variable<real[tensor::Q::size()]> dofs;
  // size is zero if Qane is not defined
  Variable<real[ZeroLengthArrayHandler(kernels::size<tensor::Qane>())]> dofsAne;
  Variable<real*> buffers;
  Variable<real*> derivatives;
  Variable<CellLocalInformation> cellInformation;
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
    tree.addVar(faceNeighbors,
                LayerMask(Ghost),
                1,
                allocationModeWP(AllocationPreset::TimedofsConstant),
                true);
    tree.addVar(
        localIntegration, LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.addVar(neighboringIntegration,
                LayerMask(Ghost),
                1,
                allocationModeWP(AllocationPreset::Constant),
                true);
    tree.addVar(material, LayerMask(Ghost), 1, AllocationMode::HostOnly, true);
    tree.addVar(
        plasticity, plasticityMask, 1, allocationModeWP(AllocationPreset::Plasticity), true);
    tree.addVar(drMapping, LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.addVar(
        boundaryMapping, LayerMask(Ghost), 1, allocationModeWP(AllocationPreset::Constant), true);
    tree.addVar(
        pstrain, plasticityMask, PagesizeHeap, allocationModeWP(AllocationPreset::Plasticity));
    tree.addVar(faceDisplacements, LayerMask(Ghost), PagesizeHeap, AllocationMode::HostOnly, true);

    tree.addBucket(
        buffersDerivatives, PagesizeHeap, allocationModeWP(AllocationPreset::Timebucket));
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
};

} // namespace seissol::initializer
#endif
