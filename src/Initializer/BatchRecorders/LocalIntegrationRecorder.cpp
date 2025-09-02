// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Interface.h"
#include "Recorders.h"
#include <Common/Constants.h>
#include <DataTypes/ConditionalKey.h>
#include <Initializer/BasicTypedefs.h>
#include <Kernels/Precision.h>
#include <Kernels/Solver.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <array>
#include <cassert>
#include <cstddef>
#include <vector>
#include <yateto.h>

#include "DataTypes/Condition.h"
#include "DataTypes/ConditionalTable.h"
#include "DataTypes/EncodedConstants.h"

using namespace device;
using namespace seissol::initializer;
using namespace seissol::initializer::recording;

void LocalIntegrationRecorder::record(LTS& handler, Layer& layer) {
  kernels::LocalData::Loader loader, loaderHost;
  loader.load(handler, layer, AllocationPlace::Device);
  loaderHost.load(handler, layer, AllocationPlace::Host);
  setUpContext(handler, layer, loader, loaderHost);
  idofsAddressRegistry.clear();

  recordTimeAndVolumeIntegrals();
  recordFreeSurfaceGravityBc();
  recordDirichletBc();
  recordAnalyticalBc(handler, layer);
  recordLocalFluxIntegral();
  recordDisplacements();
}

void LocalIntegrationRecorder::recordTimeAndVolumeIntegrals() {
  real* integratedDofsScratch = static_cast<real*>(
      currentLayer->var(currentHandler->integratedDofsScratch, AllocationPlace::Device));
  real* derivativesScratch = static_cast<real*>(
      currentLayer->var(currentHandler->derivativesScratch, AllocationPlace::Device));

  const auto size = currentLayer->size();
  if (size > 0) {
    std::vector<real*> dofsPtrs(size, nullptr);
    std::vector<real*> dofsAnePtrs(size, nullptr);
    std::vector<real*> dofsExtPtrs(size, nullptr);
    std::vector<real*> localPtrs(size, nullptr);
    std::vector<real*> idofsPtrs{};
    std::vector<real*> idofsAnePtrs(size, nullptr);
    std::vector<real*> derivativesAnePtrs(size, nullptr);
    std::vector<real*> derivativesExtPtrs(size, nullptr);
    std::vector<real*> ltsBuffers{};
    std::vector<real*> idofsForLtsBuffers{};

    idofsPtrs.reserve(size);
    dQPtrs.resize(size);

    real** derivatives = currentLayer->var(currentHandler->derivativesDevice);
    real** buffers = currentLayer->var(currentHandler->buffersDevice);

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);
      auto dataHost = currentLoaderHost->entry(cell);

      // dofs
      dofsPtrs[cell] = static_cast<real*>(data.dofs());

      // idofs
      real* nextIdofPtr = &integratedDofsScratch[integratedDofsAddressCounter];
      const bool isBuffersProvided = ((dataHost.cellInformation().ltsSetup >> 8) % 2) == 1;
      const bool isLtsBuffers = ((dataHost.cellInformation().ltsSetup >> 10) % 2) == 1;

      if (isBuffersProvided) {
        if (isLtsBuffers) {
          // lts buffers may require either accumulation or overriding (in case of reset command)
          idofsPtrs.push_back(nextIdofPtr);

          idofsForLtsBuffers.push_back(nextIdofPtr);
          ltsBuffers.push_back(buffers[cell]);

          idofsAddressRegistry[cell] = nextIdofPtr;
          integratedDofsAddressCounter += tensor::I::size();
        } else {
          // gts buffers have to be always overridden
          idofsPtrs.push_back(buffers[cell]);
          idofsAddressRegistry[cell] = buffers[cell];
        }
      } else {
        idofsPtrs.push_back(nextIdofPtr);
        idofsAddressRegistry[cell] = nextIdofPtr;
        integratedDofsAddressCounter += tensor::I::size();
      }

      // stars
      localPtrs[cell] = reinterpret_cast<real*>(&data.localIntegration());
#ifdef USE_VISCOELASTIC2
      auto* dofsAne = currentLayer->var(currentHandler->dofsAne, AllocationPlace::Device);
      dofsAnePtrs[cell] = dofsAne[cell];

      auto* idofsAne = currentLayer->var(currentHandler->idofsAneScratch, AllocationPlace::Device);
      idofsAnePtrs[cell] = static_cast<real*>(idofsAne) + tensor::Iane::size() * cell;

      auto* derivativesExt =
          currentLayer->var(currentHandler->derivativesExtScratch, AllocationPlace::Device);
      derivativesExtPtrs[cell] = static_cast<real*>(derivativesExt) +
                                 (tensor::dQext::size(1) + tensor::dQext::size(2)) * cell;

      auto* derivativesAne =
          currentLayer->var(currentHandler->derivativesAneScratch, AllocationPlace::Device);
      derivativesAnePtrs[cell] = static_cast<real*>(derivativesAne) +
                                 (tensor::dQane::size(1) + tensor::dQane::size(2)) * cell;

      auto* dofsExt = currentLayer->var(currentHandler->dofsExtScratch, AllocationPlace::Device);
      dofsExtPtrs[cell] = static_cast<real*>(dofsExt) + tensor::Qext::size() * cell;
#endif

      // derivatives
      const bool isDerivativesProvided = ((dataHost.cellInformation().ltsSetup >> 9) % 2) == 1;
      if (isDerivativesProvided) {
        dQPtrs[cell] = derivatives[cell];

      } else {
        dQPtrs[cell] = &derivativesScratch[derivativesAddressCounter];
        derivativesAddressCounter += seissol::kernels::Solver::DerivativesSize;
      }
    }
    // just to be sure that we took all branches while filling in idofsPtrs vector
    assert(dofsPtrs.size() == idofsPtrs.size());

    const ConditionalKey key(KernelNames::Time || KernelNames::Volume);
    checkKey(key);

    (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::LocalIntegrationData, localPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, idofsPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::Derivatives, dQPtrs);

    if (!idofsForLtsBuffers.empty()) {
      const ConditionalKey key(*KernelNames::Time, *ComputationKind::WithLtsBuffers);

      (*currentTable)[key].set(inner_keys::Wp::Id::Buffers, ltsBuffers);
      (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, idofsForLtsBuffers);
    }

#ifdef USE_VISCOELASTIC2
    (*currentTable)[key].set(inner_keys::Wp::Id::DofsAne, dofsAnePtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::DofsExt, dofsExtPtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::IdofsAne, idofsAnePtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::DerivativesAne, derivativesAnePtrs);
    (*currentTable)[key].set(inner_keys::Wp::Id::DerivativesExt, derivativesExtPtrs);
#endif
  }
}

void LocalIntegrationRecorder::recordLocalFluxIntegral() {
  const auto size = currentLayer->size();
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    std::vector<real*> idofsPtrs{};
    std::vector<real*> dofsPtrs{};
    std::vector<real*> localPtrs{};

    std::vector<real*> dofsExtPtrs{};

    idofsPtrs.reserve(size);
    dofsPtrs.reserve(size);
    localPtrs.reserve(size);

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);
      auto dataHost = currentLoaderHost->entry(cell);

      // no element local contribution in the case of dynamic rupture boundary conditions
      if (dataHost.cellInformation().faceTypes[face] != FaceType::DynamicRupture) {
        idofsPtrs.push_back(idofsAddressRegistry[cell]);
        dofsPtrs.push_back(static_cast<real*>(data.dofs()));
        localPtrs.push_back(reinterpret_cast<real*>(&data.localIntegration()));
#ifdef USE_VISCOELASTIC2
        auto* dofsExt = currentLayer->var(currentHandler->dofsExtScratch, AllocationPlace::Device);
        dofsExtPtrs.push_back(static_cast<real*>(dofsExt) + tensor::Qext::size() * cell);
#endif
      }
    }

    // NOTE: we can check any container, but we must check that a set is not empty!
    if (!dofsPtrs.empty()) {
      const ConditionalKey key(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);
      checkKey(key);
      (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, idofsPtrs);
      (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
      (*currentTable)[key].set(inner_keys::Wp::Id::LocalIntegrationData, localPtrs);
#ifdef USE_VISCOELASTIC2
      (*currentTable)[key].set(inner_keys::Wp::Id::DofsExt, dofsExtPtrs);
#endif
    }
  }
}

void LocalIntegrationRecorder::recordDisplacements() {
  real*(*faceDisplacements)[4] = currentLayer->var(currentHandler->faceDisplacementsDevice);
  std::array<std::vector<real*>, 4> iVelocitiesPtrs{{}};
  std::array<std::vector<real*>, 4> displacementsPtrs{};

  const auto size = currentLayer->size();
  for (std::size_t cell = 0; cell < size; ++cell) {
    auto dataHost = currentLoaderHost->entry(cell);

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      auto isRequired = faceDisplacements[cell][face] != nullptr;
      auto notFreeSurfaceGravity =
          dataHost.cellInformation().faceTypes[face] != FaceType::FreeSurfaceGravity;

      if (isRequired && notFreeSurfaceGravity) {
        auto iview = init::I::view::create(idofsAddressRegistry[cell]);
        // NOTE: velocity components are between 6th and 8th columns
        constexpr unsigned FirstVelocityComponent{6};
        iVelocitiesPtrs[face].push_back(
            &multisim::multisimWrap(iview, 0, 0, FirstVelocityComponent));
        displacementsPtrs[face].push_back(faceDisplacements[cell][face]);
      }
    }
  }

  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    if (!displacementsPtrs[face].empty()) {
      const ConditionalKey key(*KernelNames::FaceDisplacements, *ComputationKind::None, face);
      checkKey(key);
      (*currentTable)[key].set(inner_keys::Wp::Id::Ivelocities, iVelocitiesPtrs[face]);
      (*currentTable)[key].set(inner_keys::Wp::Id::FaceDisplacement, displacementsPtrs[face]);
    }
  }
}

void LocalIntegrationRecorder::recordFreeSurfaceGravityBc() {
  const auto size = currentLayer->size();
  constexpr size_t NodalAvgDisplacementsSize = tensor::averageNormalDisplacement::size();

  real* nodalAvgDisplacements = static_cast<real*>(
      currentLayer->var(currentHandler->nodalAvgDisplacements, AllocationPlace::Device));

  real* rotateDisplacementToFaceNormalScratch = static_cast<real*>(currentLayer->var(
      currentHandler->rotateDisplacementToFaceNormalScratch, AllocationPlace::Device));
  real* rotateDisplacementToGlobalScratch = static_cast<real*>(currentLayer->var(
      currentHandler->rotateDisplacementToGlobalScratch, AllocationPlace::Device));
  real* rotatedFaceDisplacementScratch = static_cast<real*>(
      currentLayer->var(currentHandler->rotatedFaceDisplacementScratch, AllocationPlace::Device));
  real* dofsFaceNodalScratch = static_cast<real*>(
      currentLayer->var(currentHandler->dofsFaceNodalScratch, AllocationPlace::Device));
  real* prevCoefficientsScratch = static_cast<real*>(
      currentLayer->var(currentHandler->prevCoefficientsScratch, AllocationPlace::Device));

  real* dofsFaceBoundaryNodalScratch = static_cast<real*>(
      currentLayer->var(currentHandler->dofsFaceBoundaryNodalScratch, AllocationPlace::Device));

  if (size > 0) {
    std::array<std::vector<unsigned>, 4> cellIndices{};
    std::array<std::vector<real*>, 4> nodalAvgDisplacementsPtrs{};
    std::array<std::vector<real*>, 4> displacementsPtrs{};

    std::array<std::vector<real*>, 4> derivatives{};
    std::array<std::vector<real*>, 4> dofsPtrs{};
    std::array<std::vector<real*>, 4> idofsPtrs{};
    std::array<std::vector<real*>, 4> neighPtrs{};
    std::array<std::vector<real*>, 4> t{};
    std::array<std::vector<real*>, 4> tInv{};

    std::array<std::vector<inner_keys::Material::DataType>, 4> rhos;
    std::array<std::vector<inner_keys::Material::DataType>, 4> lambdas;

    std::array<std::vector<real*>, 4> rotateDisplacementToFaceNormalPtrs{};
    std::array<std::vector<real*>, 4> rotateDisplacementToGlobalPtrs{};
    std::array<std::vector<real*>, 4> rotatedFaceDisplacementPtrs{};
    std::array<std::vector<real*>, 4> dofsFaceNodalPtrs{};
    std::array<std::vector<real*>, 4> prevCoefficientsPtrs{};
    std::array<std::vector<double>, 4> invImpedances{};
    std::array<std::vector<real*>, 4> dofsFaceBoundaryNodalPtrs{};

    std::array<std::size_t, 4> counter{};

    size_t nodalAvgDisplacementsCounter{0};

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);
      auto dataHost = currentLoaderHost->entry(cell);

      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (dataHost.cellInformation().faceTypes[face] == FaceType::FreeSurfaceGravity) {
          assert(dataHost.faceDisplacementsDevice()[face] != nullptr);
          cellIndices[face].push_back(cell);

          derivatives[face].push_back(dQPtrs[cell]);
          dofsPtrs[face].push_back(static_cast<real*>(data.dofs()));
          idofsPtrs[face].push_back(idofsAddressRegistry[cell]);

          neighPtrs[face].push_back(reinterpret_cast<real*>(&data.neighboringIntegration()));
          displacementsPtrs[face].push_back(dataHost.faceDisplacementsDevice()[face]);
          t[face].push_back(dataHost.boundaryMappingDevice()[face].dataT);
          tInv[face].push_back(dataHost.boundaryMappingDevice()[face].dataTinv);

          rhos[face].push_back(dataHost.material().local->getDensity());
          lambdas[face].push_back(dataHost.material().local->getLambdaBar());

          real* displ{&nodalAvgDisplacements[nodalAvgDisplacementsCounter]};
          nodalAvgDisplacementsPtrs[face].push_back(displ);
          nodalAvgDisplacementsCounter += NodalAvgDisplacementsSize;

          rotateDisplacementToFaceNormalPtrs[face].push_back(
              rotateDisplacementToFaceNormalScratch +
              counter[face] * init::displacementRotationMatrix::Size);
          rotateDisplacementToGlobalPtrs[face].push_back(
              rotateDisplacementToGlobalScratch +
              counter[face] * init::displacementRotationMatrix::Size);
          rotatedFaceDisplacementPtrs[face].push_back(
              rotatedFaceDisplacementScratch + counter[face] * init::rotatedFaceDisplacement::Size);
          dofsFaceBoundaryNodalPtrs[face].push_back(dofsFaceBoundaryNodalScratch +
                                                    counter[face] * tensor::INodal::size());
          dofsFaceNodalPtrs[face].push_back(dofsFaceNodalScratch +
                                            counter[face] * tensor::INodal::size());
          prevCoefficientsPtrs[face].push_back(prevCoefficientsScratch +
                                               counter[face] * nodal::tensor::nodes2D::Shape[0]);
          invImpedances[face].push_back(0);

          ++counter[face];
        }
      }
    }

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (!cellIndices[face].empty()) {
        const ConditionalKey key(
            *KernelNames::BoundaryConditions, *ComputationKind::FreeSurfaceGravity, face);
        checkKey(key);
        (*currentIndicesTable)[key].set(inner_keys::Indices::Id::Cells, cellIndices[face]);

        (*currentTable)[key].set(inner_keys::Wp::Id::Derivatives, derivatives[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, idofsPtrs[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::NeighborIntegrationData, neighPtrs[face]);

        (*currentTable)[key].set(inner_keys::Wp::Id::T, t[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::Tinv, tInv[face]);

        (*currentTable)[key].set(inner_keys::Wp::Id::FaceDisplacement, displacementsPtrs[face]);
        (*currentMaterialTable)[key].set(inner_keys::Material::Id::Rho, rhos[face]);
        (*currentMaterialTable)[key].set(inner_keys::Material::Id::Lambda, lambdas[face]);

        (*currentTable)[key].set(inner_keys::Wp::Id::NodalAvgDisplacements,
                                 nodalAvgDisplacementsPtrs[face]);

        (*currentTable)[key].set(inner_keys::Wp::Id::RotateDisplacementToFaceNormal,
                                 rotateDisplacementToFaceNormalPtrs[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::RotateDisplacementToGlobal,
                                 rotateDisplacementToGlobalPtrs[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::RotatedFaceDisplacement,
                                 rotatedFaceDisplacementPtrs[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::DofsFaceNodal, dofsFaceNodalPtrs[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::PrevCoefficients, prevCoefficientsPtrs[face]);
        (*currentMaterialTable)[key].set(inner_keys::Material::Id::InvImpedances,
                                         invImpedances[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::DofsFaceBoundaryNodal,
                                 dofsFaceBoundaryNodalPtrs[face]);
      }
    }
  }
}

void LocalIntegrationRecorder::recordDirichletBc() {
  const auto size = currentLayer->size();
  if (size > 0) {
    std::array<std::vector<real*>, 4> dofsPtrs{};
    std::array<std::vector<real*>, 4> idofsPtrs{};
    std::array<std::vector<real*>, 4> tInv{};
    std::array<std::vector<real*>, 4> neighPtrs{};

    std::array<std::vector<real*>, 4> easiBoundaryMapPtrs{};
    std::array<std::vector<real*>, 4> easiBoundaryConstantPtrs{};

    std::array<std::vector<real*>, 4> dofsFaceBoundaryNodalPtrs{};

    std::array<std::size_t, 4> counter{};

    real* dofsFaceBoundaryNodalScratch = static_cast<real*>(
        currentLayer->var(currentHandler->dofsFaceBoundaryNodalScratch, AllocationPlace::Device));

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto data = currentLoader->entry(cell);
      auto dataHost = currentLoaderHost->entry(cell);

      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (dataHost.cellInformation().faceTypes[face] == FaceType::Dirichlet) {

          dofsPtrs[face].push_back(static_cast<real*>(data.dofs()));
          idofsPtrs[face].push_back(idofsAddressRegistry[cell]);

          tInv[face].push_back(dataHost.boundaryMappingDevice()[face].dataTinv);
          neighPtrs[face].push_back(reinterpret_cast<real*>(&data.neighboringIntegration()));

          easiBoundaryMapPtrs[face].push_back(
              dataHost.boundaryMappingDevice()[face].easiBoundaryMap);
          easiBoundaryConstantPtrs[face].push_back(
              dataHost.boundaryMappingDevice()[face].easiBoundaryConstant);

          dofsFaceBoundaryNodalPtrs[face].push_back(dofsFaceBoundaryNodalScratch +
                                                    counter[face] * tensor::INodal::size());
          ++counter[face];
        }
      }
    }

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (!dofsPtrs[face].empty()) {
        const ConditionalKey key(
            *KernelNames::BoundaryConditions, *ComputationKind::Dirichlet, face);
        checkKey(key);
        (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::Idofs, idofsPtrs[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::NeighborIntegrationData, neighPtrs[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::Tinv, tInv[face]);

        (*currentTable)[key].set(inner_keys::Wp::Id::EasiBoundaryMap, easiBoundaryMapPtrs[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::EasiBoundaryConstant,
                                 easiBoundaryConstantPtrs[face]);

        (*currentTable)[key].set(inner_keys::Wp::Id::DofsFaceBoundaryNodal,
                                 dofsFaceBoundaryNodalPtrs[face]);
      }
    }
  }
}

void LocalIntegrationRecorder::recordAnalyticalBc(LTS& handler, Layer& layer) {
  const auto size = currentLayer->size();
  if (size > 0) {
    std::array<std::vector<real*>, 4> dofsPtrs{};
    std::array<std::vector<real*>, 4> neighPtrs{};
    std::array<std::vector<unsigned>, 4> cellIndices{};
    std::array<std::vector<real*>, 4> analytical{};

    real* analyticScratch =
        reinterpret_cast<real*>(layer.var(handler.analyticScratch, AllocationPlace::Device));

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto dataHost = currentLoaderHost->entry(cell);
      auto data = currentLoader->entry(cell);

      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        if (dataHost.cellInformation().faceTypes[face] == FaceType::Analytical) {
          cellIndices[face].push_back(cell);
          dofsPtrs[face].push_back(data.dofs());
          neighPtrs[face].push_back(reinterpret_cast<real*>(&data.neighboringIntegration()));
          analytical[face].push_back(analyticScratch + cell * tensor::INodal::size());
        }
      }
    }

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (!cellIndices[face].empty()) {
        const ConditionalKey key(
            *KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
        checkKey(key);
        (*currentIndicesTable)[key].set(inner_keys::Indices::Id::Cells, cellIndices[face]);

        (*currentTable)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs[face]);
        (*currentTable)[key].set(inner_keys::Wp::Id::NeighborIntegrationData, neighPtrs[face]);

        (*currentTable)[key].set(inner_keys::Wp::Id::Analytical, analytical[face]);
      }
    }
  }
}
