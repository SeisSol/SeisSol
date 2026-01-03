// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Constants.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/BatchRecorders/DataTypes/Condition.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalKey.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Initializer/BatchRecorders/DataTypes/EncodedConstants.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Recorders.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cassert>
#include <cstddef>
#include <vector>
#include <yateto.h>

using namespace device;
using namespace seissol::initializer;
using namespace seissol::recording;

void LocalIntegrationRecorder::record(LTS::Layer& layer) {
  setUpContext(layer);
  idofsAddressRegistry_.clear();

  recordTimeAndVolumeIntegrals();
  recordFreeSurfaceGravityBc();
  recordDirichletBc();
  recordAnalyticalBc(layer);
  recordLocalFluxIntegral();
  recordDisplacements();
}

void LocalIntegrationRecorder::recordTimeAndVolumeIntegrals() {
  real* integratedDofsScratch =
      static_cast<real*>(currentLayer_->var<LTS::IntegratedDofsScratch>(AllocationPlace::Device));
  real* derivativesScratch =
      static_cast<real*>(currentLayer_->var<LTS::DerivativesScratch>(AllocationPlace::Device));

  const auto size = currentLayer_->size();
  if (size > 0) {
    std::vector<real*> dofsPtrs(size, nullptr);
    std::vector<real*> dofsAnePtrs(size, nullptr);
    std::vector<real*> dofsExtPtrs(size, nullptr);
    std::vector<real*> localPtrs(size, nullptr);
    std::vector<real*> idofsPtrs;
    std::vector<real*> idofsAnePtrs(size, nullptr);
    std::vector<real*> derivativesAnePtrs(size, nullptr);
    std::vector<real*> derivativesExtPtrs(size, nullptr);
    std::vector<real*> ltsBuffers;
    std::vector<real*> idofsForLtsBuffers;

    idofsPtrs.reserve(size);
    dQPtrs_.resize(size);

    real** derivatives = currentLayer_->var<LTS::DerivativesDevice>();
    real** buffers = currentLayer_->var<LTS::BuffersDevice>();

    for (unsigned cell = 0; cell < size; ++cell) {
      auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);
      auto dataHost = currentLayer_->cellRef(cell, AllocationPlace::Host);

      // dofs
      dofsPtrs[cell] = static_cast<real*>(data.get<LTS::Dofs>());

      // idofs
      real* nextIdofPtr = &integratedDofsScratch[integratedDofsAddressCounter_];
      const bool isBuffersProvided = dataHost.get<LTS::CellInformation>().ltsSetup.hasBuffers();
      const bool isLtsBuffers = dataHost.get<LTS::CellInformation>().ltsSetup.accumulateBuffers();

      if (isBuffersProvided) {
        if (isLtsBuffers) {
          // lts buffers may require either accumulation or overriding (in case of reset command)
          idofsPtrs.push_back(nextIdofPtr);

          idofsForLtsBuffers.push_back(nextIdofPtr);
          ltsBuffers.push_back(buffers[cell]);

          idofsAddressRegistry_[cell] = nextIdofPtr;
          integratedDofsAddressCounter_ += tensor::I::size();
        } else {
          // gts buffers have to be always overridden
          idofsPtrs.push_back(buffers[cell]);
          idofsAddressRegistry_[cell] = buffers[cell];
        }
      } else {
        idofsPtrs.push_back(nextIdofPtr);
        idofsAddressRegistry_[cell] = nextIdofPtr;
        integratedDofsAddressCounter_ += tensor::I::size();
      }

      // stars
      localPtrs[cell] = reinterpret_cast<real*>(&data.get<LTS::LocalIntegration>());
#ifdef USE_VISCOELASTIC2
      auto* dofsAne = currentLayer_->var<LTS::DofsAne>(AllocationPlace::Device);
      dofsAnePtrs[cell] = dofsAne[cell];

      auto* idofsAne = currentLayer_->var<LTS::IDofsAneScratch>(AllocationPlace::Device);
      idofsAnePtrs[cell] = static_cast<real*>(idofsAne) + tensor::Iane::size() * cell;

      auto* derivativesExt =
          currentLayer_->var<LTS::DerivativesExtScratch>(AllocationPlace::Device);
      derivativesExtPtrs[cell] = static_cast<real*>(derivativesExt) +
                                 (tensor::dQext::size(1) + tensor::dQext::size(2)) * cell;

      auto* derivativesAne =
          currentLayer_->var<LTS::DerivativesAneScratch>(AllocationPlace::Device);
      derivativesAnePtrs[cell] = static_cast<real*>(derivativesAne) +
                                 (tensor::dQane::size(1) + tensor::dQane::size(2)) * cell;

      auto* dofsExt = currentLayer_->var<LTS::DofsExtScratch>(AllocationPlace::Device);
      dofsExtPtrs[cell] = static_cast<real*>(dofsExt) + tensor::Qext::size() * cell;
#endif

      // derivatives
      const bool isDerivativesProvided =
          dataHost.get<LTS::CellInformation>().ltsSetup.hasDerivatives();
      if (isDerivativesProvided) {
        dQPtrs_[cell] = derivatives[cell];

      } else {
        dQPtrs_[cell] = &derivativesScratch[derivativesAddressCounter_];
        derivativesAddressCounter_ += seissol::kernels::Solver::DerivativesSize;
      }
    }
    // just to be sure that we took all branches while filling in idofsPtrs vector
    assert(dofsPtrs.size() == idofsPtrs.size());

    const ConditionalKey key(KernelNames::Time || KernelNames::Volume);
    checkKey(key);

    (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::LocalIntegrationData, localPtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, idofsPtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::Derivatives, dQPtrs_);

    if (!idofsForLtsBuffers.empty()) {
      const ConditionalKey key(*KernelNames::Time, *ComputationKind::WithLtsBuffers);

      (*currentTable_)[key].set(inner_keys::Wp::Id::Buffers, ltsBuffers);
      (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, idofsForLtsBuffers);
    }

#ifdef USE_VISCOELASTIC2
    (*currentTable_)[key].set(inner_keys::Wp::Id::DofsAne, dofsAnePtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::DofsExt, dofsExtPtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::IdofsAne, idofsAnePtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::DerivativesAne, derivativesAnePtrs);
    (*currentTable_)[key].set(inner_keys::Wp::Id::DerivativesExt, derivativesExtPtrs);
#endif
  }
}

void LocalIntegrationRecorder::recordLocalFluxIntegral() {
  const auto size = currentLayer_->size();
  for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
    std::vector<real*> idofsPtrs;
    std::vector<real*> dofsPtrs;
    std::vector<real*> localPtrs;

    std::vector<real*> dofsExtPtrs;

    idofsPtrs.reserve(size);
    dofsPtrs.reserve(size);
    localPtrs.reserve(size);

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);
      auto dataHost = currentLayer_->cellRef(cell, AllocationPlace::Host);

      // no element local contribution in the case of dynamic rupture boundary conditions
      if (dataHost.get<LTS::CellInformation>().faceTypes[face] != FaceType::DynamicRupture) {
        idofsPtrs.push_back(idofsAddressRegistry_[cell]);
        dofsPtrs.push_back(static_cast<real*>(data.get<LTS::Dofs>()));
        localPtrs.push_back(reinterpret_cast<real*>(&data.get<LTS::LocalIntegration>()));
#ifdef USE_VISCOELASTIC2
        auto* dofsExt = currentLayer_->var<LTS::DofsExtScratch>(AllocationPlace::Device);
        dofsExtPtrs.push_back(static_cast<real*>(dofsExt) + tensor::Qext::size() * cell);
#endif
      }
    }

    // NOTE: we can check any container, but we must check that a set is not empty!
    if (!dofsPtrs.empty()) {
      const ConditionalKey key(*KernelNames::LocalFlux, !FaceKinds::DynamicRupture, face);
      checkKey(key);
      (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, idofsPtrs);
      (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs);
      (*currentTable_)[key].set(inner_keys::Wp::Id::LocalIntegrationData, localPtrs);
#ifdef USE_VISCOELASTIC2
      (*currentTable_)[key].set(inner_keys::Wp::Id::DofsExt, dofsExtPtrs);
#endif
    }
  }
}

void LocalIntegrationRecorder::recordDisplacements() {
  real*(*faceDisplacements)[4] = currentLayer_->var<LTS::FaceDisplacementsDevice>();
  std::array<std::vector<real*>, 4> iVelocitiesPtrs{{}};
  std::array<std::vector<real*>, 4> displacementsPtrs{};

  const auto size = currentLayer_->size();
  for (std::size_t cell = 0; cell < size; ++cell) {
    auto dataHost = currentLayer_->cellRef(cell, AllocationPlace::Host);

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      auto isRequired = faceDisplacements[cell][face] != nullptr;
      auto notFreeSurfaceGravity =
          dataHost.get<LTS::CellInformation>().faceTypes[face] != FaceType::FreeSurfaceGravity;

      if (isRequired && notFreeSurfaceGravity) {
        auto iview = init::I::view::create(idofsAddressRegistry_[cell]);
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
      (*currentTable_)[key].set(inner_keys::Wp::Id::Ivelocities, iVelocitiesPtrs[face]);
      (*currentTable_)[key].set(inner_keys::Wp::Id::FaceDisplacement, displacementsPtrs[face]);
    }
  }
}

void LocalIntegrationRecorder::recordFreeSurfaceGravityBc() {
  const auto size = currentLayer_->size();
  constexpr size_t NodalAvgDisplacementsSize = tensor::averageNormalDisplacement::size();

  real* nodalAvgDisplacements =
      static_cast<real*>(currentLayer_->var<LTS::NodalAvgDisplacements>(AllocationPlace::Device));

  real* rotateDisplacementToFaceNormalScratch = static_cast<real*>(
      currentLayer_->var<LTS::RotateDisplacementToFaceNormalScratch>(AllocationPlace::Device));
  real* rotateDisplacementToGlobalScratch = static_cast<real*>(
      currentLayer_->var<LTS::RotateDisplacementToGlobalScratch>(AllocationPlace::Device));
  real* rotatedFaceDisplacementScratch = static_cast<real*>(
      currentLayer_->var<LTS::RotatedFaceDisplacementScratch>(AllocationPlace::Device));
  real* dofsFaceNodalScratch =
      static_cast<real*>(currentLayer_->var<LTS::DofsFaceNodalScratch>(AllocationPlace::Device));
  real* prevCoefficientsScratch =
      static_cast<real*>(currentLayer_->var<LTS::PrevCoefficientsScratch>(AllocationPlace::Device));

  real* dofsFaceBoundaryNodalScratch = static_cast<real*>(
      currentLayer_->var<LTS::DofsFaceBoundaryNodalScratch>(AllocationPlace::Device));

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
      auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);
      auto dataHost = currentLayer_->cellRef(cell, AllocationPlace::Host);

      for (std::size_t face = 0; face < 4; ++face) {
        if (dataHost.get<LTS::CellInformation>().faceTypes[face] == FaceType::FreeSurfaceGravity) {
          assert(dataHost.get<LTS::FaceDisplacementsDevice>()[face] != nullptr);
          cellIndices[face].push_back(cell);

          derivatives[face].push_back(dQPtrs_[cell]);
          dofsPtrs[face].push_back(static_cast<real*>(data.get<LTS::Dofs>()));
          idofsPtrs[face].push_back(idofsAddressRegistry_[cell]);

          neighPtrs[face].push_back(
              reinterpret_cast<real*>(&data.get<LTS::NeighboringIntegration>()));
          displacementsPtrs[face].push_back(dataHost.get<LTS::FaceDisplacementsDevice>()[face]);
          t[face].push_back(dataHost.get<LTS::BoundaryMappingDevice>()[face].dataT);
          tInv[face].push_back(dataHost.get<LTS::BoundaryMappingDevice>()[face].dataTinv);

          rhos[face].push_back(dataHost.get<LTS::Material>().local->getDensity());
          lambdas[face].push_back(dataHost.get<LTS::Material>().local->getLambdaBar());

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
          prevCoefficientsPtrs[face].push_back(
              prevCoefficientsScratch +
              counter[face] * nodal::tensor::nodes2D::Shape[multisim::BasisFunctionDimension]);
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
        (*currentIndicesTable_)[key].set(inner_keys::Indices::Id::Cells, cellIndices[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::Derivatives, derivatives[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, idofsPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::NeighborIntegrationData, neighPtrs[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::T, t[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Tinv, tInv[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::FaceDisplacement, displacementsPtrs[face]);
        (*currentMaterialTable_)[key].set(inner_keys::Material::Id::Rho, rhos[face]);
        (*currentMaterialTable_)[key].set(inner_keys::Material::Id::Lambda, lambdas[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::NodalAvgDisplacements,
                                  nodalAvgDisplacementsPtrs[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::RotateDisplacementToFaceNormal,
                                  rotateDisplacementToFaceNormalPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::RotateDisplacementToGlobal,
                                  rotateDisplacementToGlobalPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::RotatedFaceDisplacement,
                                  rotatedFaceDisplacementPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::DofsFaceNodal, dofsFaceNodalPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::PrevCoefficients, prevCoefficientsPtrs[face]);
        (*currentMaterialTable_)[key].set(inner_keys::Material::Id::InvImpedances,
                                          invImpedances[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::DofsFaceBoundaryNodal,
                                  dofsFaceBoundaryNodalPtrs[face]);
      }
    }
  }
}

void LocalIntegrationRecorder::recordDirichletBc() {
  const auto size = currentLayer_->size();
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
        currentLayer_->var<LTS::DofsFaceBoundaryNodalScratch>(AllocationPlace::Device));

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);
      auto dataHost = currentLayer_->cellRef(cell, AllocationPlace::Host);

      for (std::size_t face = 0; face < 4; ++face) {
        if (dataHost.get<LTS::CellInformation>().faceTypes[face] == FaceType::Dirichlet) {

          dofsPtrs[face].push_back(static_cast<real*>(data.get<LTS::Dofs>()));
          idofsPtrs[face].push_back(idofsAddressRegistry_[cell]);

          tInv[face].push_back(dataHost.get<LTS::BoundaryMappingDevice>()[face].dataTinv);
          neighPtrs[face].push_back(
              reinterpret_cast<real*>(&data.get<LTS::NeighboringIntegration>()));

          easiBoundaryMapPtrs[face].push_back(
              dataHost.get<LTS::BoundaryMappingDevice>()[face].easiBoundaryMap);
          easiBoundaryConstantPtrs[face].push_back(
              dataHost.get<LTS::BoundaryMappingDevice>()[face].easiBoundaryConstant);

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
        (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Idofs, idofsPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::NeighborIntegrationData, neighPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::Tinv, tInv[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::EasiBoundaryMap, easiBoundaryMapPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::EasiBoundaryConstant,
                                  easiBoundaryConstantPtrs[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::DofsFaceBoundaryNodal,
                                  dofsFaceBoundaryNodalPtrs[face]);
      }
    }
  }
}

void LocalIntegrationRecorder::recordAnalyticalBc(LTS::Layer& layer) {
  const auto size = currentLayer_->size();
  if (size > 0) {
    std::array<std::vector<real*>, 4> dofsPtrs{};
    std::array<std::vector<real*>, 4> neighPtrs{};
    std::array<std::vector<unsigned>, 4> cellIndices{};
    std::array<std::vector<real*>, 4> analytical{};

    real* analyticScratch =
        reinterpret_cast<real*>(layer.var<LTS::AnalyticScratch>(AllocationPlace::Device));

    for (std::size_t cell = 0; cell < size; ++cell) {
      auto data = currentLayer_->cellRef(cell, AllocationPlace::Device);
      auto dataHost = currentLayer_->cellRef(cell, AllocationPlace::Host);

      for (std::size_t face = 0; face < 4; ++face) {
        if (dataHost.get<LTS::CellInformation>().faceTypes[face] == FaceType::Analytical) {
          cellIndices[face].push_back(cell);
          dofsPtrs[face].push_back(data.get<LTS::Dofs>());
          neighPtrs[face].push_back(
              reinterpret_cast<real*>(&data.get<LTS::NeighboringIntegration>()));
          analytical[face].push_back(analyticScratch + cell * tensor::INodal::size());
        }
      }
    }

    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      if (!cellIndices[face].empty()) {
        const ConditionalKey key(
            *KernelNames::BoundaryConditions, *ComputationKind::Analytical, face);
        checkKey(key);
        (*currentIndicesTable_)[key].set(inner_keys::Indices::Id::Cells, cellIndices[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::Dofs, dofsPtrs[face]);
        (*currentTable_)[key].set(inner_keys::Wp::Id::NeighborIntegrationData, neighPtrs[face]);

        (*currentTable_)[key].set(inner_keys::Wp::Id::Analytical, analytical[face]);
      }
    }
  }
}
