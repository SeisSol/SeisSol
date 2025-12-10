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
#include "Kernels/Interface.h"
#include "Kernels/Precision.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Recorders.h"

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
  idofsAddressRegistry.clear();

  recordTimeAndVolumeIntegrals();
  recordFreeSurfaceGravityBc();
  recordDirichletBc();
  recordAnalyticalBc(layer);
  recordLocalFluxIntegral();
  recordDisplacements();
}

void LocalIntegrationRecorder::recordTimeAndVolumeIntegrals() {
  currentLayer->wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;

    real* integratedDofsScratch = static_cast<real*>(
        currentLayer->var<LTS::IntegratedDofsScratch>(cfg, AllocationPlace::Device));
    real* derivativesScratch = static_cast<real*>(
        currentLayer->var<LTS::DerivativesScratch>(cfg, AllocationPlace::Device));

    const auto size = currentLayer->size();
    if (size > 0) {
      std::vector<void*> dofsPtrs(size, nullptr);
      std::vector<void*> dofsAnePtrs(size, nullptr);
      std::vector<void*> dofsExtPtrs(size, nullptr);
      std::vector<void*> localPtrs(size, nullptr);
      std::vector<void*> idofsPtrs{};
      std::vector<void*> idofsAnePtrs(size, nullptr);
      std::vector<void*> derivativesAnePtrs(size, nullptr);
      std::vector<void*> derivativesExtPtrs(size, nullptr);
      std::vector<void*> ltsBuffers{};
      std::vector<void*> idofsForLtsBuffers{};

      idofsPtrs.reserve(size);
      dQPtrs.resize(size);

      real** derivatives = currentLayer->var<LTS::DerivativesDevice>(cfg);
      real** buffers = currentLayer->var<LTS::BuffersDevice>(cfg);

      for (unsigned cell = 0; cell < size; ++cell) {
        auto data = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Device);
        auto dataHost = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Host);

        // dofs
        dofsPtrs[cell] = static_cast<real*>(data.template get<LTS::Dofs>());

        // idofs
        real* nextIdofPtr = &integratedDofsScratch[integratedDofsAddressCounter];
        const bool isBuffersProvided =
            dataHost.template get<LTS::CellInformation>().ltsSetup.hasBuffers();
        const bool isLtsBuffers =
            dataHost.template get<LTS::CellInformation>().ltsSetup.accumulateBuffers();

        if (isBuffersProvided) {
          if (isLtsBuffers) {
            // lts buffers may require either accumulation or overriding (in case of reset command)
            idofsPtrs.push_back(nextIdofPtr);

            idofsForLtsBuffers.push_back(nextIdofPtr);
            ltsBuffers.push_back(buffers[cell]);

            idofsAddressRegistry[cell] = nextIdofPtr;
            integratedDofsAddressCounter += tensor::I<Cfg>::size();
          } else {
            // gts buffers have to be always overridden
            idofsPtrs.push_back(buffers[cell]);
            idofsAddressRegistry[cell] = buffers[cell];
          }
        } else {
          idofsPtrs.push_back(nextIdofPtr);
          idofsAddressRegistry[cell] = nextIdofPtr;
          integratedDofsAddressCounter += tensor::I<Cfg>::size();
        }

        // stars
        localPtrs[cell] = reinterpret_cast<real*>(&data.template get<LTS::LocalIntegration>());
#ifdef USE_VISCOELASTIC2
        auto* dofsAne = currentLayer->var<LTS::DofsAne>(cfg, AllocationPlace::Device);
        dofsAnePtrs[cell] = dofsAne[cell];

        auto* idofsAne = currentLayer->var<LTS::IDofsAneScratch>(cfg, AllocationPlace::Device);
        idofsAnePtrs[cell] = static_cast<real*>(idofsAne) + tensor::Iane<Cfg>::size() * cell;

        auto* derivativesExt =
            currentLayer->var<LTS::DerivativesExtScratch>(cfg, AllocationPlace::Device);
        derivativesExtPtrs[cell] =
            static_cast<real*>(derivativesExt) +
            (tensor::dQext<Cfg>::size(1) + tensor::dQext<Cfg>::size(2)) * cell;

        auto* derivativesAne =
            currentLayer->var<LTS::DerivativesAneScratch>(cfg, AllocationPlace::Device);
        derivativesAnePtrs[cell] =
            static_cast<real*>(derivativesAne) +
            (tensor::dQane<Cfg>::size(1) + tensor::dQane<Cfg>::size(2)) * cell;

        auto* dofsExt = currentLayer->var<LTS::DofsExtScratch>(cfg, AllocationPlace::Device);
        dofsExtPtrs[cell] = static_cast<real*>(dofsExt) + tensor::Qext<Cfg>::size() * cell;
#endif

        // derivatives
        const bool isDerivativesProvided =
            dataHost.template get<LTS::CellInformation>().ltsSetup.hasDerivatives();
        if (isDerivativesProvided) {
          dQPtrs[cell] = derivatives[cell];

        } else {
          dQPtrs[cell] = &derivativesScratch[derivativesAddressCounter];
          derivativesAddressCounter += seissol::kernels::Solver<Cfg>::template DerivativesSize<Cfg>;
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
  });
}

void LocalIntegrationRecorder::recordLocalFluxIntegral() {
  currentLayer->wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;

    const auto size = currentLayer->size();
    for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
      std::vector<void*> idofsPtrs{};
      std::vector<void*> dofsPtrs{};
      std::vector<void*> localPtrs{};

      std::vector<void*> dofsExtPtrs{};

      idofsPtrs.reserve(size);
      dofsPtrs.reserve(size);
      localPtrs.reserve(size);

      for (std::size_t cell = 0; cell < size; ++cell) {
        auto data = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Device);
        auto dataHost = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Host);

        // no element local contribution in the case of dynamic rupture boundary conditions
        if (dataHost.template get<LTS::CellInformation>().faceTypes[face] !=
            FaceType::DynamicRupture) {
          idofsPtrs.push_back(idofsAddressRegistry[cell]);
          dofsPtrs.push_back(static_cast<real*>(data.template get<LTS::Dofs>()));
          localPtrs.push_back(reinterpret_cast<real*>(&data.template get<LTS::LocalIntegration>()));
#ifdef USE_VISCOELASTIC2
          auto* dofsExt = currentLayer->var<LTS::DofsExtScratch>(cfg, AllocationPlace::Device);
          dofsExtPtrs.push_back(static_cast<real*>(dofsExt) + tensor::Qext<Cfg>::size() * cell);
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
  });
}

void LocalIntegrationRecorder::recordDisplacements() {
  currentLayer->wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;

    real*(*faceDisplacements)[4] = currentLayer->var<LTS::FaceDisplacementsDevice>(cfg);
    std::array<std::vector<void*>, 4> iVelocitiesPtrs{{}};
    std::array<std::vector<void*>, 4> displacementsPtrs{};

    const auto size = currentLayer->size();
    for (std::size_t cell = 0; cell < size; ++cell) {
      auto dataHost = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Host);

      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        auto isRequired = faceDisplacements[cell][face] != nullptr;
        auto notFreeSurfaceGravity =
            dataHost.template get<LTS::CellInformation>().faceTypes[face] !=
            FaceType::FreeSurfaceGravity;

        if (isRequired && notFreeSurfaceGravity) {
          auto iview =
              init::I<Cfg>::view::create(reinterpret_cast<real*>(idofsAddressRegistry[cell]));
          // NOTE: velocity components are between 6th and 8th columns
          constexpr unsigned FirstVelocityComponent{6};
          iVelocitiesPtrs[face].push_back(
              &multisim::multisimWrap<Cfg>(iview, 0, 0, FirstVelocityComponent));
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
  });
}

void LocalIntegrationRecorder::recordFreeSurfaceGravityBc() {
  currentLayer->wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;

    const auto size = currentLayer->size();
    constexpr size_t NodalAvgDisplacementsSize = tensor::averageNormalDisplacement<Cfg>::size();

    real* nodalAvgDisplacements = static_cast<real*>(
        currentLayer->var<LTS::NodalAvgDisplacements>(cfg, AllocationPlace::Device));

    real* rotateDisplacementToFaceNormalScratch =
        static_cast<real*>(currentLayer->var<LTS::RotateDisplacementToFaceNormalScratch>(
            cfg, AllocationPlace::Device));
    real* rotateDisplacementToGlobalScratch = static_cast<real*>(
        currentLayer->var<LTS::RotateDisplacementToGlobalScratch>(cfg, AllocationPlace::Device));
    real* rotatedFaceDisplacementScratch = static_cast<real*>(
        currentLayer->var<LTS::RotatedFaceDisplacementScratch>(cfg, AllocationPlace::Device));
    real* dofsFaceNodalScratch = static_cast<real*>(
        currentLayer->var<LTS::DofsFaceNodalScratch>(cfg, AllocationPlace::Device));
    real* prevCoefficientsScratch = static_cast<real*>(
        currentLayer->var<LTS::PrevCoefficientsScratch>(cfg, AllocationPlace::Device));

    real* dofsFaceBoundaryNodalScratch = static_cast<real*>(
        currentLayer->var<LTS::DofsFaceBoundaryNodalScratch>(cfg, AllocationPlace::Device));

    if (size > 0) {
      std::array<std::vector<unsigned>, 4> cellIndices{};
      std::array<std::vector<void*>, 4> nodalAvgDisplacementsPtrs{};
      std::array<std::vector<void*>, 4> displacementsPtrs{};

      std::array<std::vector<void*>, 4> derivatives{};
      std::array<std::vector<void*>, 4> dofsPtrs{};
      std::array<std::vector<void*>, 4> idofsPtrs{};
      std::array<std::vector<void*>, 4> neighPtrs{};
      std::array<std::vector<void*>, 4> t{};
      std::array<std::vector<void*>, 4> tInv{};

      std::array<std::vector<inner_keys::Material::DataType>, 4> rhos;
      std::array<std::vector<inner_keys::Material::DataType>, 4> lambdas;

      std::array<std::vector<void*>, 4> rotateDisplacementToFaceNormalPtrs{};
      std::array<std::vector<void*>, 4> rotateDisplacementToGlobalPtrs{};
      std::array<std::vector<void*>, 4> rotatedFaceDisplacementPtrs{};
      std::array<std::vector<void*>, 4> dofsFaceNodalPtrs{};
      std::array<std::vector<void*>, 4> prevCoefficientsPtrs{};
      std::array<std::vector<double>, 4> invImpedances{};
      std::array<std::vector<void*>, 4> dofsFaceBoundaryNodalPtrs{};

      std::array<std::size_t, 4> counter{};

      size_t nodalAvgDisplacementsCounter{0};

      for (std::size_t cell = 0; cell < size; ++cell) {
        auto data = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Device);
        auto dataHost = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Host);

        for (std::size_t face = 0; face < 4; ++face) {
          if (dataHost.template get<LTS::CellInformation>().faceTypes[face] ==
              FaceType::FreeSurfaceGravity) {
            assert(dataHost.template get<LTS::FaceDisplacementsDevice>()[face] != nullptr);
            cellIndices[face].push_back(cell);

            derivatives[face].push_back(dQPtrs[cell]);
            dofsPtrs[face].push_back(static_cast<real*>(data.template get<LTS::Dofs>()));
            idofsPtrs[face].push_back(idofsAddressRegistry[cell]);

            neighPtrs[face].push_back(
                reinterpret_cast<real*>(&data.template get<LTS::NeighboringIntegration>()));
            displacementsPtrs[face].push_back(
                dataHost.template get<LTS::FaceDisplacementsDevice>()[face]);
            t[face].push_back(dataHost.template get<LTS::BoundaryMappingDevice>()[face].dataT);
            tInv[face].push_back(
                dataHost.template get<LTS::BoundaryMappingDevice>()[face].dataTinv);

            rhos[face].push_back(dataHost.template get<LTS::Material>().local->getDensity());
            lambdas[face].push_back(dataHost.template get<LTS::Material>().local->getLambdaBar());

            real* displ{&nodalAvgDisplacements[nodalAvgDisplacementsCounter]};
            nodalAvgDisplacementsPtrs[face].push_back(displ);
            nodalAvgDisplacementsCounter += NodalAvgDisplacementsSize;

            rotateDisplacementToFaceNormalPtrs[face].push_back(
                rotateDisplacementToFaceNormalScratch +
                counter[face] * init::displacementRotationMatrix<Cfg>::Size);
            rotateDisplacementToGlobalPtrs[face].push_back(
                rotateDisplacementToGlobalScratch +
                counter[face] * init::displacementRotationMatrix<Cfg>::Size);
            rotatedFaceDisplacementPtrs[face].push_back(
                rotatedFaceDisplacementScratch +
                counter[face] * init::rotatedFaceDisplacement<Cfg>::Size);
            dofsFaceBoundaryNodalPtrs[face].push_back(dofsFaceBoundaryNodalScratch +
                                                      counter[face] * tensor::INodal<Cfg>::size());
            dofsFaceNodalPtrs[face].push_back(dofsFaceNodalScratch +
                                              counter[face] * tensor::INodal<Cfg>::size());
            prevCoefficientsPtrs[face].push_back(
                prevCoefficientsScratch +
                counter[face] * nodal::tensor::nodes2D<Cfg>::Shape[multisim::BasisDim<Cfg>]);
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
          (*currentTable)[key].set(inner_keys::Wp::Id::PrevCoefficients,
                                   prevCoefficientsPtrs[face]);
          (*currentMaterialTable)[key].set(inner_keys::Material::Id::InvImpedances,
                                           invImpedances[face]);
          (*currentTable)[key].set(inner_keys::Wp::Id::DofsFaceBoundaryNodal,
                                   dofsFaceBoundaryNodalPtrs[face]);
        }
      }
    }
  });
}

void LocalIntegrationRecorder::recordDirichletBc() {
  currentLayer->wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;

    const auto size = currentLayer->size();
    if (size > 0) {
      std::array<std::vector<void*>, 4> dofsPtrs{};
      std::array<std::vector<void*>, 4> idofsPtrs{};
      std::array<std::vector<void*>, 4> tInv{};
      std::array<std::vector<void*>, 4> neighPtrs{};

      std::array<std::vector<void*>, 4> easiBoundaryMapPtrs{};
      std::array<std::vector<void*>, 4> easiBoundaryConstantPtrs{};

      std::array<std::vector<void*>, 4> dofsFaceBoundaryNodalPtrs{};

      std::array<std::size_t, 4> counter{};

      real* dofsFaceBoundaryNodalScratch = static_cast<real*>(
          currentLayer->var<LTS::DofsFaceBoundaryNodalScratch>(cfg, AllocationPlace::Device));

      for (std::size_t cell = 0; cell < size; ++cell) {
        auto data = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Device);
        auto dataHost = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Host);

        for (std::size_t face = 0; face < 4; ++face) {
          if (dataHost.template get<LTS::CellInformation>().faceTypes[face] ==
              FaceType::Dirichlet) {

            dofsPtrs[face].push_back(static_cast<real*>(data.template get<LTS::Dofs>()));
            idofsPtrs[face].push_back(idofsAddressRegistry[cell]);

            tInv[face].push_back(
                dataHost.template get<LTS::BoundaryMappingDevice>()[face].dataTinv);
            neighPtrs[face].push_back(
                reinterpret_cast<real*>(&data.template get<LTS::NeighboringIntegration>()));

            easiBoundaryMapPtrs[face].push_back(
                dataHost.template get<LTS::BoundaryMappingDevice>()[face].easiBoundaryMap);
            easiBoundaryConstantPtrs[face].push_back(
                dataHost.template get<LTS::BoundaryMappingDevice>()[face].easiBoundaryConstant);

            dofsFaceBoundaryNodalPtrs[face].push_back(dofsFaceBoundaryNodalScratch +
                                                      counter[face] * tensor::INodal<Cfg>::size());
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
  });
}

void LocalIntegrationRecorder::recordAnalyticalBc(LTS::Layer& layer) {
  currentLayer->wrap([&](auto cfg) {
    using Cfg = decltype(cfg);
    using real = Real<Cfg>;

    const auto size = currentLayer->size();
    if (size > 0) {
      std::array<std::vector<void*>, 4> dofsPtrs{};
      std::array<std::vector<void*>, 4> neighPtrs{};
      std::array<std::vector<unsigned>, 4> cellIndices{};
      std::array<std::vector<void*>, 4> analytical{};

      real* analyticScratch =
          reinterpret_cast<real*>(layer.var<LTS::AnalyticScratch>(AllocationPlace::Device));

      for (std::size_t cell = 0; cell < size; ++cell) {
        auto data = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Device);
        auto dataHost = currentLayer->cellRef<Cfg>(cell, AllocationPlace::Host);

        for (std::size_t face = 0; face < 4; ++face) {
          if (dataHost.template get<LTS::CellInformation>().faceTypes[face] ==
              FaceType::Analytical) {
            cellIndices[face].push_back(cell);
            dofsPtrs[face].push_back(data.template get<LTS::Dofs>());
            neighPtrs[face].push_back(
                reinterpret_cast<real*>(&data.template get<LTS::NeighboringIntegration>()));
            analytical[face].push_back(analyticScratch + cell * tensor::INodal<Cfg>::size());
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
  });
}
