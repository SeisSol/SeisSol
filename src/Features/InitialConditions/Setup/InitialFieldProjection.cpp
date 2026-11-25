// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "InitialFieldProjection.h"

#include "Alignment.h"
#include "Core/Common/ConfigDispatch.h"
#include "Core/Common/Constants.h"
#include "Core/Common/MultipleSimulations.h"
#include "Core/Common/PreProcessorMacros.h"
#include "Core/Common/Real.h"
#include "Core/Data/SimulationFiles.h"
#include "Core/Geometry/CellTransform.h"
#include "Core/Geometry/MeshReader.h"
#include "Core/IO/Easi/EasiContext.h"
#include "Core/Kernels/Common.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Numerical/Quadrature.h"
#include "Features/InitialConditions/EasiFields.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "WavePropagation/Physics/Datastructures.h"
#include "WavePropagation/Physics/InitialField.h"

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

GENERATE_HAS_MEMBER(Qane)

namespace seissol::initializer {

namespace {

/// Projects the initial fields onto the cells of `layer`, which compute in the configuration `Cfg`.
template <typename Cfg>
void projectInitialFieldOnLayer(
    const std::vector<std::unique_ptr<physics::InitialField>>& iniFields,
    const seissol::geometry::MeshReader& meshReader,
    LTS::Layer& layer) {
  using real = Real<Cfg>;
  constexpr auto Variant = configIdOf<Cfg>();
  // Looked up rather than named: a configuration without anelastic unknowns has no Qane.
  const auto* anelasticLayout = runtime::tensorTable(Variant).find("Qane", {});

  constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const auto rule = seissol::quadrature::simplexRule<Cell::Dim>(QuadPolyDegree);
  const auto& quadraturePoints = rule.first;

#if !NVHPC_AVOID_OMP
#pragma omp parallel
#endif
  {
    alignas(Alignment) real iniCondData[tensor::iniCond<Cfg>::size()] = {};
    auto iniCond = init::iniCond<Cfg>::view::create(iniCondData);

    std::vector<std::array<double, Cell::Dim>> quadraturePointsXyz;
    quadraturePointsXyz.resize(NumQuadPoints);

    runtime::kernel::projectIniCond krnl;
    krnl.iniCond = runtime::init::iniCond::view(Variant, iniCondData);

    const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    const auto* material = layer.var<LTS::Material>();
    auto* dofs = layer.var<LTS::Dofs>(Cfg());

#if !NVHPC_AVOID_OMP
#pragma omp for schedule(static)
#endif
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      const auto meshId = secondaryInformation[cell].meshId;
      const auto transform = seissol::geometry::AffineTransform::fromMeshCell(meshId, meshReader);
      transform.refToSpace(
          quadraturePoints.data(), quadraturePointsXyz.data(), quadraturePoints.size());

      const CellMaterialData& materialData = material[cell];
      // the copy of the cell holds the simulations from copyId * NumSimulations on
      const auto firstSimulation = secondaryInformation[cell].copyId * Cfg::NumSimulations;
      for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
        auto sub = multisim::simtensor<Cfg>(iniCond, s);
        iniFields[(firstSimulation + s) % iniFields.size()]->evaluate(
            0.0, quadraturePointsXyz.data(), quadraturePointsXyz.size(), materialData, sub);
      }

      krnl.Q = runtime::init::Q::view(Variant, dofs[cell]);
      if constexpr (kernels::HasSize<tensor::Qane<Cfg>>::Value) {
        set_Qane(krnl, yateto::viewOf(anelasticLayout, layer.var<LTS::DofsAne>(Cfg())[cell]));
      }
      krnl.execute(Variant);
    }
  }
}

/// Projects the values `data` of the easi fields, as `projectEasiFields<Cfg>` gives them, onto the
/// cells of `layer`, which compute in the configuration `Cfg`: every simulation takes the field
/// `fieldOfSimulation` names for it.
template <typename Cfg>
void projectEasiFieldsOnLayer(const std::vector<double>& data,
                              std::size_t fieldCount,
                              const std::vector<std::size_t>& fieldOfSimulation,
                              LTS::Layer& layer) {
  using real = Real<Cfg>;
  constexpr auto Variant = configIdOf<Cfg>();
  // Looked up rather than named: a configuration without anelastic unknowns has no Qane.
  const auto* anelasticLayout = runtime::tensorTable(Variant).find("Qane", {});
  constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const auto quantityCount = model::MaterialOf<Cfg>::Quantities.size();
  const auto dataStride = NumQuadPoints * fieldCount * quantityCount;

#if !NVHPC_AVOID_OMP
#pragma omp parallel
#endif
  {
    alignas(Alignment) real iniCondData[tensor::iniCond<Cfg>::size()] = {};
    auto iniCond = init::iniCond<Cfg>::view::create(iniCondData);

    runtime::kernel::projectIniCond krnl;
    krnl.iniCond = runtime::init::iniCond::view(Variant, iniCondData);

    const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
    auto* dofs = layer.var<LTS::Dofs>(Cfg());

#if !NVHPC_AVOID_OMP
#pragma omp for schedule(static)
#endif
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      const auto meshId = secondaryInformation[cell].meshId;
      // the copy of the cell holds the simulations from copyId * NumSimulations on; the fused
      // simulations which only fill up the last copy repeat the first simulations
      const auto firstSimulation = secondaryInformation[cell].copyId * Cfg::NumSimulations;

      for (std::size_t s = 0; s < Cfg::NumSimulations; s++) {
        const auto field = fieldOfSimulation[(firstSimulation + s) % fieldOfSimulation.size()];
        auto sub = multisim::simtensor<Cfg>(iniCond, s);
        for (std::size_t i = 0; i < NumQuadPoints; ++i) {
          for (std::size_t j = 0; j < quantityCount; ++j) {
            sub(i, j) = data.at(meshId * dataStride + (quantityCount * i + j) * fieldCount + field);
          }
        }
      }

      krnl.Q = runtime::init::Q::view(Variant, dofs[cell]);
      if constexpr (kernels::HasSize<tensor::Qane<Cfg>>::Value) {
        set_Qane(krnl, yateto::viewOf(anelasticLayout, layer.var<LTS::DofsAne>(Cfg())[cell]));
      }
      krnl.execute(Variant);
    }
  }
}

} // namespace

void projectInitialField(
    const std::vector<std::vector<std::unique_ptr<physics::InitialField>>>& iniFields,
    const seissol::geometry::MeshReader& meshReader,
    LTS::Storage& storage) {
  for (auto& layer : storage.leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      projectInitialFieldOnLayer<decltype(cfg)>(
          iniFields.at(layer.getIdentifier().config), meshReader, layer);
    });
  }
}

void projectEasiInitialField(const parameters::SimulationFiles& iniFields,
                             const seissol::geometry::MeshReader& meshReader,
                             LTS::Storage& storage,
                             bool needsTime,
                             const data::EasiContext& easi) {
  // the fields are sampled at the points of each configuration, once for all of its layers
  for (const auto config : storage.configs()) {
    dispatchConfig(config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      const auto data = projectEasiFields<Cfg>(iniFields.files, 0, meshReader, needsTime, easi);
      for (auto& layer : storage.leaves(Ghost)) {
        if (layer.getIdentifier().config == config) {
          projectEasiFieldsOnLayer<Cfg>(
              data, iniFields.files.size(), iniFields.fileOfSimulation, layer);
        }
      }
    });
  }
}

} // namespace seissol::initializer
