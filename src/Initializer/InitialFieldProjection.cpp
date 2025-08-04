// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "InitialFieldProjection.h"

#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/MemoryManager.h"
#include "Numerical/Quadrature.h"
#include "Numerical/Transformation.h"
#include "ParameterDB.h"

#include "Initializer/PreProcessorMacros.h"
#include <Alignment.h>
#include <Common/Constants.h>
#include <Equations/Datastructures.h>
#include <Geometry/MeshReader.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Layer.h>
#include <Physics/InitialField.h>
#include <Solver/MultipleSimulations.h>

#include <GeneratedCode/init.h>
#include <array>
#include <cstddef>
#include <easi/Query.h>
#include <easi/ResultAdapter.h>
#include <easi/YAMLParser.h>
#include <memory>
#include <string>
#include <vector>

#ifdef USE_ASAGI
#include <Reader/AsagiReader.h>
#include <easi/util/AsagiReader.h>
#endif

// time-dependent conditions require easi version 1.5.0 or higher
#ifdef EASI_VERSION_MAJOR
#if EASI_VERSION_MINOR >= 5 || EASI_VERSION_MAJOR > 1
#define SUPPORTS_EASI_TIME
#endif
#endif

#ifdef SUPPORTS_EASI_TIME
#include <set>
#else
#include <utils/logger.h>
#endif

GENERATE_HAS_MEMBER(selectAneFull)
GENERATE_HAS_MEMBER(selectElaFull)
GENERATE_HAS_MEMBER(Values)
GENERATE_HAS_MEMBER(Qane)

namespace seissol::init {
template<typename> class selectAneFull;
template<typename> class selectElaFull;
} // namespace seissol::init

#ifndef USE_ASAGI
namespace easi {
class AsagiReader {};
} // namespace easi
#endif

namespace {
struct EasiLoader {
  bool hasTime;
  std::vector<std::unique_ptr<easi::Component>> components;
  std::unique_ptr<easi::AsagiReader> asagiReader;
  std::unique_ptr<easi::YAMLParser> parser;
  EasiLoader(bool hasTime, const std::vector<std::string>& files) : hasTime(hasTime) {
#ifdef USE_ASAGI
    asagiReader = std::make_unique<seissol::asagi::AsagiReader>();
#else
    asagiReader.reset();
#endif

    // NOTE: easi currently sorts the dimension names lexicographically (due to using std::set)
    // hence: if we have time as a dimension, it will come first
#ifdef SUPPORTS_EASI_TIME
    const auto dimensionNames =
        hasTime ? std::set<std::string>{"t", "x", "y", "z"} : std::set<std::string>{"x", "y", "z"};
    parser = std::make_unique<easi::YAMLParser>(dimensionNames, asagiReader.get());
#else
    // ignore time
    if (hasTime) {
      logError() << "easi is too old for time-dependent initial conditions. You need at least "
                    "version 1.5.0.";
    }
    parser = std::make_unique<easi::YAMLParser>(3, asagiReader.get(), 'x');
#endif
    components.resize(files.size());
    for (std::size_t i = 0; i < files.size(); ++i) {
      components[i] = std::unique_ptr<easi::Component>(parser->parse(files.at(i)));
    }
  }
};
} // namespace

namespace seissol::initializer {

void projectInitialField(const std::vector<std::unique_ptr<physics::InitialField>>& iniFields,
                         const GlobalData& globalData,
                         const seissol::geometry::MeshReader& meshReader,
                         seissol::initializer::MemoryManager& memoryManager,
                         LTS::Storage& storage) {
  const auto& vertices = meshReader.getVertices();
  const auto& elements = meshReader.getElements();

  for (auto& layer : storage.leaves(Ghost)) {
    layer.wrap([&](auto cfg) {
      using Cfg = decltype(cfg);

      const auto& global = globalData.get<Cfg>();

      constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
      constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

      double quadraturePoints[NumQuadPoints][Cell::Dim];
      double quadratureWeights[NumQuadPoints];
      seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, QuadPolyDegree);
#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel
#endif
    {
      alignas(Alignment) real iniCondData[tensor::iniCond<Cfg>::size()] = {};
      auto iniCond = init::iniCond<Cfg>::view::create(iniCondData);

      std::vector<std::array<double, Cell::Dim>> quadraturePointsXyz;
      quadraturePointsXyz.resize(NumQuadPoints);

      kernel::projectIniCond<Cfg> krnl;
      krnl.projectQP = global.projectQPMatrix;
      krnl.iniCond = iniCondData;
      kernels::set_selectAneFull(krnl, kernels::get_static_ptr_Values<init::selectAneFull<Cfg>>());
      kernels::set_selectElaFull(krnl, kernels::get_static_ptr_Values<init::selectElaFull<Cfg>>());

      const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
      const auto* material = layer.var<LTS::Material>();
      auto* dofs = layer.var<LTS::Dofs>(cfg);
      auto* dofsAne = layer.var<LTS::DofsAne>(cfg);

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp for schedule(static)
#endif
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        const auto meshId = secondaryInformation[cell].meshId;
        const double* elementCoords[Cell::NumVertices];
        for (size_t v = 0; v < Cell::NumVertices; ++v) {
          elementCoords[v] = vertices[elements[meshId].vertices[v]].coords;
        }
        for (size_t i = 0; i < NumQuadPoints; ++i) {
          seissol::transformations::tetrahedronReferenceToGlobal(elementCoords[0],
                                                                 elementCoords[1],
                                                                 elementCoords[2],
                                                                 elementCoords[3],
                                                                 quadraturePoints[i],
                                                                 quadraturePointsXyz[i].data());
        }

        const CellMaterialData& materialData = material[cell];
        for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
          auto sub = multisim::simtensor(iniCond, s);
          iniFields[s % iniFields.size()]->evaluate(
              0.0, quadraturePointsXyz.data(), quadraturePointsXyz.size(), materialData, sub);
        }

        krnl.Q = dofs[cell];
        if constexpr (kernels::HasSize<tensor::Qane<Cfg>>::Value) {
          kernels::set_Qane(krnl, dofsAne[cell]);
        }
        krnl.execute();
      }
    }
  });
  }
}

std::vector<double> projectEasiFields(const std::vector<std::string>& iniFields,
                                      double time,
                                      const seissol::geometry::MeshReader& meshReader,
                                      bool needsTime) {
  const auto& vertices = meshReader.getVertices();
  const auto& elements = meshReader.getElements();

  constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const int dimensions = needsTime ? (Cell::Dim + 1) : Cell::Dim;
  const int spaceStart = needsTime ? 1 : 0;
  easi::Query query(elements.size() * NumQuadPoints, dimensions);

  {
    double quadraturePoints[NumQuadPoints][Cell::Dim];
    double quadratureWeights[NumQuadPoints];
    seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, QuadPolyDegree);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t elem = 0; elem < elements.size(); ++elem) {
      const double* elementCoords[Cell::NumVertices];
      for (size_t v = 0; v < Cell::NumVertices; ++v) {
        elementCoords[v] = vertices[elements[elem].vertices[v]].coords;
      }
      for (size_t i = 0; i < NumQuadPoints; ++i) {
        std::array<double, Cell::Dim> transformed;
        seissol::transformations::tetrahedronReferenceToGlobal(elementCoords[0],
                                                               elementCoords[1],
                                                               elementCoords[2],
                                                               elementCoords[3],
                                                               quadraturePoints[i],
                                                               transformed.data());
        for (std::size_t d = 0; d < Cell::Dim; ++d) {
          query.x(elem * NumQuadPoints + i, spaceStart + d) = transformed[d];
        }
        if (needsTime) {
          query.x(elem * NumQuadPoints + i, 0) = 0;
        }
        query.group(elem * NumQuadPoints + i) = elements[elem].group;
      }
    }
  }

  std::vector<double> data(NumQuadPoints * iniFields.size() * model::MaterialTT<Cfg>::Quantities.size() *
                           elements.size());
  const auto dataPointStride = iniFields.size() * model::MaterialTT<Cfg>::Quantities.size();
  {
    auto models = EasiLoader(needsTime, iniFields);
    for (std::size_t i = 0; i < iniFields.size(); ++i) {
      auto adapter = easi::ArraysAdapter();
      for (std::size_t j = 0; j < model::MaterialTT<Cfg>::Quantities.size(); ++j) {
        const auto& quantity = model::MaterialTT<Cfg>::Quantities.at(j);
        const std::size_t bindOffset = i + j * iniFields.size();
        adapter.addBindingPoint(quantity, data.data() + bindOffset, dataPointStride);
      }
      models.components.at(i)->evaluate(query, adapter);
    }
  }

  return data;
}

void projectEasiInitialField(const std::vector<std::string>& iniFields,
                             const GlobalData& globalData,
                             const seissol::geometry::MeshReader& meshReader,
                             seissol::initializer::MemoryManager& memoryManager,
                             LTS::Storage& storage,
                             bool needsTime) {
  const auto data = projectEasiFields(iniFields, 0, meshReader, needsTime);

  for (auto& layer : storage.leaves(Ghost)) {
    layer.wrap([&](auto cfg) {
      using Cfg = decltype(cfg);
      constexpr auto QuadPolyDegree = Cfg::ConvergenceOrder + 1;
      constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;
      const auto dataStride = NumQuadPoints * iniFields.size() * model::MaterialTT<Cfg>::Quantities.size();
      const auto quantityCount = model::MaterialTT<Cfg>::Quantities.size();

      const auto& global = globalData.get<Cfg>();

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel
#endif
    {
      alignas(Alignment) real iniCondData[tensor::iniCond<Cfg>::size()] = {};
      auto iniCond = init::iniCond<Cfg>::view::create(iniCondData);

      std::vector<std::array<double, 3>> quadraturePointsXyz;
      quadraturePointsXyz.resize(NumQuadPoints);

      kernel::projectIniCond<Cfg> krnl;
      krnl.projectQP = global.projectQPMatrix;
      krnl.iniCond = iniCondData;
      kernels::set_selectAneFull(krnl, kernels::get_static_ptr_Values<init::selectAneFull<Cfg>>());
      kernels::set_selectElaFull(krnl, kernels::get_static_ptr_Values<init::selectElaFull<Cfg>>());

      const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
      auto* dofs = layer.var<LTS::Dofs>(cfg);
      auto* dofsAne = layer.var<LTS::DofsAne>(cfg);

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp for schedule(static)
#endif
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        const auto meshId = secondaryInformation[cell].meshId;
        // TODO: multisim loop

        for (std::size_t s = 0; s < seissol::multisim::NumSimulations; s++) {
          auto sub = multisim::simtensor(iniCond, s);
          for (std::size_t i = 0; i < NumQuadPoints; ++i) {
            for (std::size_t j = 0; j < quantityCount; ++j) {
              sub(i, j) = data.at(meshId * dataStride + quantityCount * i + j);
            }
          }
        }

        krnl.Q = dofs[cell];
        if constexpr (kernels::HasSize<tensor::Qane<Cfg>>::Value) {
          kernels::set_Qane(krnl, dofsAne[cell]);
        }
        krnl.execute();
      }
    }
  });
  }
}

} // namespace seissol::initializer
