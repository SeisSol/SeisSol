// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "InitialFieldProjection.h"

#include "Memory/Tree/LTSSync.h"

#include "Initializer/MemoryManager.h"
#include "Numerical/Quadrature.h"
#include "Numerical/Transformation.h"
#include "ParameterDB.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"

#include "Initializer/PreProcessorMacros.h"
#include <Common/Constants.h>
#include <Equations/Datastructures.h>
#include <Geometry/MeshReader.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/Tree/Lut.h>
#include <Physics/InitialField.h>
#include <array>
#include <cstddef>
#include <easi/Query.h>
#include <easi/ResultAdapter.h>
#include <easi/YAMLParser.h>
#include <init.h>
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
class selectAneFull;
class selectElaFull;
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
    asagiReader = std::make_unique<seissol::asagi::AsagiReader>("SEISSOL_ASAGI");
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
                         LTS const& lts,
                         const Lut& ltsLut) {
  const auto& vertices = meshReader.getVertices();
  const auto& elements = meshReader.getElements();

  constexpr auto QuadPolyDegree = ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  double quadraturePoints[NumQuadPoints][3];
  double quadratureWeights[NumQuadPoints];
  seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, QuadPolyDegree);

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel
  {
#endif
    alignas(Alignment) real iniCondData[tensor::iniCond::size()] = {};
    auto iniCond = init::iniCond::view::create(iniCondData);

    std::vector<std::array<double, 3>> quadraturePointsXyz;
    quadraturePointsXyz.resize(NumQuadPoints);

    kernel::projectIniCond krnl;
    krnl.projectQP = globalData.projectQPMatrix;
    krnl.iniCond = iniCondData;
    kernels::set_selectAneFull(krnl, kernels::get_static_ptr_Values<init::selectAneFull>());
    kernels::set_selectElaFull(krnl, kernels::get_static_ptr_Values<init::selectElaFull>());

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp for schedule(static)
#endif
    for (unsigned int meshId = 0; meshId < elements.size(); ++meshId) {
      const double* elementCoords[4];
      for (size_t v = 0; v < 4; ++v) {
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

      const CellMaterialData& material = ltsLut.lookup(lts.material, meshId);
#ifdef MULTIPLE_SIMULATIONS
      for (int s = 0; s < MULTIPLE_SIMULATIONS; ++s) {
        auto sub = iniCond.subtensor(s, yateto::slice<>(), yateto::slice<>());
        iniFields[s % iniFields.size()]->evaluate(0.0, quadraturePointsXyz, material, sub);
      }
#else
    iniFields[0]->evaluate(0.0, quadraturePointsXyz, material, iniCond);
#endif

      krnl.Q = ltsLut.lookup(lts.dofs, meshId);
      if constexpr (kernels::HasSize<tensor::Qane>::Value) {
        kernels::set_Qane(krnl, &ltsLut.lookup(lts.dofsAne, meshId)[0]);
      }
      krnl.execute();
    }
#if defined(_OPENMP) && !NVHPC_AVOID_OMP
  }
#endif

  seissol::initializer::synchronizeLTSTreeDuplicates(lts.dofs, memoryManager);
  if (kernels::size<tensor::Qane>() > 0) {
    seissol::initializer::synchronizeLTSTreeDuplicates(lts.dofsAne, memoryManager);
  }
}

std::vector<double> projectEasiFields(const std::vector<std::string>& iniFields,
                                      double time,
                                      const seissol::geometry::MeshReader& meshReader,
                                      bool needsTime) {
  const auto& vertices = meshReader.getVertices();
  const auto& elements = meshReader.getElements();

  constexpr auto QuadPolyDegree = ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const int dimensions = needsTime ? 4 : 3;
  const int spaceStart = needsTime ? 1 : 0;
  easi::Query query(elements.size() * NumQuadPoints, dimensions);

  {
    double quadraturePoints[NumQuadPoints][3];
    double quadratureWeights[NumQuadPoints];
    seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, QuadPolyDegree);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (std::size_t elem = 0; elem < elements.size(); ++elem) {
      const double* elementCoords[4];
      for (size_t v = 0; v < 4; ++v) {
        elementCoords[v] = vertices[elements[elem].vertices[v]].coords;
      }
      for (size_t i = 0; i < NumQuadPoints; ++i) {
        std::array<double, 3> transformed;
        seissol::transformations::tetrahedronReferenceToGlobal(elementCoords[0],
                                                               elementCoords[1],
                                                               elementCoords[2],
                                                               elementCoords[3],
                                                               quadraturePoints[i],
                                                               transformed.data());
        query.x(elem * NumQuadPoints + i, spaceStart + 0) = transformed[0];
        query.x(elem * NumQuadPoints + i, spaceStart + 1) = transformed[1];
        query.x(elem * NumQuadPoints + i, spaceStart + 2) = transformed[2];
        if (needsTime) {
          query.x(elem * NumQuadPoints + i, 0) = 0;
        }
        query.group(elem * NumQuadPoints + i) = elements[elem].group;
      }
    }
  }

  std::vector<double> data(NumQuadPoints * iniFields.size() * model::MaterialT::Quantities.size() *
                           elements.size());
  const auto dataPointStride = iniFields.size() * model::MaterialT::Quantities.size();
  {
    auto models = EasiLoader(needsTime, iniFields);
    for (std::size_t i = 0; i < iniFields.size(); ++i) {
      auto adapter = easi::ArraysAdapter();
      for (std::size_t j = 0; j < model::MaterialT::Quantities.size(); ++j) {
        const auto& quantity = model::MaterialT::Quantities.at(j);
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
                             LTS const& lts,
                             const Lut& ltsLut,
                             bool needsTime) {
  const auto& vertices = meshReader.getVertices();
  const auto& elements = meshReader.getElements();

  constexpr auto QuadPolyDegree = ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const auto data = projectEasiFields(iniFields, 0, meshReader, needsTime);

  const auto dataStride = NumQuadPoints * iniFields.size() * model::MaterialT::Quantities.size();
  const auto quantityCount = model::MaterialT::Quantities.size();

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp parallel
  {
#endif
    alignas(Alignment) real iniCondData[tensor::iniCond::size()] = {};
    auto iniCond = init::iniCond::view::create(iniCondData);

    std::vector<std::array<double, 3>> quadraturePointsXyz;
    quadraturePointsXyz.resize(NumQuadPoints);

    kernel::projectIniCond krnl;
    krnl.projectQP = globalData.projectQPMatrix;
    krnl.iniCond = iniCondData;
    kernels::set_selectAneFull(krnl, kernels::get_static_ptr_Values<init::selectAneFull>());
    kernels::set_selectElaFull(krnl, kernels::get_static_ptr_Values<init::selectElaFull>());

#if defined(_OPENMP) && !NVHPC_AVOID_OMP
#pragma omp for schedule(static)
#endif
    for (unsigned int meshId = 0; meshId < elements.size(); ++meshId) {
      // TODO: multisim loop
      for (std::size_t i = 0; i < NumQuadPoints; ++i) {
        for (std::size_t j = 0; j < quantityCount; ++j) {
          iniCond(i, j) = data.at(meshId * dataStride + quantityCount * i + j);
        }
      }

      krnl.Q = ltsLut.lookup(lts.dofs, meshId);
      if constexpr (kernels::HasSize<tensor::Qane>::Value) {
        kernels::set_Qane(krnl, &ltsLut.lookup(lts.dofsAne, meshId)[0]);
      }
      krnl.execute();
    }
#if defined(_OPENMP) && !NVHPC_AVOID_OMP
  }
#endif

  seissol::initializer::synchronizeLTSTreeDuplicates(lts.dofs, memoryManager);
  if (kernels::size<tensor::Qane>() > 0) {
    seissol::initializer::synchronizeLTSTreeDuplicates(lts.dofsAne, memoryManager);
  }
}

} // namespace seissol::initializer
