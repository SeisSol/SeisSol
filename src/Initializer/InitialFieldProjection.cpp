// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#include "InitialFieldProjection.h"

#include "Alignment.h"
#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/MeshReader.h"
#include "Initializer/PreProcessorMacros.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Numerical/Quadrature.h"
#include "Numerical/Transformation.h"
#include "ParameterDB.h"
#include "Physics/InitialField.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/ReaderBuilder.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cstddef>
#include <exception>
#include <memory>
#include <string>
#include <utils/logger.h>
#include <vector>

GENERATE_HAS_MEMBER(selectAneFull)
GENERATE_HAS_MEMBER(selectElaFull)
GENERATE_HAS_MEMBER(Values)
GENERATE_HAS_MEMBER(Qane)

namespace seissol::init {
struct selectAneFull;
struct selectElaFull;
} // namespace seissol::init

namespace seissol::initializer {

void projectInitialField(const std::vector<std::unique_ptr<physics::InitialField>>& iniFields,
                         const GlobalData& globalData,
                         const seissol::geometry::MeshReader& meshReader,
                         LTS::Storage& storage) {
  const auto& vertices = meshReader.getVertices();
  const auto& elements = meshReader.getElements();

  constexpr auto QuadPolyDegree = ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  double quadraturePoints[NumQuadPoints][Cell::Dim];
  double quadratureWeights[NumQuadPoints];
  seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, QuadPolyDegree);

  for (auto& layer : storage.leaves(Ghost)) {
#if !NVHPC_AVOID_OMP
#pragma omp parallel
#endif
    {
      alignas(Alignment) real iniCondData[tensor::iniCond::size()] = {};
      auto iniCond = init::iniCond::view::create(iniCondData);

      std::vector<std::array<double, Cell::Dim>> quadraturePointsXyz;
      quadraturePointsXyz.resize(NumQuadPoints);

      kernel::projectIniCond krnl;
      krnl.projectQP = globalData.projectQPMatrix;
      krnl.iniCond = iniCondData;
      set_selectAneFull(krnl, get_static_ptr_Values<init::selectAneFull>());
      set_selectElaFull(krnl, get_static_ptr_Values<init::selectElaFull>());

      const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
      const auto* material = layer.var<LTS::Material>();
      auto* dofs = layer.var<LTS::Dofs>();
      auto* dofsAne = layer.var<LTS::DofsAne>();

#if !NVHPC_AVOID_OMP
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
        if constexpr (kernels::HasSize<tensor::Qane>::Value) {
          set_Qane(krnl, dofsAne[cell]);
        }
        krnl.execute();
      }
    }
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

  std::vector<std::array<double, Cell::Dim>> points(elements.size() * NumQuadPoints);

  const auto inVars = needsTime ? std::vector<std::string>{"t", "x", "y", "z"}
                                : std::vector<std::string>{"x", "y", "z"};

  double quadraturePoints[NumQuadPoints][Cell::Dim];
  double quadratureWeights[NumQuadPoints];
  seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, QuadPolyDegree);

#pragma omp parallel for schedule(static)
  for (std::size_t elem = 0; elem < elements.size(); ++elem) {
    const double* elementCoords[Cell::NumVertices];
    for (size_t v = 0; v < Cell::NumVertices; ++v) {
      elementCoords[v] = vertices[elements[elem].vertices[v]].coords;
    }
    for (size_t i = 0; i < NumQuadPoints; ++i) {
      const auto idx = elem * NumQuadPoints + i;

      seissol::transformations::tetrahedronReferenceToGlobal(elementCoords[0],
                                                             elementCoords[1],
                                                             elementCoords[2],
                                                             elementCoords[3],
                                                             quadraturePoints[i],
                                                             points[idx].data());
    }
  }

  std::vector<double> data(NumQuadPoints * iniFields.size() * model::MaterialT::Quantities.size() *
                           elements.size());
  const auto dataPointStride = iniFields.size() * model::MaterialT::Quantities.size();

  for (std::size_t i = 0; i < iniFields.size(); ++i) {
    reader::scripting::DataTable table(elements.size() * NumQuadPoints);
    table.bindComputed("__group", [&](std::size_t index) -> std::int32_t {
      return elements[index / NumQuadPoints].group;
    });
    table.bindComputed("t", [&](std::size_t) -> double { return time; });
    table.bindComputed("x", [&](std::size_t index) -> double { return points[index][0]; });
    table.bindComputed("y", [&](std::size_t index) -> double { return points[index][1]; });
    table.bindComputed("z", [&](std::size_t index) -> double { return points[index][2]; });
    for (std::size_t j = 0; j < model::MaterialT::Quantities.size(); ++j) {
      const auto& quantity = model::MaterialT::Quantities.at(j);
      const std::size_t bindOffset = i + j * iniFields.size();
      table.bindView(
          quantity, reader::scripting::Direction::Out, data.data(), dataPointStride, bindOffset);
    }

    const auto reader = reader::scripting::buildReader(iniFields[0], inVars);
    reader->call(table);
  }

  return data;
}

void projectEasiInitialField(const std::vector<std::string>& iniFields,
                             const GlobalData& globalData,
                             const seissol::geometry::MeshReader& meshReader,
                             LTS::Storage& storage,
                             bool needsTime) {
  constexpr auto QuadPolyDegree = ConvergenceOrder + 1;
  constexpr auto NumQuadPoints = QuadPolyDegree * QuadPolyDegree * QuadPolyDegree;

  const auto data = projectEasiFields(iniFields, 0, meshReader, needsTime);

  const auto dataStride = NumQuadPoints * iniFields.size() * model::MaterialT::Quantities.size();
  const auto quantityCount = model::MaterialT::Quantities.size();

  for (auto& layer : storage.leaves(Ghost)) {

#if !NVHPC_AVOID_OMP
#pragma omp parallel
#endif
    {
      alignas(Alignment) real iniCondData[tensor::iniCond::size()] = {};
      auto iniCond = init::iniCond::view::create(iniCondData);

      std::vector<std::array<double, 3>> quadraturePointsXyz;
      quadraturePointsXyz.resize(NumQuadPoints);

      kernel::projectIniCond krnl;
      krnl.projectQP = globalData.projectQPMatrix;
      krnl.iniCond = iniCondData;
      set_selectAneFull(krnl, get_static_ptr_Values<init::selectAneFull>());
      set_selectElaFull(krnl, get_static_ptr_Values<init::selectElaFull>());

      const auto* secondaryInformation = layer.var<LTS::SecondaryInformation>();
      auto* dofs = layer.var<LTS::Dofs>();
      auto* dofsAne = layer.var<LTS::DofsAne>();

#if !NVHPC_AVOID_OMP
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
        if constexpr (kernels::HasSize<tensor::Qane>::Value) {
          set_Qane(krnl, dofsAne[cell]);
        }
        krnl.execute();
      }
    }
  }
}

} // namespace seissol::initializer
