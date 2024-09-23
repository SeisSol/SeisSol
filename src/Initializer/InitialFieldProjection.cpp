// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 */

#include "InitialFieldProjection.h"

#include "Initializer/Tree/LTSSync.h"

#include "Initializer/MemoryManager.h"
#include "Numerical/Quadrature.h"
#include "Numerical/Transformation.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"

#include "Initializer/PreProcessorMacros.h"
#include <Common/Constants.h>
#include <Geometry/MeshReader.h>
#include <Initializer/LTS.h>
#include <Initializer/Tree/Lut.h>
#include <Initializer/Typedefs.h>
#include <Kernels/Common.h>
#include <Kernels/Precision.h>
#include <Physics/InitialField.h>
#include <array>
#include <cstddef>
#include <init.h>
#include <memory>
#include <vector>

GENERATE_HAS_MEMBER(selectAneFull)
GENERATE_HAS_MEMBER(selectElaFull)
GENERATE_HAS_MEMBER(Values)
GENERATE_HAS_MEMBER(Qane)

namespace seissol::init {
class selectAneFull;
class selectElaFull;
} // namespace seissol::init

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

} // namespace seissol::initializer
