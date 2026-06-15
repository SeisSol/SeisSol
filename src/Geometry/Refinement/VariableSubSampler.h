// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_GEOMETRY_REFINEMENT_VARIABLESUBSAMPLER_H_
#define SEISSOL_SRC_GEOMETRY_REFINEMENT_VARIABLESUBSAMPLER_H_

#include "Geometry/MeshReader.h"
#include "Numerical/BasisFunction.h"
#include "RefinerUtils.h"

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cstddef>

namespace seissol::refinement {

//------------------------------------------------------------------------------

template <class T>
class VariableSubsampler {
  private:
  std::vector<basisFunction::SampledBasisFunctions<T>> basisFunctions_;

  /** The original number of cells (without refinement) */
  unsigned int numCells_;

  std::size_t subCellsPerCell_;
  std::size_t numVariables_;
  std::size_t numAlignedDOF_;

  bool nodal_{false};

  std::size_t
      getInVarOffset(std::size_t cell, std::size_t variable, const unsigned int* cellMap) const {
    return (cellMap[cell] * numVariables_ + variable) * numAlignedDOF_;
  }

  [[nodiscard]] std::size_t getOutVarOffset(std::size_t cell, std::size_t subcell) const {
    return subCellsPerCell_ * cell + subcell;
  }

  public:
  VariableSubsampler(std::size_t numCells,
                     const TetrahedronRefiner<T>& tetRefiner,
                     unsigned int order,
                     std::size_t numVariables,
                     std::size_t numAlignedDOF,
                     bool nodal);

  // NOLINTNEXTLINE
  void get(const real* inData, const unsigned int* cellMap, int variable, real* outData) const;
};

//------------------------------------------------------------------------------

template <typename T>
VariableSubsampler<T>::VariableSubsampler(std::size_t numCells,
                                          const TetrahedronRefiner<T>& tetRefiner,
                                          unsigned int order,
                                          std::size_t numVariables,
                                          std::size_t numAlignedDOF,
                                          bool nodal)
    : numCells_(numCells), subCellsPerCell_(tetRefiner.getDivisionCount()),
      numVariables_(numVariables), numAlignedDOF_(numAlignedDOF), nodal_(nodal) {
  // Generate cell centerpoints in the reference or unit tetrahedron.
  auto* subCells = new Tetrahedron<T>[subCellsPerCell_];
  auto* additionalVertices = new Eigen::Matrix<T, 3, 1>[tetRefiner.additionalVerticesPerCell()];

  tetRefiner.refine(Tetrahedron<T>::unitTetrahedron(), 0, subCells, additionalVertices);

  // Generate sampled basicfunctions
  for (unsigned int i = 0; i < subCellsPerCell_; i++) {
    const Eigen::Matrix<T, 3, 1> pnt = subCells[i].center();
    basisFunctions_.push_back(
        basisFunction::SampledBasisFunctions<T>(order, pnt(0), pnt(1), pnt(2)));
  }

  delete[] subCells;
  delete[] additionalVertices;
}

//------------------------------------------------------------------------------

template <typename T>
void VariableSubsampler<T>::get(const real* inData,
                                const unsigned int* cellMap,
                                int variable,
                                // NOLINTNEXTLINE
                                real* outData) const {

#pragma omp parallel for schedule(static)
  // Iterate over original Cells
  for (unsigned int c = 0; c < numCells_; ++c) {
    for (unsigned int sc = 0; sc < subCellsPerCell_; ++sc) {
      const real* __restrict inCellData = &inData[getInVarOffset(c, variable, cellMap)];
      alignas(Alignment) real modalBuffer[tensor::modalVar::Size]{};
      if (nodal_) {
        kernel::plOutput krnl{};
        krnl.nodalVar = inCellData;
        krnl.modalVar = modalBuffer;
        krnl.vInv = init::vInv::Values;
        krnl.execute();

        inCellData = modalBuffer;
      }
      outData[getOutVarOffset(c, sc)] = basisFunctions_[sc].evalWithCoeffs(inCellData);
    }
  }
}

//------------------------------------------------------------------------------

} // namespace seissol::refinement

#endif // SEISSOL_SRC_GEOMETRY_REFINEMENT_VARIABLESUBSAMPLER_H_
