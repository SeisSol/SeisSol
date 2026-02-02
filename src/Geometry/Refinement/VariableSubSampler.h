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
  std::vector<basisFunction::SampledBasisFunctions<T>> BasisFunctions_;

  /** The original number of cells (without refinement) */
  unsigned int mNumCells_;

  std::size_t kSubCellsPerCell_;
  std::size_t kNumVariables_;
  std::size_t kNumAlignedDOF_;

  bool nodal_{false};

  std::size_t
      getInVarOffset(std::size_t cell, std::size_t variable, const unsigned int* cellMap) const {
    return (cellMap[cell] * kNumVariables_ + variable) * kNumAlignedDOF_;
  }

  [[nodiscard]] std::size_t getOutVarOffset(std::size_t cell, std::size_t subcell) const {
    return kSubCellsPerCell_ * cell + subcell;
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
    : mNumCells_(numCells), kSubCellsPerCell_(tetRefiner.getDivisionCount()),
      kNumVariables_(numVariables), kNumAlignedDOF_(numAlignedDOF), nodal_(nodal) {
  // Generate cell centerpoints in the reference or unit tetrahedron.
  auto* subCells = new Tetrahedron<T>[kSubCellsPerCell_];
  auto* additionalVertices = new Eigen::Matrix<T, 3, 1>[tetRefiner.additionalVerticesPerCell()];

  tetRefiner.refine(Tetrahedron<T>::unitTetrahedron(), 0, subCells, additionalVertices);

  // Generate sampled basicfunctions
  for (unsigned int i = 0; i < kSubCellsPerCell_; i++) {
    const Eigen::Matrix<T, 3, 1> pnt = subCells[i].center();
    BasisFunctions_.push_back(
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
  for (unsigned int c = 0; c < mNumCells_; ++c) {
    for (unsigned int sc = 0; sc < kSubCellsPerCell_; ++sc) {
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
      outData[getOutVarOffset(c, sc)] = BasisFunctions_[sc].evalWithCoeffs(inCellData);
    }
  }
}

//------------------------------------------------------------------------------

} // namespace seissol::refinement

#endif // SEISSOL_SRC_GEOMETRY_REFINEMENT_VARIABLESUBSAMPLER_H_
