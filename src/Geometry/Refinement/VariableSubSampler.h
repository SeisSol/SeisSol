// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_GEOMETRY_REFINEMENT_VARIABLESUBSAMPLER_H_
#define SEISSOL_SRC_GEOMETRY_REFINEMENT_VARIABLESUBSAMPLER_H_

#include <algorithm>
#include <cassert>

#include <Eigen/Dense>
#include <cstddef>

#include "Geometry/MeshReader.h"
#include "Numerical/BasisFunction.h"
#include "RefinerUtils.h"

namespace seissol::refinement {

//------------------------------------------------------------------------------

template <class T>
class VariableSubsampler {
  private:
  std::vector<basisFunction::SampledBasisFunctions<T>> m_BasisFunctions;

  /** The original number of cells (without refinement) */
  unsigned int mNumCells;

  unsigned int kSubCellsPerCell;
  unsigned int kNumVariables;
  unsigned int kNumAlignedDOF;

  std::size_t
      getInVarOffset(unsigned int cell, unsigned int variable, const unsigned int* cellMap) const {
    return static_cast<std::size_t>((cellMap[cell] * kNumVariables + variable) * kNumAlignedDOF);
  }

  [[nodiscard]] std::size_t getOutVarOffset(unsigned cell, unsigned int subcell) const {
    return kSubCellsPerCell * cell + subcell;
  }

  public:
  VariableSubsampler(unsigned int numCells,
                     const TetrahedronRefiner<T>& tetRefiner,
                     unsigned int order,
                     unsigned int numVariables,
                     unsigned int numAlignedDOF);

  template <typename RealT>
  void get(const RealT* inData, const unsigned int* cellMap, int variable, RealT* outData) const;
};

//------------------------------------------------------------------------------

template <typename T>
VariableSubsampler<T>::VariableSubsampler(unsigned int numCells,
                                          const TetrahedronRefiner<T>& tetRefiner,
                                          unsigned int order,
                                          unsigned int numVariables,
                                          unsigned int numAlignedDOF)
    : mNumCells(numCells), kSubCellsPerCell(tetRefiner.getDivisionCount()),
      kNumVariables(numVariables), kNumAlignedDOF(numAlignedDOF) {
  // Generate cell centerpoints in the reference or unit tetrahedron.
  auto* subCells = new Tetrahedron<T>[kSubCellsPerCell];
  auto* additionalVertices = new Eigen::Matrix<T, 3, 1>[tetRefiner.additionalVerticesPerCell()];

  tetRefiner.refine(Tetrahedron<T>::unitTetrahedron(), 0, subCells, additionalVertices);

  // Generate sampled basicfunctions
  for (unsigned int i = 0; i < kSubCellsPerCell; i++) {
    const Eigen::Matrix<T, 3, 1> pnt = subCells[i].center();
    m_BasisFunctions.push_back(
        basisFunction::SampledBasisFunctions<T>(order, pnt(0), pnt(1), pnt(2)));
  }

  delete[] subCells;
  delete[] additionalVertices;
}

//------------------------------------------------------------------------------

template <typename T>
template <typename RealT>
void VariableSubsampler<T>::get(const RealT* inData,
                                const unsigned int* cellMap,
                                int variable,
                                RealT* outData) const {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  // Iterate over original Cells
  for (unsigned int c = 0; c < mNumCells; ++c) {
    for (unsigned int sc = 0; sc < kSubCellsPerCell; ++sc) {
      outData[getOutVarOffset(c, sc)] =
          m_BasisFunctions[sc].evalWithCoeffs(&inData[getInVarOffset(c, variable, cellMap)]);
    }
  }
}

//------------------------------------------------------------------------------

} // namespace seissol::refinement

#endif // SEISSOL_SRC_GEOMETRY_REFINEMENT_VARIABLESUBSAMPLER_H_
