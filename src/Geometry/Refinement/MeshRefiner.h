// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_GEOMETRY_REFINEMENT_MESHREFINER_H_
#define SEISSOL_SRC_GEOMETRY_REFINEMENT_MESHREFINER_H_

#include "Geometry/MeshReader.h"
#include "RefinerUtils.h"

#include <cstddef>
#include <cstring>

namespace seissol::refinement {

//------------------------------------------------------------------------------

template <typename T>
class MeshRefiner {
  private:
  // cells_ contains the indices of the cells
  unsigned int* cells_;
  T* vertices_;

  size_t numSubCells_;
  size_t numVertices_;

  static const unsigned int KIndicesPerCell = 4;

  const unsigned int kSubCellsPerCell;

  public:
  MeshRefiner(const seissol::geometry::MeshReader& meshReader,
              const TetrahedronRefiner<T>& tetRefiner);

  MeshRefiner(const std::vector<const Element*>& subElements,
              const std::vector<const Vertex*>& subVertices,
              const std::map<int, int>& oldToNewVertexMap,
              const TetrahedronRefiner<T>& tetRefiner);

  ~MeshRefiner();

  auto operator=(const MeshRefiner&) = delete;
  auto operator=(MeshRefiner&&) = delete;
  MeshRefiner(const MeshRefiner&) = delete;
  MeshRefiner(MeshRefiner&&) = delete;

  [[nodiscard]] const unsigned int* getCellData() const;
  const T* getVertexData() const;
  [[nodiscard]] std::size_t getkSubCellsPerCell() const;

  [[nodiscard]] std::size_t getNumCells() const;
  [[nodiscard]] std::size_t getNumVertices() const;
};

//------------------------------------------------------------------------------

template <typename T>
MeshRefiner<T>::MeshRefiner(const seissol::geometry::MeshReader& meshReader,
                            const TetrahedronRefiner<T>& tetRefiner)
    : kSubCellsPerCell(tetRefiner.getDivisionCount())

{
  using std::size_t;

  const size_t kInVertexCount = meshReader.getVertices().size();
  const size_t kInCellCount = meshReader.getElements().size();
  numSubCells_ = kInCellCount * kSubCellsPerCell;

  const unsigned int additionalVertices = tetRefiner.additionalVerticesPerCell();
  numVertices_ = kInVertexCount + kInCellCount * additionalVertices;

  cells_ = new unsigned int[numSubCells_ * KIndicesPerCell];
  vertices_ = new T[numVertices_ * 3];

  const std::vector<Vertex>& kVertices = meshReader.getVertices();
  const std::vector<Element>& kElements = meshReader.getElements();

  // Copy original vertices

#pragma omp parallel for
  for (unsigned int i = 0; i < kInVertexCount; i++) {
    memcpy(&vertices_[static_cast<size_t>(i * 3)], kVertices[i].coords, sizeof(double) * 3);
  }

  // The pointer to the new vertices
  T* newVertices = &vertices_[kInVertexCount * 3];

  // Start the actual cell-refinement

#pragma omp parallel
  {
    auto* newVerticesTmp = new Eigen::Matrix<T, 3, 1>[additionalVertices];
    auto* newTetsTmp = new Tetrahedron<T>[kSubCellsPerCell];

#pragma omp for schedule(static) nowait
    for (size_t c = 0; c < kInCellCount; ++c) {
      // Build a Terahedron containing the coordinates of the vertices.
      const Tetrahedron<T> inTet = Tetrahedron<T>(kVertices[kElements[c].vertices[0]].coords,
                                                  kVertices[kElements[c].vertices[1]].coords,
                                                  kVertices[kElements[c].vertices[2]].coords,
                                                  kVertices[kElements[c].vertices[3]].coords,
                                                  kElements[c].vertices[0],
                                                  kElements[c].vertices[1],
                                                  kElements[c].vertices[2],
                                                  kElements[c].vertices[3]);

      // Generate the tets
      tetRefiner.refine(inTet, kInVertexCount + c * additionalVertices, newTetsTmp, newVerticesTmp);

      // Copy new vertices
      for (unsigned int i = 0; i < additionalVertices; i++) {
        memcpy(&newVertices[(c * additionalVertices + i) * 3],
               newVerticesTmp[i].data(),
               sizeof(T) * 3);
      }

      // Copy tets
      for (unsigned int i = 0; i < kSubCellsPerCell; i++) {
        cells_[(c * kSubCellsPerCell + i) * 4] = newTetsTmp[i].i;
        cells_[(c * kSubCellsPerCell + i) * 4 + 1] = newTetsTmp[i].j;
        cells_[(c * kSubCellsPerCell + i) * 4 + 2] = newTetsTmp[i].k;
        cells_[(c * kSubCellsPerCell + i) * 4 + 3] = newTetsTmp[i].l;
      }
    }

    delete[] newVerticesTmp;
    delete[] newTetsTmp;
  }
}

template <typename T>
MeshRefiner<T>::MeshRefiner(const std::vector<const Element*>& subElements,
                            const std::vector<const Vertex*>& subVertices,
                            const std::map<int, int>& oldToNewVertexMap,
                            const TetrahedronRefiner<T>& tetRefiner)
    : kSubCellsPerCell(tetRefiner.getDivisionCount())

{
  using std::size_t;

  const size_t kInVertexCount = subVertices.size();
  const size_t kInCellCount = subElements.size();
  numSubCells_ = kInCellCount * kSubCellsPerCell;

  const unsigned int additionalVertices = tetRefiner.additionalVerticesPerCell();
  numVertices_ = kInVertexCount + kInCellCount * additionalVertices;

  cells_ = new unsigned int[numSubCells_ * KIndicesPerCell];
  vertices_ = new T[numVertices_ * 3];

  const std::vector<const Vertex*>& kVertices = subVertices;
  const std::vector<const Element*>& kElements = subElements;

  // Copy original vertices
#pragma omp parallel for
  for (unsigned int i = 0; i < kInVertexCount; i++) {
    memcpy(&vertices_[static_cast<size_t>(i * 3)], kVertices[i]->coords, sizeof(double) * 3);
  }

  // The pointer to the new vertices
  T* newVertices = &vertices_[kInVertexCount * 3];

  // Start the actual cell-refinement
#pragma omp parallel shared(oldToNewVertexMap)
  {
    auto* newVerticesTmp = new Eigen::Matrix<T, 3, 1>[additionalVertices];
    auto* newTetsTmp = new Tetrahedron<T>[kSubCellsPerCell];

#pragma omp for schedule(static) nowait
    for (size_t c = 0; c < kInCellCount; ++c) {
      // Build a Terahedron containing the coordinates of the vertices.
      const Tetrahedron<T> inTet =
          Tetrahedron<T>(kVertices[oldToNewVertexMap.at(kElements[c]->vertices[0])]->coords,
                         kVertices[oldToNewVertexMap.at(kElements[c]->vertices[1])]->coords,
                         kVertices[oldToNewVertexMap.at(kElements[c]->vertices[2])]->coords,
                         kVertices[oldToNewVertexMap.at(kElements[c]->vertices[3])]->coords,
                         oldToNewVertexMap.at(kElements[c]->vertices[0]),
                         oldToNewVertexMap.at(kElements[c]->vertices[1]),
                         oldToNewVertexMap.at(kElements[c]->vertices[2]),
                         oldToNewVertexMap.at(kElements[c]->vertices[3]));

      // Generate the tets
      tetRefiner.refine(inTet, kInVertexCount + c * additionalVertices, newTetsTmp, newVerticesTmp);

      // Copy new vertices
      for (unsigned int i = 0; i < additionalVertices; i++) {
        memcpy(&newVertices[(c * additionalVertices + i) * 3],
               newVerticesTmp[i].data(),
               sizeof(T) * 3);
      }

      // Copy tets
      for (unsigned int i = 0; i < kSubCellsPerCell; i++) {
        cells_[(c * kSubCellsPerCell + i) * 4] = newTetsTmp[i].i;
        cells_[(c * kSubCellsPerCell + i) * 4 + 1] = newTetsTmp[i].j;
        cells_[(c * kSubCellsPerCell + i) * 4 + 2] = newTetsTmp[i].k;
        cells_[(c * kSubCellsPerCell + i) * 4 + 3] = newTetsTmp[i].l;
      }
    }

    delete[] newVerticesTmp;
    delete[] newTetsTmp;
  }
}

template <typename T>
MeshRefiner<T>::~MeshRefiner() {
  delete[] cells_;
  delete[] vertices_;
}

//------------------------------------------------------------------------------

template <typename T>
const unsigned int* MeshRefiner<T>::getCellData() const {
  return &cells_[0];
}

template <typename T>
std::size_t MeshRefiner<T>::getkSubCellsPerCell() const {
  return kSubCellsPerCell;
}

//------------------------------------------------------------------------------

template <typename T>
const T* MeshRefiner<T>::getVertexData() const {
  return &vertices_[0];
}

//------------------------------------------------------------------------------

template <typename T>
std::size_t MeshRefiner<T>::getNumCells() const {
  return numSubCells_;
}

//------------------------------------------------------------------------------

template <typename T>
std::size_t MeshRefiner<T>::getNumVertices() const {
  return numVertices_;
}

//------------------------------------------------------------------------------

} // namespace seissol::refinement

#endif // SEISSOL_SRC_GEOMETRY_REFINEMENT_MESHREFINER_H_
