// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_GEOMETRY_REFINEMENT_MESHREFINER_H_
#define SEISSOL_SRC_GEOMETRY_REFINEMENT_MESHREFINER_H_

#include <cstring>

#include "Geometry/MeshReader.h"
#include "RefinerUtils.h"

namespace seissol {
namespace refinement {

//------------------------------------------------------------------------------

template <typename T>
class MeshRefiner {
  private:
  // m_cells contains the indices of the cells
  unsigned int* m_cells;
  T* m_vertices;

  size_t m_numSubCells;
  size_t m_numVertices;

  static const unsigned int kIndicesPerCell = 4;

  const unsigned int kSubCellsPerCell;

  public:
  MeshRefiner(const seissol::geometry::MeshReader& meshReader,
              const TetrahedronRefiner<T>& tetRefiner);

  MeshRefiner(const std::vector<const Element*>& subElements,
              const std::vector<const Vertex*>& subVertices,
              const std::map<int, int>& oldToNewVertexMap,
              const TetrahedronRefiner<T>& tetRefiner);

  ~MeshRefiner();

  const unsigned int* getCellData() const;
  const T* getVertexData() const;
  std::size_t getkSubCellsPerCell() const;

  std::size_t getNumCells() const;
  std::size_t getNumVertices() const;
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
  m_numSubCells = kInCellCount * kSubCellsPerCell;

  const unsigned int additionalVertices = tetRefiner.additionalVerticesPerCell();
  m_numVertices = kInVertexCount + kInCellCount * additionalVertices;

  m_cells = new unsigned int[m_numSubCells * kIndicesPerCell];
  m_vertices = new T[m_numVertices * 3];

  const std::vector<Vertex>& kVertices = meshReader.getVertices();
  const std::vector<Element>& kElements = meshReader.getElements();

  // Copy original vertices
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
  for (unsigned int i = 0; i < kInVertexCount; i++) {
    memcpy(&m_vertices[i * 3], kVertices[i].coords, sizeof(double) * 3);
  }

  // The pointer to the new vertices
  T* newVertices = &m_vertices[kInVertexCount * 3];

  // Start the actual cell-refinement
#ifdef _OPENMP
#pragma omp parallel
  {
#endif // _OPENMPI
    Eigen::Matrix<T, 3, 1>* newVerticesTmp = new Eigen::Matrix<T, 3, 1>[additionalVertices];
    Tetrahedron<T>* newTetsTmp = new Tetrahedron<T>[kSubCellsPerCell];

#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif // _OPENMP
    for (size_t c = 0; c < kInCellCount; ++c) {
      // Build a Terahedron containing the coordinates of the vertices.
      Tetrahedron<T> inTet = Tetrahedron<T>(kVertices[kElements[c].vertices[0]].coords,
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
        m_cells[(c * kSubCellsPerCell + i) * 4] = newTetsTmp[i].i;
        m_cells[(c * kSubCellsPerCell + i) * 4 + 1] = newTetsTmp[i].j;
        m_cells[(c * kSubCellsPerCell + i) * 4 + 2] = newTetsTmp[i].k;
        m_cells[(c * kSubCellsPerCell + i) * 4 + 3] = newTetsTmp[i].l;
      }
    }

    delete[] newVerticesTmp;
    delete[] newTetsTmp;
#ifdef _OPENMP
  }
#endif
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
  m_numSubCells = kInCellCount * kSubCellsPerCell;

  const unsigned int additionalVertices = tetRefiner.additionalVerticesPerCell();
  m_numVertices = kInVertexCount + kInCellCount * additionalVertices;

  m_cells = new unsigned int[m_numSubCells * kIndicesPerCell];
  m_vertices = new T[m_numVertices * 3];

  const std::vector<const Vertex*>& kVertices = subVertices;
  const std::vector<const Element*>& kElements = subElements;

  // Copy original vertices
#ifdef _OPENMP
#pragma omp parallel for
#endif // _OPENMP
  for (unsigned int i = 0; i < kInVertexCount; i++) {
    memcpy(&m_vertices[i * 3], kVertices[i]->coords, sizeof(double) * 3);
  }

  // The pointer to the new vertices
  T* newVertices = &m_vertices[kInVertexCount * 3];

  // Start the actual cell-refinement
#ifdef _OPENMP
#pragma omp parallel shared(oldToNewVertexMap)
  {
#endif // _OPENMPI
    Eigen::Matrix<T, 3, 1>* newVerticesTmp = new Eigen::Matrix<T, 3, 1>[additionalVertices];
    Tetrahedron<T>* newTetsTmp = new Tetrahedron<T>[kSubCellsPerCell];

#ifdef _OPENMP
#pragma omp for schedule(static) nowait
#endif // _OPENMP
    for (size_t c = 0; c < kInCellCount; ++c) {
      // Build a Terahedron containing the coordinates of the vertices.
      Tetrahedron<T> inTet =
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
        m_cells[(c * kSubCellsPerCell + i) * 4] = newTetsTmp[i].i;
        m_cells[(c * kSubCellsPerCell + i) * 4 + 1] = newTetsTmp[i].j;
        m_cells[(c * kSubCellsPerCell + i) * 4 + 2] = newTetsTmp[i].k;
        m_cells[(c * kSubCellsPerCell + i) * 4 + 3] = newTetsTmp[i].l;
      }
    }

    delete[] newVerticesTmp;
    delete[] newTetsTmp;
#ifdef _OPENMP
  }
#endif
}

template <typename T>
MeshRefiner<T>::~MeshRefiner() {
  delete[] m_cells;
  delete[] m_vertices;
}

//------------------------------------------------------------------------------

template <typename T>
const unsigned int* MeshRefiner<T>::getCellData() const {
  return &m_cells[0];
}

template <typename T>
std::size_t MeshRefiner<T>::getkSubCellsPerCell() const {
  return kSubCellsPerCell;
}

//------------------------------------------------------------------------------

template <typename T>
const T* MeshRefiner<T>::getVertexData() const {
  return &m_vertices[0];
}

//------------------------------------------------------------------------------

template <typename T>
std::size_t MeshRefiner<T>::getNumCells() const {
  return m_numSubCells;
}

//------------------------------------------------------------------------------

template <typename T>
std::size_t MeshRefiner<T>::getNumVertices() const {
  return m_numVertices;
}

//------------------------------------------------------------------------------

} // namespace refinement
} // namespace seissol

#endif // SEISSOL_SRC_GEOMETRY_REFINEMENT_MESHREFINER_H_
