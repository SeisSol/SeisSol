#ifndef MESH_REFINER_H_
#define MESH_REFINER_H_

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include "Geometry/MeshReader.h"
#include "RefinerUtils.h"

namespace refinement
{

//------------------------------------------------------------------------------

template<typename T>
class MeshRefiner
{
private:
    // m_cells contains the indices of the cells
    std::vector<unsigned int> m_cells;
    std::vector<double> m_vertices;

    static const unsigned int kCoordsPerCell = 4 * 3;
    static const unsigned int kIndicesPerCell = 4;

    const unsigned int kSubCellsPerCell;

public:
    MeshRefiner(const MeshReader& meshReader,
            const TetrahedronRefiner<T>& tetRefiner);

    const unsigned int* getCellData() const;
    const T* getVertexData() const;

    std::size_t getNumCells() const;
    std::size_t getNumVertices() const;
};

//------------------------------------------------------------------------------

template<typename T>
MeshRefiner<T>::MeshRefiner(
        const MeshReader& meshReader,
        const TetrahedronRefiner<T>& tetRefiner
        ) : kSubCellsPerCell(tetRefiner.getDivisionCount())
{
    using std::size_t;

    const size_t kInCellCount = meshReader.getElements().size();
    const size_t kSubCellCount = kInCellCount * kSubCellsPerCell;

    m_cells.resize(kSubCellCount * kIndicesPerCell);
    m_vertices.resize(kSubCellCount * kCoordsPerCell);

    const std::vector<Vertex>& kVertices = meshReader.getVertices();
    const std::vector<Element>& kElements = meshReader.getElements();

    // Start the actual cell-refinement
    std::vector<Tetrahedron<T> > tetBuffer(kSubCellsPerCell);
#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp for schedule(static) firstprivate(tetBuffer) nowait
#endif
        for (size_t c = 0; c < kInCellCount; ++c)
        {
            const size_t kCellIndexInVertexMap = c * kSubCellsPerCell * kCoordsPerCell;

            // Build a Terahedron containing the coordinates of the vertices.
            Tetrahedron<T> inTet = Tetrahedron<T>(
                    kVertices[kElements[c].vertices[0]].coords,
                    kVertices[kElements[c].vertices[1]].coords,
                    kVertices[kElements[c].vertices[2]].coords,
                    kVertices[kElements[c].vertices[3]].coords);
            T* out = m_vertices.data() + kCellIndexInVertexMap;
            tetRefiner(inTet, tetBuffer.data());
            for (int i = 0; i < kSubCellsPerCell; ++i)
                tetBuffer[i].dumpData(out);
        }

        // Set the vertex indices (vid) of each tetrahedron.
        // Each tetrahedron has 4 indices and produces exactly 4 vertices.
        // Bijectivity is given and the order is well known.
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (unsigned int vid  = 0; vid < m_cells.size(); ++vid) {
            m_cells[vid] = vid*3;
        }
#ifdef _OPENMP
    }
#endif
};

//------------------------------------------------------------------------------

template<typename T>
const unsigned int* MeshRefiner<T>::getCellData() const {
    return m_cells.data();
}

//------------------------------------------------------------------------------

template<typename T>
const T* MeshRefiner<T>::getVertexData() const {
    return m_vertices.data();
}

//------------------------------------------------------------------------------

template<typename T>
std::size_t MeshRefiner<T>::getNumCells() const {
    return m_cells.size() / 4;
}

//------------------------------------------------------------------------------

template<typename T>
std::size_t MeshRefiner<T>::getNumVertices() const {
    return m_vertices.size() / 3;
}

//------------------------------------------------------------------------------

} // namespace

#endif // MESH_REFINER_H_
