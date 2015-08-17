#ifndef VARIABLE_SUBSAMPLER_H_
#define VARIABLE_SUBSAMPLER_H_

#ifdef _OPENMP
#include "omp.h"
#endif // _OPENMP

#include <cassert>
#include <algorithm>

#include "Geometry/MeshReader.h"
#include "Numerical_aux/BasicFunction.h"
#include "RefinerUtils.h"

namespace refinement
{

//------------------------------------------------------------------------------

template<class T>
class VariableSubsampler
{
private:
    std::vector<BasicFunction::SampledBasicFunctions<T> > m_BasicFunctions;

    std::vector<Vec3D<T> > m_CenterPoints;

    const unsigned int kSubCellsPerCell;
    const unsigned int kNumVariables;
    const unsigned int kNumAlignedDOF;


    std::size_t getInVarOffset(unsigned int cell, unsigned int variable) const {
        return cell*kNumVariables*kNumAlignedDOF + kNumAlignedDOF*variable;
    }
    std::size_t getOutVarOffset(unsigned cell, unsigned int subcell,
            unsigned int variable, unsigned int selectedvars) const {
        return selectedvars * (kSubCellsPerCell * cell + subcell) + variable;
    }

public:
    VariableSubsampler(
            const TetrahedronRefiner<T>& tetRefiner,
            unsigned int order,
            unsigned int numVariables,
            unsigned int numAlignedDOF
            );

    void getSingle(const double* inData,  const unsigned int* cellMap,
            int variable, std::size_t numCells, double* outData) const;

    void getSelection(const double* inData, const unsigned int* cellMap,
            const std::vector<bool>& selection, std::size_t numCells,
            double* outData) const;
};

//------------------------------------------------------------------------------

template<typename T>
VariableSubsampler<T>::VariableSubsampler(
        const TetrahedronRefiner<T>& tetRefiner,
        unsigned int order,
        unsigned int numVariables,
        unsigned int numAlignedDOF
        ) : kSubCellsPerCell(tetRefiner.getDivisionCount()),
    kNumVariables(numVariables), kNumAlignedDOF(numAlignedDOF)
{
    // Generate cell centerpoints in the reference or unit tetrahedron.
    std::vector<Tetrahedron<T> > subCells(kSubCellsPerCell);
    tetRefiner(Tetrahedron<T>::unitTetrahedron(), subCells.data());
    m_CenterPoints.resize(kSubCellsPerCell);
    for (unsigned int t = 0; t < kSubCellsPerCell; t++)
        m_CenterPoints[t] = subCells[t].center();

    // Generate sampled basicfunctions
    for (unsigned int i = 0; i < kSubCellsPerCell; i++) {
        const Vec3D<T>& pnt = m_CenterPoints[i];
        m_BasicFunctions.push_back(
                BasicFunction::SampledBasicFunctions<T>::sampleAt(
                    order, pnt.x, pnt.y, pnt.z));
    }
}

//------------------------------------------------------------------------------

template<typename T>
void VariableSubsampler<T>::getSingle(const double* inData,  const unsigned int* cellMap,
        int variable, std::size_t numCells, double* outData) const
{
    std::vector<bool> selected(kNumVariables, false);
    selected[variable] = true;
    getSelection(inData, cellMap, selected, numCells, outData);
}

//------------------------------------------------------------------------------

template<typename T>
void VariableSubsampler<T>::getSelection(const double* inData, const unsigned int* cellMap,
        const std::vector<bool>& selection, std::size_t numCells, double* outData) const
{
    using std::size_t;
    // Enforce input assumptions
    assert(kNumVariables <= 255); // We use char to compute prefix
    assert(inData != NULL);
    assert(outData != NULL);
    assert(numCells > 0);
    assert(selection.size() == kNumVariables);

    // Calculate prefix sum of variables this denotes their position in the
    // output
    std::vector<char> prefix(kNumVariables);
    char numSelVars = 0;
    for (size_t i = 0; i < selection.size(); ++i) {
        prefix[i] = numSelVars;
        if (selection[i])
            ++numSelVars;
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    // Iterate over original Cells
    for (unsigned int c = 0; c < numCells; ++c) {
        // Process variables in selection
        for (unsigned int var = 0; var < kNumVariables; ++var) {
            // Skip unselected variables
            if (!selection[var]) continue;
            for (unsigned int sc = 0; sc < kSubCellsPerCell; ++sc) {
                outData[getOutVarOffset(c, sc, prefix[var], numSelVars)] = 
                m_BasicFunctions[sc].evalWithCoefs(
                        &inData[getInVarOffset(c, var)]
                        );
            }
        }
    }
}

//------------------------------------------------------------------------------

} // namespace

#endif // VARIABLE_SUBSAMPLER_H_
