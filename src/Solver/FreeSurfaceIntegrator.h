/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 */

#ifndef FREE_SURFACE_INTEGRATOR_H
#define FREE_SURFACE_INTEGRATOR_H

#include <Geometry/MeshReader.h>
#include <Geometry/refinement/TriangleRefiner.h>
#include <Kernels/precision.hpp>
#include <Kernels/equations.hpp>
#include <Initializer/LTS.h>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Lut.hpp>

#define FREESURFACE_MAX_REFINEMENT 3
#define FREESURFACE_NUMBER_OF_COMPONENTS 3

namespace seissol
{
  namespace solver
  {
    class FreeSurfaceIntegrator;
  }
}

class seissol::solver::FreeSurfaceIntegrator {
private:
  struct SurfaceLTS {
    seissol::initializers::Variable<real*> dofs;
    seissol::initializers::Variable<real*> displacementDofs;
    seissol::initializers::Variable<unsigned> side;
    seissol::initializers::Variable<unsigned> meshId;
    
    void addTo(seissol::initializers::LTSTree& surfaceLtsTree);
  };

  real* projectionMatrixMemory;
  real* projectionMatrix[4];
  unsigned numberOfSubTriangles;
  unsigned numberOfAlignedSubTriangles;
  
  bool m_enabled;
  
  void initializeProjectionMatrices(unsigned maxRefinementDepth);
  void computeSubTriangleAverages(  real* projectionMatrixRow,
                                    double const (*bfPoints)[3],
                                    double const* weights,
                                    unsigned numQuadraturePoints  );
  void initializeSurfaceLTSTree(  seissol::initializers::LTS* lts,
                                  seissol::initializers::LTSTree* ltsTree,
                                  seissol::initializers::Lut* ltsLut );
  
public:  
  real* velocities[FREESURFACE_NUMBER_OF_COMPONENTS];
  real* displacements[FREESURFACE_NUMBER_OF_COMPONENTS];
  unsigned totalNumberOfFreeSurfaces;
  unsigned totalNumberOfTriangles;

  SurfaceLTS surfaceLts;
  seissol::initializers::LTSTree surfaceLtsTree;
  seissol::refinement::TriangleRefiner triRefiner;
  
  explicit FreeSurfaceIntegrator();
  ~FreeSurfaceIntegrator();
  
  void initialize(  unsigned maxRefinementDepth,
                    seissol::initializers::LTS* lts,
                    seissol::initializers::LTSTree* ltsTree,
                    seissol::initializers::Lut* ltsLut );

  void calculateOutput();
  
  bool enabled() const { return m_enabled; }
};

#endif // FREE_SURFACE_INTEGRATOR_H
