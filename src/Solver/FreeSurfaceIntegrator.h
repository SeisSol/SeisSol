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

#include <memory>

#include <Geometry/MeshReader.h>
#include <Geometry/refinement/TriangleRefiner.h>
#include <Kernels/precision.hpp>
#include <Kernels/equations.hpp>
#include <Initializer/LTS.h>
#include <Initializer/tree/LTSTree.hpp>
#include <Initializer/tree/Lut.hpp>

#define FREESURFACE_MAX_REFINEMENT 3
#define FREESURFACE_NUMBER_OF_COMPONENTS 3

namespace seissol::solver { class FreeSurfaceIntegrator; }

class seissol::solver::FreeSurfaceIntegrator {
private:
  enum class LocationFlag {
    Elastic = 0,
    Acoustic = 1,
    FreeSurface = 2,
    FreeSurfaceWithGravity = 3
  };
  struct SurfaceLTS {
    seissol::initializers::Variable<real*> dofs;
    seissol::initializers::Variable<real*> displacementDofs;
    seissol::initializers::Variable<unsigned> side;
    seissol::initializers::Variable<unsigned> meshId;
    seissol::initializers::Variable<CellBoundaryMapping*> boundaryMapping;

    void addTo(seissol::initializers::LTSTree& surfaceLtsTree);
  };

  std::unique_ptr<real> projectionMatrixMemory;
  real* projectionMatrix[4];
  std::unique_ptr<real> projectionMatrixFromFace;
  unsigned numberOfSubTriangles;
  unsigned numberOfAlignedSubTriangles;

  static constexpr auto polyDegree = CONVERGENCE_ORDER-1;
  static constexpr auto numQuadraturePoints = polyDegree*polyDegree;
  bool m_enabled;
  
  void initializeProjectionMatrices(unsigned maxRefinementDepth);
  void computeSubTriangleAverages(real* projectionMatrixRow,
                                  const std::array<std::array<double, 3>,numQuadraturePoints>& bfPoints,
                                  double const* weights) const;
  void computeSubTriangleAveragesFromFaces(real* projectionMatrixFromFaceRow,
                                           const std::array<std::array<double, 2>,numQuadraturePoints>& bfPoints,
                                           double const* weights) const;
  void initializeSurfaceLTSTree(  seissol::initializers::LTS* lts,
                                  seissol::initializers::LTSTree* ltsTree,
                                  seissol::initializers::Lut* ltsLut );

  static LocationFlag getLocationFlag(CellMaterialData materialData, FaceType faceType, unsigned face);
public:
  real* velocities[FREESURFACE_NUMBER_OF_COMPONENTS];
  real* displacements[FREESURFACE_NUMBER_OF_COMPONENTS];

public:
  std::vector<double> locationFlags;
  unsigned totalNumberOfFreeSurfaces;
  unsigned totalNumberOfTriangles;

  SurfaceLTS surfaceLts;
  seissol::initializers::LTSTree surfaceLtsTree;
  seissol::refinement::TriangleRefiner triRefiner;
  
  explicit FreeSurfaceIntegrator();
  ~FreeSurfaceIntegrator();
  
  void initialize(  unsigned maxRefinementDepth,
                    GlobalData* globalData,
                    seissol::initializers::LTS* lts,
                    seissol::initializers::LTSTree* ltsTree,
                    seissol::initializers::Lut* ltsLut );

  void calculateOutput();
  
  bool enabled() const { return m_enabled; }
};

#endif // FREE_SURFACE_INTEGRATOR_H
