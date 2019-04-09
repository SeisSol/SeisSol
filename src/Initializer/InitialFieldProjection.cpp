/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2019, SeisSol Group
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
 * 
 **/

#include "InitialFieldProjection.h"

#include <Numerical_aux/Quadrature.h>
#include <Numerical_aux/BasisFunction.h>
#include <Numerical_aux/Transformation.h>
#include <generated_code/kernels.h>

void seissol::initializers::projectInitialField(  physics::InitialField const&  iniField,
                                                  GlobalData const&             globalData,
                                                  MeshReader const&             meshReader,                                                    
                                                  LTS const&                    lts,
                                                  Lut const&                    ltsLut )
{
  auto const& vertices = meshReader.getVertices();
  auto const& elements = meshReader.getElements();

  constexpr auto quadPolyDegree = CONVERGENCE_ORDER+1;
  constexpr auto numQuadPoints = quadPolyDegree * quadPolyDegree * quadPolyDegree;

  double quadraturePoints[numQuadPoints][3];
  double quadratureWeights[numQuadPoints];
  seissol::quadrature::TetrahedronQuadrature(quadraturePoints, quadratureWeights, quadPolyDegree);

  assert(model::dofsQP::rows == numQuadPoints);

#ifdef _OPENMP
  #pragma omp parallel
  {
#endif
  real dofsQPdata[model::dofsQP::reals] __attribute__((aligned(ALIGNMENT)));
  MatrixView dofsQP(dofsQPdata, sizeof(dofsQPdata)/sizeof(real), &colMjrIndex<model::dofsQP::ld>);

  std::vector<std::array<double, 3>> quadraturePointsXyz;
  quadraturePointsXyz.resize(numQuadPoints);

#ifdef _OPENMP
  #pragma omp for schedule(static)
#endif
  for (unsigned int meshId = 0; meshId < elements.size(); ++meshId) {
    double const* elementCoords[4];
    for (size_t v = 0; v < 4; ++v) {
      elementCoords[v] = vertices[elements[meshId].vertices[ v ] ].coords;
    }
    for (size_t i = 0; i < numQuadPoints; ++i) {
      seissol::transformations::tetrahedronReferenceToGlobal(elementCoords[0], elementCoords[1], elementCoords[2], elementCoords[3], quadraturePoints[i], quadraturePointsXyz[i].data());
    }

    iniField.evaluate(0.0, quadraturePointsXyz, dofsQP);
    generatedKernels::projectQP(  dofsQPdata,
                                  globalData.projectQPMatrix,
                                  ltsLut.lookup(lts.dofs, meshId) );
  }
#ifdef _OPENMP
  }
#endif
}
