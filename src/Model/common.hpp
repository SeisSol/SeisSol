/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Wolf (wolf.sebastian AT in.tum.de,
 *https://www5.in.tum.de/wiki/index.php/Sebastian_Wolf,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2015 - 2020, SeisSol Group
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
 **/

#ifndef MODEL_COMMON_HPP_
#define MODEL_COMMON_HPP_

#include <Eigen/Dense>

#include "utils/logger.h"
#include "Initializer/typedefs.hpp"
#include "generated_code/init.h"
#include "Geometry/MeshTools.h"
#include "Numerical_aux/Transformation.h"

#include "Model/common_datastructures.hpp"

// for now
#include "Equations/viscoelastic/Model/datastructures.hpp"
#include "Equations/poroelastic/Model/datastructures.hpp"

namespace seissol::model {
template <typename T, typename Tmatrix, typename Tarray1, typename Tarray2>
void setBlocks(T QgodLocal, Tmatrix S, Tarray1 traction_indices, Tarray2 velocity_indices) {
  // set lower left block
  int col = 0;
  for (auto& t : traction_indices) {
    int row = 0;
    for (auto& v : velocity_indices) {
      QgodLocal(t, v) = S(row, col);
      row++;
    }
    col++;
  }

  // set lower right block
  for (auto& v : velocity_indices) {
    QgodLocal(v, v) = 1.0;
  }
}

bool testIfAcoustic(real mu);

template <typename Tmaterial, typename Tmatrix>
void getTransposedCoefficientMatrix(Tmaterial const& i_material, unsigned i_dim, Tmatrix& o_M) {
  o_M.setZero();
}

template <typename Tmaterial, typename T>
void getTransposedSourceCoefficientTensor(Tmaterial const& material, T& E) {}

template <typename Tmaterial, typename Tloc, typename Tneigh>
void getTransposedGodunovState(Tmaterial const& local,
                               Tmaterial const& neighbor,
                               FaceType faceType,
                               Tloc& QgodLocal,
                               Tneigh& QgodNeighbor);

template <typename Tmaterial, typename T, typename Tmatrix>
void getTransposedFreeSurfaceGodunovState(Tmaterial const& local,
                                          T& QgodLocal,
                                          T& QgodNeighbor,
                                          Tmatrix& R);

template <typename T>
void getPlaneWaveOperator(
    T const& material,
    double const n[3],
    std::complex<double> Mdata[T::NumberOfQuantities * T::NumberOfQuantities]);

template <typename Config>
void initializeSpecificLocalData(const typename Config::MaterialT&,
                                 double timeStepWidth,
                                 LocalSpecificData<Config>* LocalData) {}

template <typename Config>
void initializeSpecificNeighborData(const typename Config::MaterialT&,
                                    NeighborSpecificData<Config>* NeighborData) {}

/*
 * Calculates the so called Bond matrix. Anisotropic materials are characterized by
 * 21 different material parameters. Due to the directional dependence of anisotropic
 * materials the parameters are not independet of the choice of the coordinate system.
 * The Bond matrix transforms materials from one orthogonal coordinate system to
 * another one.
 * c.f. 10.1111/j.1365-246X.2007.03381.x
 */
void getBondMatrix(VrtxCoords const i_normal,
                   VrtxCoords const i_tangent1,
                   VrtxCoords const i_tangent2,
                   real* o_N);

template <typename MaterialT>
void getFaceRotationMatrix(VrtxCoords const i_normal,
                           VrtxCoords const i_tangent1,
                           VrtxCoords const i_tangent2,
                           init::T::view::type& o_T,
                           init::Tinv::view::type& o_Tinv);

template <typename MaterialT>
void getFaceRotationMatrix(Eigen::Vector3d const i_normal,
                           Eigen::Vector3d const i_tangent1,
                           Eigen::Vector3d const i_tangent2,
                           init::T::view::type& o_T,
                           init::Tinv::view::type& o_Tinv) {
  VrtxCoords n = {i_normal(0), i_normal(1), i_normal(2)};
  VrtxCoords s = {i_tangent1(0), i_tangent1(1), i_tangent1(2)};
  VrtxCoords t = {i_tangent2(0), i_tangent2(1), i_tangent2(2)};
  getFaceRotationMatrix<MaterialT>(n, s, t, o_T, o_Tinv);
}

template <typename MaterialT, typename RealT>
MaterialT getRotatedMaterialCoefficients(RealT rotationParameters[36], MaterialT&& material) {
  return std::forward<MaterialT>(material);
}
} // namespace seissol::model

template <typename T>
void seissol::model::getPlaneWaveOperator(
    T const& material,
    double const n[3],
    std::complex<double> Mdata[T::NumberOfQuantities * T::NumberOfQuantities]) {
  yateto::DenseTensorView<2, std::complex<double>> M(
      Mdata, {T::NumberOfQuantities, T::NumberOfQuantities});
  M.setZero();

  double data[T::NumberOfQuantities * T::NumberOfQuantities];
  yateto::DenseTensorView<2, double> Coeff(data, {T::NumberOfQuantities, T::NumberOfQuantities});

  for (unsigned d = 0; d < 3; ++d) {
    Coeff.setZero();
    getTransposedCoefficientMatrix(material, d, Coeff);

    for (unsigned i = 0; i < T::NumberOfQuantities; ++i) {
      for (unsigned j = 0; j < T::NumberOfQuantities; ++j) {
        M(i, j) += n[d] * Coeff(j, i);
      }
    }
  }
  Coeff.setZero();
  getTransposedSourceCoefficientTensor(material, Coeff);

  for (unsigned i = 0; i < T::NumberOfQuantities; ++i) {
    for (unsigned j = 0; j < T::NumberOfQuantities; ++j) {
      M(i, j) -= std::complex<real>(0.0, Coeff(j, i));
    }
  }
}

#endif
