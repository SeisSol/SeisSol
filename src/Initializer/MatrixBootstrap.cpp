// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "MatrixBootstrap.h"

#include "Common/Constants.h"

#include <Eigen/Dense>
#include <GeneratedCode/kernel.h>
#include <yateto/TensorView.h>

namespace seissol::initializer {

namespace {

template <typename RealT>
Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic>
    toEigen(const yateto::DenseTensorView<2, RealT>& tensorView) {
  Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic> eigenMatrix;
  eigenMatrix.resize(tensorView.shape(0), tensorView.shape(1));
  for (std::size_t i = 0; i < tensorView.shape(0); ++i) {
    for (std::size_t j = 0; j < tensorView.shape(1); ++j) {
      eigenMatrix(i, j) = tensorView(i, j);
    }
  }
  return eigenMatrix;
}

template <typename RealT>
Eigen::FullPivLU<Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic>>
    invert(const Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic>& matrix) {
  return Eigen::FullPivLU<Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic>>(matrix);
}

// compute matC = inv(matA) @ matB
template <typename ResT, typename LeftT, typename RealT, typename IntT, unsigned N>
void invertAgainst(real* matC,
                   const LeftT& matA,
                   yateto::DenseTensorView<N, RealT, IntT>& matB,
                   bool transpose,
                   bool antitranspose = false) {
  if constexpr (N == 3) {
    for (std::size_t i = 0; i < matB.shape(2); ++i) {
      auto matBView = matB.subtensor(yateto::slice<>(), yateto::slice<>(), i);
      invertAgainst<ResT>(matC, matA, matBView, transpose);
      matC += ResT::size(i);
    }
  } else {
    auto resview = ResT::template view<0>::create(matC);
    const auto eigenB = toEigen(matB);
    const auto solveB = transpose ? eigenB.transpose() : eigenB;
    const Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic> result = matA.solve(solveB);
    if (transpose && !antitranspose) {
      for (std::size_t i = 0; i < resview.shape(0); ++i) {
        for (std::size_t j = 0; j < resview.shape(1); ++j) {
          resview(i, j) = result(j, i);
        }
      }
    } else {
      for (std::size_t i = 0; i < resview.shape(0); ++i) {
        for (std::size_t j = 0; j < resview.shape(1); ++j) {
          resview(i, j) = result(i, j);
        }
      }
    }
  }
}

} // namespace

void GlobalMatrixPointers::bootstrapMatrices() {
  real dataM[seissol::init::ew_M::Size]{};
  real datafM[seissol::init::ew_fM::Size]{};
  real datak[seissol::init::ew_k::Size]{};
  real datakT[seissol::init::ew_kT::Size]{};
  real datar[seissol::init::ew_r::Size]{};
  real datar2[seissol::init::ew_r::Size]{};
  auto matM = seissol::init::ew_M::view::create(dataM);
  auto matfM = seissol::init::ew_fM::view::create(datafM);
  auto matk = seissol::init::ew_k::view::create(datak);
  auto matkT = seissol::init::ew_kT::view::create(datakT);
  auto matr = seissol::init::ew_r::view::create(datar);
  auto matr2 = seissol::init::ew_r::view::create(datar2);
  seissol::kernel::bootstrap bootstrap;

  bootstrap.ew_collocate_f_vv = seissol::init::ew_collocate_f_vv::Values;
  bootstrap.ew_collocate_df_vv = seissol::init::ew_collocate_df_vv::Values;
  bootstrap.ew_collocate_f_vf = seissol::init::ew_collocate_f_vf::Values;
  bootstrap.ew_collocate_f_ff = seissol::init::ew_collocate_f_ff::Values;
  bootstrap.ew_quad_weights_v = seissol::init::ew_quad_weights_v::Values;
  bootstrap.ew_quad_weights_f = seissol::init::ew_quad_weights_f::Values;

  bootstrap.ew_extraweight_M = ew_extraweight_M;
  bootstrap.ew_extraweight_fM = ew_extraweight_fM;
  bootstrap.ew_extraweight_k = ew_extraweight_k;
  bootstrap.ew_extraweight_r = ew_extraweight_r;

  bootstrap.ew_M = dataM;
  bootstrap.ew_fM = datafM;
  bootstrap.ew_k = datak;
  bootstrap.ew_kT = datakT;
  bootstrap.ew_r = datar;

  bootstrap.execute();

  auto divM = invert(toEigen(matM));
  invertAgainst<seissol::init::globalMkDivM>(kDivM, divM, matk, false);
  invertAgainst<seissol::init::globalMkDivMT>(kDivMT, divM, matkT, false);
  invertAgainst<seissol::init::globalMrDivM>(datar2, divM, matr, false);
  real* rDivMPtr = rDivM;
  real* rTPtr = rT;
  for (std::size_t f = 0; f < Cell::NumFaces; ++f) {
    auto divfM = invert(toEigen(matfM.subtensor(yateto::slice<>(), yateto::slice<>(), f)));
    auto subRDiv = matr2.subtensor(yateto::slice<>(), yateto::slice<>(), f);
    invertAgainst<seissol::init::globalMrDivM>(rDivMPtr, divfM, subRDiv, true);
    rDivMPtr += seissol::init::globalMrDivM::size(f);

    auto subR = matr.subtensor(yateto::slice<>(), yateto::slice<>(), f);
    invertAgainst<seissol::init::globalMrT>(rTPtr, divfM, subR, true, true);
    rTPtr += seissol::init::globalMrT::size(f);
  }

  // rT is reverse
  // fMrT is averse

  const auto nextfM = seissol::init::ew_fM::size() / seissol::init::ew_fM::Shape[2];
  for (std::size_t i = 0; i < seissol::init::ew_fM::size(); ++i) {
    const auto fscale = facescale[i / nextfM];
    fMrT[i] = datafM[i] / fscale;
  }

  // const auto nextR = seissol::init::ew_r::size() / seissol::init::ew_r::Shape[2];
  for (std::size_t i = 0; i < seissol::init::ew_r::size(); ++i) {
    // const auto fscale = facescale[i / nextR];
    rT[i] /= volscale;
  }
}

} // namespace seissol::initializer
