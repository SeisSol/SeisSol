#include "MatrixBootstrap.h"

#include <Eigen/Dense>
#include <kernel.h>
#include <yateto/TensorView.h>

namespace {

template<typename RealT>
Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic> toEigen(const yateto::DenseTensorView<2, RealT>& tensorView) {
    Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic> eigenMatrix;
    eigenMatrix.resize(tensorView.shape(0), tensorView.shape(1));
    for (int i = 0; i < tensorView.shape(0); ++i) {
        for (int j = 0; j < tensorView.shape(1); ++j) {
            eigenMatrix(i,j) = tensorView(i,j);
        }
    }
    return eigenMatrix;
}

template<typename RealT>
Eigen::FullPivLU<Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic>> invert(const Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic>& matrix) {
    return Eigen::FullPivLU<Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic>>(matrix);
}

// compute matC = inv(matA) @ matB
template<typename ResT, typename LeftT, typename RealT, typename IntT, unsigned N>
void invertAgainst(real* matC, const LeftT& matA, yateto::DenseTensorView<N, RealT, IntT>& matB) {
    if constexpr (N == 3) {
        for (int i = 0; i < matB.shape(2); ++i) {
            auto matBView = matB.subtensor(yateto::slice<>(), yateto::slice<>(), i);
            invertAgainst<ResT>(matC, matA, matBView);
            matC += ResT::size(i);
        }
    }
    else {
        auto resview = ResT::template view<0>::create(matC);
        Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic> result = matA.solve(toEigen(matB));
        for (int i = 0; i < resview.shape(0); ++i) {
            for (int j = 0; j < resview.shape(1); ++j) {
                resview(i, j) = result(i, j);
            }
        }
    }
}

} // namespace

namespace seissol::initializer {

void GlobalMatrixPointers::bootstrapMatrices() {
    real dataM[seissol::init::ew_M::Size];
    real datak[seissol::init::ew_k::Size];
    real datakT[seissol::init::ew_kT::Size];
    real datar[seissol::init::ew_r::Size];
    auto matM = seissol::init::ew_M::view::create(dataM);
    auto matk = seissol::init::ew_k::view::create(datak);
    auto matkT = seissol::init::ew_kT::view::create(datakT);
    auto matr = seissol::init::ew_r::view::create(datar);
    seissol::kernel::bootstrap bootstrap;

    bootstrap.ew_collocate_f_vv = seissol::init::ew_collocate_f_vv::Values;
    bootstrap.ew_collocate_df_vv = seissol::init::ew_collocate_df_vv::Values;
    bootstrap.ew_collocate_f_vf = seissol::init::ew_collocate_f_vf::Values;
    bootstrap.ew_collocate_f_ff = seissol::init::ew_collocate_f_ff::Values;
    bootstrap.ew_quad_weights_v = seissol::init::ew_quad_weights_v::Values;
    bootstrap.ew_quad_weights_f = seissol::init::ew_quad_weights_f::Values;

    bootstrap.ew_extraweight_M = ew_extraweight_M;
    bootstrap.ew_extraweight_k = ew_extraweight_k;
    bootstrap.ew_extraweight_r = ew_extraweight_r;

    bootstrap.ew_M = dataM;
    bootstrap.ew_k = datak;
    bootstrap.ew_kT = datakT;
    bootstrap.ew_r = datar;
    
    bootstrap.execute();
    auto divM = invert(toEigen(matM));
    invertAgainst<seissol::init::globalMkDivM>(kDivM, divM, matk);
    invertAgainst<seissol::init::globalMkDivMT>(kDivMT, divM, matkT);
    invertAgainst<seissol::init::globalMrDivM>(rDivM, divM, matr);
}

} // namespace seissol::initializer
