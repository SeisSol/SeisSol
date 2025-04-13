#pragma once

#include <Eigen/Dense>
#include <Kernels/Precision.h>
#include <init.h>

namespace seissol::initializer {

struct GlobalMatrixPointers {
    real* kDivMT;
    real* kDivM;
    real* rDivM;

    real ew_extraweight_M[seissol::init::ew_extraweight_M::Size];
    real ew_extraweight_k[seissol::init::ew_extraweight_k::Size];
    real ew_extraweight_r[seissol::init::ew_extraweight_r::Size];

    template<typename TransformT>
    void sampleBasis(const TransformT& transform) {
        auto matM = seissol::init::ew_extraweight_M::view::create(ew_extraweight_M);
        auto matk = seissol::init::ew_extraweight_k::view::create(ew_extraweight_k);
        auto matr = seissol::init::ew_extraweight_r::view::create(ew_extraweight_r);
        auto volumePoints = seissol::init::ew_quad_nodes_vv::view::create(const_cast<real*>(seissol::init::ew_quad_nodes_vv::Values));
        for (int i = 0; i < seissol::init::ew_quad_nodes_vv::Shape[0]; ++i) {
            Eigen::Vector3d point(volumePoints(i, 0), volumePoints(i, 1), volumePoints(i, 2));
            Eigen::Matrix3d dvol = transform.mapDVolume(point);
            const auto dvoldet = dvol.determinant();
            matM(i) = dvoldet;

            // reminder: inverse function derivative:
            // f^{-1}'(y) = 1 / (f'(f^{-1}(y)))
            // thus: dphi(map^{-1})/dx = Jphi @ dmap^{-1}/dx = Jphi @ (1 / (dmap/dx(map^{-1})))
            // after the variable transformation, the map^{-1} cancels out.
            // (may need some double checking)
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    // TODO: check ordering (transpose?)
                    matk(i, j, k) = dvoldet / dvol(j, k);
                }
            }
        }

        auto facePoints = seissol::init::ew_quad_nodes_ff::view::create(const_cast<real*>(seissol::init::ew_quad_nodes_ff::Values));
        for (int f = 0; f < 4; ++f) {
            for (int i = 0; i < seissol::init::ew_quad_nodes_ff::Shape[0]; ++i) {
                Eigen::Vector2d point(facePoints(i, 0), facePoints(i, 1));
                Eigen::Vector3d dface = transform.faceDirection(f, point);
                matr(i, f) = dface.norm();
            }
        }
    }

    void bootstrapMatrices();
};

} // namespace seissol::initializer
