#pragma once

#include <Eigen/Dense>
#include <Kernels/Precision.h>
#include <Numerical/Transformation.h>
namespace seissol::initializer {

class AffineTransform {
public:
    AffineTransform(real x[4], real y[4], real z[4]) {
        Eigen::Vector3d v1(x[0], y[0], z[0]);
        Eigen::Vector3d v2(x[1], y[1], z[1]);
        Eigen::Vector3d v3(x[2], y[2], z[2]);
        Eigen::Vector3d v4(x[3], y[3], z[3]);

        Eigen::Vector3d v21 = v2 - v1;
        Eigen::Vector3d v31 = v3 - v1;
        Eigen::Vector3d v41 = v4 - v1;

        Eigen::Matrix3d imat;
        imat(0, 0) = v21(0);
        imat(1, 0) = v21(1);
        imat(2, 0) = v21(2);
        imat(0, 1) = v31(0);
        imat(1, 1) = v31(1);
        imat(2, 1) = v31(2);
        imat(0, 2) = v41(0);
        imat(1, 2) = v41(1);
        imat(2, 2) = v41(2);

        offset = v1;
        transform = imat.inverse();
    }

    [[nodiscard]] Eigen::Matrix3d mapDVolume(const Eigen::Vector3d& point) const {
        return transform;
    }

    [[nodiscard]] Eigen::Vector3d mapVolume(const Eigen::Vector3d& point) const {
        return transform * point + offset;
    }

    [[nodiscard]] Eigen::Vector3d faceDirection(int face, const Eigen::Vector2d& point) const {
        return transform * baseFaceDirection(face);
    }

    static Eigen::Vector3d baseFaceDirection(int face) {
        if (face == 0) {
            return Eigen::Vector3d(1, 0, 0);
        }
        else if (face == 1) {
            return Eigen::Vector3d(0, 1, 0);
        }
        else if (face == 2) {
            return Eigen::Vector3d(0, 0, 1);
        }
        else if (face == 3) {
            return Eigen::Vector3d(1 / std::sqrt(3), 1 / std::sqrt(3), 1 / std::sqrt(3));
        }
        else {
            return Eigen::Vector3d(0,0,0);
        }
    }

private:
    Eigen::Matrix3d transform;
    Eigen::Vector3d offset;
};

} // namespace seissol::initializer
