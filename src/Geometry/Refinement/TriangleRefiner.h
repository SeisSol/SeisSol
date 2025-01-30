// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_GEOMETRY_REFINEMENT_TRIANGLEREFINER_H_
#define SEISSOL_SRC_GEOMETRY_REFINEMENT_TRIANGLEREFINER_H_

#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace seissol {
namespace refinement {
struct Triangle;
class TriangleRefiner;
} // namespace refinement
} // namespace seissol

struct seissol::refinement::Triangle {
  Eigen::Vector2d x[3];
  double area;
};

class seissol::refinement::TriangleRefiner {
  private:
  void refine(Eigen::Vector2d x0,
              Eigen::Vector2d x1,
              Eigen::Vector2d x2,
              unsigned depth,
              unsigned maxDepth) {
    Eigen::Vector2d chi = x1 - x0;
    Eigen::Vector2d tau = x2 - x0;

    if (depth < maxDepth) {
      // Edge midpoints
      Eigen::Vector2d m0 = 0.5 * chi + 0.5 * tau + x0;
      Eigen::Vector2d m1 = 0.5 * tau + x0;
      Eigen::Vector2d m2 = 0.5 * chi + x0;

      refine(x0, m2, m1, depth + 1, maxDepth);
      refine(m2, x1, m0, depth + 1, maxDepth);
      refine(m1, m0, x2, depth + 1, maxDepth);
      refine(m0, m1, m2, depth + 1, maxDepth);
    } else {
      Triangle sub;
      sub.x[0] = x0;
      sub.x[1] = x1;
      sub.x[2] = x2;
      sub.area = fabs(chi(0) * tau(1) - chi(1) * tau(0));
      subTris.push_back(sub);
    }
  }

  public:
  std::vector<Triangle> subTris;
  unsigned maxDepth;

  explicit TriangleRefiner() {}

  void refine(unsigned maxRefinementDepth) {
    subTris.clear();
    maxDepth = maxRefinementDepth;
    refine(Eigen::Vector2d(0.0, 0.0),
           Eigen::Vector2d(1.0, 0.0),
           Eigen::Vector2d(0.0, 1.0),
           0,
           maxRefinementDepth);
  }
};

#endif // SEISSOL_SRC_GEOMETRY_REFINEMENT_TRIANGLEREFINER_H_
