// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_GEOMETRY_SAMPLER_H_
#define SEISSOL_SRC_IO_INSTANCE_GEOMETRY_SAMPLER_H_

#include <functional>
#include <vector>
namespace seissol::io::instance::geometry {

template <std::size_t Dim>
class SamplerBuilder {
  public:
  SamplerBuilder() {
    shapes.resize(1);
    shapes[0].resize(Dim + 1);

    // shapes[0][0] is implicitly all zero

    for (std::size_t i = 0; i < Dim; ++i) {
      shapes[0][i + 1][i] = 1;
    }
  }

  void project();

  SamplerBuilder& subdivide(std::vector<std::vector<std::array<double, Dim>>>& divide) {
    std::vector<std::vector<std::array<double, Dim>>> shapesNew;
    for (const auto& shape : shapes) {
      for (const auto& section : divide) {
      }
    }
    return *this;
  }

  SamplerBuilder& sample(int order) { return *this; }

  const std::function<void(void*, const void*)>& getSampler();
  const std::function<void> getPointMapper();

  private:
  std::vector<std::vector<std::array<double, Dim>>> shapes;
  std::vector<std::array<double, Dim>> points;
};

} // namespace seissol::io::instance::geometry

#endif // SEISSOL_SRC_IO_INSTANCE_GEOMETRY_SAMPLER_H_
