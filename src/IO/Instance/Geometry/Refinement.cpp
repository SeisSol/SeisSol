// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger
// SPDX-FileContributor: David Schneller

#include "Refinement.h"

#include <array>
#include <vector>

namespace {
using T3 = std::vector<std::array<double, 2>>;
using T4 = std::vector<std::array<double, 3>>;
using P2 = std::array<double, 2>;
using P3 = std::array<double, 3>;

const auto Va = P3{0, 0, 0};
const auto Vb = P3{0, 0, 1};
const auto Vc = P3{1, 0, 0};
const auto Vd = P3{0, 1, 0};

template <std::size_t Dim, typename... Points>
std::array<double, Dim> average(const Points&... points) {
  std::array<double, Dim> output{};
  for (std::size_t i = 0; i < Dim; ++i) {
    output[i] = (points[i] + ...) / static_cast<double>(sizeof...(Points));
  }
  return output;
}

const auto Vab = average<3>(Va, Vb);
const auto Vac = average<3>(Va, Vc);
const auto Vad = average<3>(Va, Vd);
const auto Vbc = average<3>(Vb, Vc);
const auto Vbd = average<3>(Vb, Vd);
const auto Vcd = average<3>(Vc, Vd);

const auto Vabcd = average<3>(Va, Vb, Vc, Vd);

} // namespace

namespace seissol::io::instance::geometry {

// don't care about the ordering here (kind of; it's copied from the old RefinerUtils file)

// divide by center point (in 3D)
const std::vector<std::vector<std::array<double, 3>>> TetrahedronRefine4{
    T4{Va, Vb, Vc, Vabcd}, T4{Va, Vb, Vd, Vabcd}, T4{Va, Vc, Vd, Vabcd}, T4{Vb, Vc, Vd, Vabcd}};

// divide edges (in 3D)
const std::vector<std::vector<std::array<double, 3>>> TetrahedronRefine8{
    T4{Va, Vab, Vac, Vad},
    T4{Vb, Vab, Vbc, Vbd},
    T4{Vc, Vac, Vbc, Vcd},
    T4{Vd, Vad, Vbd, Vcd},
    T4{Vab, Vac, Vad, Vbd},
    T4{Vab, Vac, Vbc, Vbd},
    T4{Vac, Vad, Vbd, Vcd},
    T4{Vac, Vbc, Vbd, Vcd},
};

// divide edges (in 2D)
const std::vector<std::vector<std::array<double, 2>>> TriangleRefine4{
    T3{P2{0, 0}, P2{0.5, 0}, P2{0, 0.5}},
    T3{P2{0.5, 0}, P2{1, 0}, P2{0.5, 0.5}},
    T3{P2{0, 0.5}, P2{0.5, 0.5}, P2{0, 1}},
    T3{P2{0.5, 0.5}, P2{0, 0.5}, P2{0.5, 0}},
};

} // namespace seissol::io::instance::geometry
