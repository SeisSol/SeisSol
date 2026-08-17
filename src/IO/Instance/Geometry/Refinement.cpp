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
#include <cstddef>
#include <vector>

namespace {
using T3 = std::vector<std::array<double, 2>>;
using T4 = std::vector<std::array<double, 3>>;
using P2 = std::array<double, 2>;
using P3 = std::array<double, 3>;

// The subdivision tables are written in the coordinates of the reference simplex that
// AffineMap::fromVertices reconstructs from, i.e. Va = f(0) and Vb, Vc, Vd = f(e_1), f(e_2),
// f(e_3). Via transformations::tetrahedronReferenceToGlobal that identifies Va, Vb, Vc, Vd with
// element.vertices[0], [1], [2], [3] -- the same identification the legacy
// refinement::MeshRefiner used. Do not confuse this with the vertex order of the legacy
// refinement::Tetrahedron::unitTetrahedron(), which is a different (rotated) labelling.
constexpr auto Va = P3{0, 0, 0};
constexpr auto Vb = P3{1, 0, 0};
constexpr auto Vc = P3{0, 1, 0};
constexpr auto Vd = P3{0, 0, 1};

template <std::size_t Dim, typename... Points>
constexpr std::array<double, Dim> average(const Points&... points) {
  std::array<double, Dim> output{};
  for (std::size_t i = 0; i < Dim; ++i) {
    output[i] = (points[i] + ...) / static_cast<double>(sizeof...(Points));
  }
  return output;
}

constexpr auto Vab = average<3>(Va, Vb);
constexpr auto Vac = average<3>(Va, Vc);
constexpr auto Vad = average<3>(Va, Vd);
constexpr auto Vbc = average<3>(Vb, Vc);
constexpr auto Vbd = average<3>(Vb, Vd);
constexpr auto Vcd = average<3>(Vc, Vd);

constexpr auto Vabcd = average<3>(Va, Vb, Vc, Vd);

} // namespace

namespace seissol::io::instance::geometry {

// The subcells are enumerated as in the legacy refinement::DivideTetrahedronBy4 and
// ...By8, so that the cells of a refined output mesh keep their previous identity. Their
// vertex order does not: the legacy tables enumerate the vertex subsets lexicographically,
// which alternates the parity and left every second subcell with a negative Jacobian. Every
// entry below is oriented positively instead, so all four/eight subcells cover the same
// region as before while remaining valid VTK/XDMF tetrahedra.

// divide by center point (in 3D)
const std::vector<std::vector<std::array<double, 3>>> TetrahedronRefine4{
    T4{Va, Vb, Vc, Vabcd}, T4{Va, Vb, Vabcd, Vd}, T4{Va, Vc, Vd, Vabcd}, T4{Vb, Vc, Vabcd, Vd}};

// divide edges (in 3D)
const std::vector<std::vector<std::array<double, 3>>> TetrahedronRefine8{
    T4{Va, Vab, Vac, Vad},
    T4{Vb, Vab, Vbd, Vbc},
    T4{Vc, Vac, Vbc, Vcd},
    T4{Vd, Vad, Vcd, Vbd},
    T4{Vab, Vac, Vad, Vbd},
    T4{Vab, Vac, Vbd, Vbc},
    T4{Vac, Vad, Vbd, Vcd},
    T4{Vac, Vbc, Vcd, Vbd},
};

// divide edges (in 2D); enumeration and vertex order as in the legacy
// refinement::TriangleRefiner, which is already oriented positively throughout
const std::vector<std::vector<std::array<double, 2>>> TriangleRefine4{
    T3{P2{0, 0}, P2{0.5, 0}, P2{0, 0.5}},
    T3{P2{0.5, 0}, P2{1, 0}, P2{0.5, 0.5}},
    T3{P2{0, 0.5}, P2{0.5, 0.5}, P2{0, 1}},
    T3{P2{0.5, 0.5}, P2{0, 0.5}, P2{0.5, 0}},
};

} // namespace seissol::io::instance::geometry
