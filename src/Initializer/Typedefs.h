// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_INITIALIZER_TYPEDEFS_H_
#define SEISSOL_SRC_INITIALIZER_TYPEDEFS_H_

#include "Alignment.h"
#include "BasicTypedefs.h"
#include "CellLocalInformation.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "Kernels/Data.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Dense>
#include <complex>
#include <cstddef>
#include <vector>

namespace seissol {

namespace kernels {
constexpr std::size_t NumSpaceQuadraturePoints = (ConvergenceOrder + 1) * (ConvergenceOrder + 1);
} // namespace kernels

/**
 * The generated constant matrices, as one table of pointers into the pool.
 *
 * There are two of these: one built on the image in this binary, one on a
 * copy of it in device memory. Which entries a table holds is decided by the
 * code generator, so adding a matrix no longer means touching this file.
 **/
using GlobalData = seissol::Pool;

struct CompoundGlobalData {
  GlobalData* onHost{nullptr};
  GlobalData* onDevice{nullptr};
};

// data for the cell local integration
struct alignas(Alignment) LocalIntegrationData {
  // star matrices
  real starMatrices[3][seissol::tensor::star::size(0)]{};

  // flux solver for element local contribution
  real nApNm1[4][seissol::tensor::AplusT::size()]{};

  // solver-specific data
  seissol::model::MaterialT::Solver::LocalData specific;
};

// data for the neighboring boundary integration
struct alignas(Alignment) NeighboringIntegrationData {
  // flux solver for the contribution of the neighboring elements
  real nAmNm1[4][seissol::tensor::AminusT::size()]{};

  // solver-specific data
  seissol::model::MaterialT::Solver::NeighborData specific;
};

// material constants per cell
struct CellMaterialData {
  seissol::model::Material* local{};
  seissol::model::Material* neighbor[4]{};
};

struct DRFaceInformation {
  std::size_t meshFace{};
  std::uint8_t plusSide{};
  std::uint8_t minusSide{};
  std::uint8_t faceRelation{};
  bool plusSideOnThisRank{};
};

struct DRGodunovData {
  real dataTinvT[seissol::tensor::TinvT::size()]{};
  real tractionPlusMatrix[seissol::tensor::tractionPlusMatrix::size()]{};
  real tractionMinusMatrix[seissol::tensor::tractionMinusMatrix::size()]{};
  // When integrating quantities over the fault
  // we need to integrate over each physical element.
  // The integration is effectively done in the reference element, and the scaling factor of
  // the transformation, the surface Jacobian (e.g. |n^e(\chi)| in eq. (35) of Uphoff et al. (2023))
  // is incorporated. This explains the factor 2 (doubledSurfaceArea)
  //
  // Uphoff, C., May, D. A., & Gabriel, A. A. (2023). A discontinuous Galerkin method for
  // sequences of earthquakes and aseismic slip on multiple faults using unstructured curvilinear
  // grids. Geophysical Journal International, 233(1), 586-626.
  double doubledSurfaceArea{};
};

struct DREnergyOutput {
  real slip[seissol::tensor::slipInterpolated::size()]{};
  real accumulatedSlip[seissol::dr::misc::NumPaddedPoints]{};
  real frictionalEnergy[seissol::dr::misc::NumPaddedPoints]{};
  real timeSinceSlipRateBelowThreshold[seissol::dr::misc::NumPaddedPoints]{};

  static std::vector<seissol::io::datatype::StructDatatype::MemberInfo> datatypeLayout() {
    return {
        seissol::io::datatype::StructDatatype::MemberInfo{
            "slip",
            offsetof(DREnergyOutput, slip),
            seissol::io::datatype::inferDatatype<decltype(slip)>()},
        seissol::io::datatype::StructDatatype::MemberInfo{
            "accumulatedSlip",
            offsetof(DREnergyOutput, accumulatedSlip),
            seissol::io::datatype::inferDatatype<decltype(accumulatedSlip)>()},
        seissol::io::datatype::StructDatatype::MemberInfo{
            "frictionalEnergy",
            offsetof(DREnergyOutput, frictionalEnergy),
            seissol::io::datatype::inferDatatype<decltype(frictionalEnergy)>()},
        seissol::io::datatype::StructDatatype::MemberInfo{
            "timeSinceSlipRateBelowThreshold",
            offsetof(DREnergyOutput, timeSinceSlipRateBelowThreshold),
            seissol::io::datatype::inferDatatype<decltype(timeSinceSlipRateBelowThreshold)>()},
    };
  }
};

struct CellDRMapping {
  unsigned side{};
  unsigned faceRelation{};
  real* godunov{nullptr};
  real* fluxSolver{nullptr};
};

struct BoundaryFaceInformation {
  // nodes is an array of 3d-points in global coordinates.
  real nodes[seissol::nodal::tensor::nodes2D::Shape[multisim::BasisFunctionDimension] * 3]{};
  real dataT[seissol::tensor::T::size()]{};
  real dataTinv[seissol::tensor::Tinv::size()]{};
  real easiBoundaryConstant[seissol::tensor::easiBoundaryConstant::size()]{};
  real easiBoundaryMap[seissol::tensor::easiBoundaryMap::size()]{};
};

struct CellBoundaryMapping {
  real* nodes{nullptr};
  real* dataT{nullptr};
  real* dataTinv{nullptr};
  real* easiBoundaryConstant{nullptr};
  real* easiBoundaryMap{nullptr};

  CellBoundaryMapping() = default;
  explicit CellBoundaryMapping(BoundaryFaceInformation& faceInfo)
      : nodes(faceInfo.nodes), dataT(faceInfo.dataT), dataTinv(faceInfo.dataTinv),
        easiBoundaryConstant(faceInfo.easiBoundaryConstant),
        easiBoundaryMap(faceInfo.easiBoundaryMap) {}
};

struct GravitationSetup {
  double acceleration = 9.81; // m/s
};

struct TravellingWaveParameters {
  Eigen::Vector3d origin;
  Eigen::Vector3d kVec;
  std::vector<int> varField;
  std::vector<std::complex<double>> ampField;
};

struct AcousticTravellingWaveParametersITM {
  double k{};
  double itmStartingTime{};
  double itmDuration{};
  double itmVelocityScalingFactor{};
};

struct PressureInjectionParameters {
  std::array<double, 3> origin{};
  double magnitude{};
  double width{};
};

} // namespace seissol

#endif // SEISSOL_SRC_INITIALIZER_TYPEDEFS_H_
