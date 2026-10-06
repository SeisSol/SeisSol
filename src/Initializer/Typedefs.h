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
#include "Common/Constants.h"
#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/coefficients.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "Kernels/Data.h"
#include "Model/OperatorLayout.h"
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
  // star matrices, where the cell carries them assembled
  real starMatrices[3][zeroGuard(FactoredStar ? 0 : seissol::tensor::star::size(0))]{};

  // the rows of the Jacobian and the material coefficients, where it does not. A curved cell has
  // the rows at every point the operator is formed at.
  real referenceGradients[3][zeroGuard(FactoredStar ? (Curvilinear ? 3 * OperatorPointCount : 3)
                                                    : 0)]{};
  // one per coefficient, and where the material varies inside the cell one per
  // sample point of it. The sample index is the slower one, so that a
  // coefficient's samples lie together the way the kernel reads them.
  real materialCoefficients[zeroGuard(FactoredStar ? StarCoefficientCount : 0)]
                           [zeroGuard(FactoredStar ? MaterialSampleCount : 0)]{};

  // The scalars the source term is linear in, at the sample points, where a
  // cell forms it from the material it carries rather than from one matrix.
  real sourceCoefficients[zeroGuard(NodalSource ? SourceCoefficientCount : 0)]
                         [zeroGuard(NodalSource ? MaterialSampleCount : 0)]{};
  // What a sample point deviates from the source term the cell carries for
  // itself, for a solver that factorises that term into a solve of its own.
  real sourceDeviation[zeroGuard(NodalSourceDeviation ? SourceDeviationCount : 0)]
                      [zeroGuard(NodalSourceDeviation ? MaterialSampleCount : 0)]{};

  // flux solver for element local contribution. It is filled where the flux
  // reads the material at the nodes of a face as well, although the flux
  // kernels of such a build read the scalars below instead, on the host and on
  // the device alike.
  real nApNm1[4][seissol::tensor::AplusT::size()]{};

  // Where the material varies along a face, the flux operator does too, and a
  // face carries the scalars it is built from at the nodes of that face rather
  // than the matrix they fold into. The rotation into the face coordinates
  // those scalars are stated in is the same for both sides, so a face keeps one
  // of them; the inverse follows from it inside the kernel. A face that may be
  // curved keeps one per node.
  real fluxCoefficients[zeroGuard(NodalFlux ? Cell::NumFaces : 0)][zeroGuard(
      NodalFlux ? FluxCoefficientCount : 0)][zeroGuard(NodalFlux ? FluxFaceNodes : 0)]{};
  real faceRotation[zeroGuard(NodalFlux ? Cell::NumFaces : 0)][zeroGuard(
      NodalMaterial ? (Curvilinear ? seissol::tensor::TNodes::size() : seissol::tensor::T::size())
                    : 0)]{};

  // solver-specific data
  seissol::model::MaterialT::Solver::LocalData specific;
};

// data for the neighboring boundary integration
struct alignas(Alignment) NeighboringIntegrationData {
  // flux solver for the contribution of the neighboring elements
  real nAmNm1[4][seissol::tensor::AminusT::size()]{};

  // the counterpart of LocalIntegrationData::fluxCoefficients for the operator
  // the neighbour contributes. The matrix above stays: a boundary face takes
  // its neighbour state from a nodal boundary condition, already at the nodes
  // of the face and already rotated, and applies the matrix to it.
  real fluxCoefficients[zeroGuard(NodalFlux ? Cell::NumFaces : 0)][zeroGuard(
      NodalFlux ? FluxCoefficientCount : 0)][zeroGuard(NodalFlux ? FluxFaceNodes : 0)]{};

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
  std::int8_t side{};
  std::int8_t faceRelation{};
  real* godunov{nullptr};
  real* fluxSolver{nullptr};
};

struct BoundaryFaceInformation {
  // nodes is an array of 3d-points in global coordinates.
  real nodes[seissol::nodal::tensor::nodes2D::Shape[multisim::BasisFunctionDimension] * 3]{};
  real dataT[seissol::tensor::T::size()]{};
  real dataTinv[seissol::tensor::Tinv::size()]{};
  real dirichletOffset[seissol::tensor::dirichletOffset::size()]{};
  real dirichletMap[seissol::tensor::dirichletMap::size()]{};
  real fsgData[3]{};
};

struct CellBoundaryMapping {
  real* nodes{nullptr};
  real* dataT{nullptr};
  real* dataTinv{nullptr};
  real* dirichletOffset{nullptr};
  real* dirichletMap{nullptr};
  real* fsgData{nullptr};

  CellBoundaryMapping() = default;
  explicit CellBoundaryMapping(BoundaryFaceInformation& faceInfo)
      : nodes(faceInfo.nodes), dataT(faceInfo.dataT), dataTinv(faceInfo.dataTinv),
        dirichletOffset(faceInfo.dirichletOffset), dirichletMap(faceInfo.dirichletMap),
        fsgData(faceInfo.fsgData) {}
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
