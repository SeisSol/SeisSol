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
#include "Common/Real.h"
#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "Kernels/Data.h"
#include "Kernels/SolverSelector.h"
#include "Solver/MultipleSimulations.h"

#include <Eigen/Dense>
#include <complex>
#include <cstddef>
#include <vector>

namespace seissol {

namespace kernels {
template <typename Cfg>
constexpr std::size_t NumSpaceQuadraturePoints =
    (Cfg::ConvergenceOrder + 1) * (Cfg::ConvergenceOrder + 1);
} // namespace kernels

/**
 * The generated constant matrices of the configuration `Cfg`, as one table of pointers into the
 * pool of its generator.
 *
 * There are two of these per configuration: one built on the image in this binary, one on a
 * copy of it in device memory. Which entries a table holds is decided by the code generator.
 **/
template <typename Cfg>
using GlobalData = seissol::Pool<Cfg>;

/// The global data of the configuration `Cfg` on the host and, in a GPU build, on the device.
template <typename Cfg>
struct CompoundGlobalData {
  GlobalData<Cfg>* onHost{nullptr};
  GlobalData<Cfg>* onDevice{nullptr};
};

// data for the cell local integration
template <typename Cfg>
struct alignas(Alignment) LocalIntegrationData {
  // star matrices
  Real<Cfg> starMatrices[3][seissol::tensor::star<Cfg>::size(0)]{};

  // flux solver for element local contribution
  Real<Cfg> nApNm1[4][seissol::tensor::AplusT<Cfg>::size()]{};

  // solver-specific data
  typename seissol::kernels::SolverOf<Cfg>::LocalData specific;
};

// data for the neighboring boundary integration
template <typename Cfg>
struct alignas(Alignment) NeighboringIntegrationData {
  // flux solver for the contribution of the neighboring elements
  Real<Cfg> nAmNm1[4][seissol::tensor::AminusT<Cfg>::size()]{};

  // solver-specific data
  typename seissol::kernels::SolverOf<Cfg>::NeighborData specific;
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

template <typename Cfg>
struct DRGodunovData {
  Real<Cfg> dataTinvT[seissol::tensor::TinvT<Cfg>::size()]{};
  Real<Cfg> tractionPlusMatrix[seissol::tensor::tractionPlusMatrix<Cfg>::size()]{};
  Real<Cfg> tractionMinusMatrix[seissol::tensor::tractionMinusMatrix<Cfg>::size()]{};
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

template <typename Cfg>
struct DREnergyOutput {
  Real<Cfg> slip[seissol::tensor::slipInterpolated<Cfg>::size()]{};
  Real<Cfg> accumulatedSlip[seissol::dr::misc::NumPaddedPoints<Cfg>]{};
  Real<Cfg> frictionalEnergy[seissol::dr::misc::NumPaddedPoints<Cfg>]{};
  Real<Cfg> timeSinceSlipRateBelowThreshold[seissol::dr::misc::NumPaddedPoints<Cfg>]{};

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

template <typename Cfg>
struct CellDRMapping {
  std::int8_t side{};
  std::int8_t faceRelation{};
  Real<Cfg>* godunov{nullptr};
  Real<Cfg>* fluxSolver{nullptr};
};

template <typename Cfg>
struct BoundaryFaceInformation {
  // nodes is an array of 3d-points in global coordinates.
  Real<Cfg> nodes[seissol::nodal::tensor::nodes2D<Cfg>::Shape[multisim::BasisDim<Cfg>] * 3]{};
  Real<Cfg> dataT[seissol::tensor::T<Cfg>::size()]{};
  Real<Cfg> dataTinv[seissol::tensor::Tinv<Cfg>::size()]{};
  Real<Cfg> dirichletOffset[seissol::tensor::dirichletOffset<Cfg>::size()]{};
  Real<Cfg> dirichletMap[seissol::tensor::dirichletMap<Cfg>::size()]{};
  Real<Cfg> fsgData[3]{};
};

template <typename Cfg>
struct CellBoundaryMapping {
  Real<Cfg>* nodes{nullptr};
  Real<Cfg>* dataT{nullptr};
  Real<Cfg>* dataTinv{nullptr};
  Real<Cfg>* dirichletOffset{nullptr};
  Real<Cfg>* dirichletMap{nullptr};
  Real<Cfg>* fsgData{nullptr};

  CellBoundaryMapping() = default;
  explicit CellBoundaryMapping(BoundaryFaceInformation<Cfg>& faceInfo)
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
