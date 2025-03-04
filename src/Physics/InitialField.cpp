// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <Equations/Datastructures.h>
#include <Initializer/Parameters/InitializationParameters.h>
#include <Initializer/Typedefs.h>
#include <Model/Common.h>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <init.h>
#include <limits>
#include <math.h>
#include <tensor.h>

#include "Kernels/Precision.h"
#include "Numerical/Eigenvalues.h"
#include "Physics/InitialField.h"
#include <stdexcept>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>
#include <yateto.h>

// FIXME: the following line is absolutely necessary for the plain-wave operator to work correctly
// (template specializations for the equations).
#include "Equations/Setup.h" // IWYU pragma: keep

seissol::physics::Planarwave::Planarwave(const CellMaterialData& materialData,
                                         double phase,
                                         Eigen::Vector3d kVec,
                                         std::vector<int> varField,
                                         std::vector<std::complex<double>> ampField)
    : m_varField(std::move(varField)), m_ampField(std::move(ampField)), m_phase(phase),
      m_kVec(std::move(kVec)) {
  init(materialData);
}

seissol::physics::Planarwave::Planarwave(const CellMaterialData& materialData,
                                         double phase,
                                         Eigen::Vector3d kVec)
    : m_phase(phase), m_kVec(std::move(kVec)) {

#ifndef USE_POROELASTIC
  bool isAcoustic = false;
#ifndef USE_ANISOTROPIC
  isAcoustic = materialData.local.getMuBar() <= 1e-15;
#endif
  if (isAcoustic) {
    // Acoustic materials has the following wave modes:
    // -P, N, N, N, N, N, N, N, P
    // Here we impose the P mode
    m_varField = {8};
    m_ampField = {1.0};
  } else {
    // Elastic materials have the following wave modes:
    // -P, -S2, -S1, N, N, N, S1, S2, P
    // Here we impose the -S2 and P mode
    m_varField = {1, 8};
    m_ampField = {1.0, 1.0};
  }
#else
  // Poroelastic materials have the following wave modes:
  //-P, -S2, -S1, -Ps, N, N, N, N, N, Ps, S1, S2, P
  // Here we impose -S1, -Ps and P
  m_varField = {2, 3, 12};
  m_ampField = {1.0, 1.0, 1.0};
#endif
  init(materialData);
}

void seissol::physics::Planarwave::init(const CellMaterialData& materialData) {
  assert(m_varField.size() == m_ampField.size());

  std::array<std::complex<double>,
             seissol::model::MaterialT::NumQuantities * seissol::model::MaterialT::NumQuantities>
      planeWaveOperator{};
  seissol::model::getPlaneWaveOperator(materialData.local, m_kVec.data(), planeWaveOperator.data());
  seissol::eigenvalues::Eigenpair<std::complex<double>, seissol::model::MaterialT::NumQuantities>
      eigendecomposition;
  computeEigenvalues(planeWaveOperator, eigendecomposition);
  m_lambdaA = eigendecomposition.values;
  m_eigenvectors = eigendecomposition.vectors;
}

void seissol::physics::Planarwave::evaluate(
    double time,
    const std::vector<std::array<double, 3>>& points,
    const CellMaterialData& materialData,
    yateto::DenseTensorView<2, real, unsigned>& dofsQP) const {
  dofsQP.setZero();

  auto r = yateto::DenseTensorView<2, std::complex<double>>(
      const_cast<std::complex<double>*>(m_eigenvectors.data()),
      {seissol::model::MaterialT::NumQuantities, seissol::model::MaterialT::NumQuantities});
  for (unsigned v = 0; v < m_varField.size(); ++v) {
    const auto omega = m_lambdaA[m_varField[v]];
    for (unsigned j = 0; j < dofsQP.shape(1); ++j) {
      for (size_t i = 0; i < points.size(); ++i) {
        dofsQP(i, j) +=
            (r(j, m_varField[v]) * m_ampField[v] *
             std::exp(std::complex<double>(0.0, 1.0) *
                      (omega * time - m_kVec[0] * points[i][0] - m_kVec[1] * points[i][1] -
                       m_kVec[2] * points[i][2] + std::complex<double>(m_phase, 0))))
                .real();
      }
    }
  }
}

seissol::physics::SuperimposedPlanarwave::SuperimposedPlanarwave(
    const CellMaterialData& materialData, real phase)
    : m_kVec({{{M_PI, 0.0, 0.0}, {0.0, M_PI, 0.0}, {0.0, 0.0, M_PI}}}), m_phase(phase),
      m_pw({Planarwave(materialData, phase, m_kVec.at(0)),
            Planarwave(materialData, phase, m_kVec.at(1)),
            Planarwave(materialData, phase, m_kVec.at(2))}) {}

void seissol::physics::SuperimposedPlanarwave::evaluate(
    double time,
    const std::vector<std::array<double, 3>>& points,
    const CellMaterialData& materialData,
    yateto::DenseTensorView<2, real, unsigned>& dofsQP) const {
  dofsQP.setZero();

  real dofsPwData[tensor::dofsQP::size()];
  yateto::DenseTensorView<2, real, unsigned> dofsPW = init::dofsQP::view::create(dofsPwData);

  for (int pw = 0; pw < 3; pw++) {
    // evaluate each planarwave
    m_pw.at(pw).evaluate(time, points, materialData, dofsPW);
    // and add results together
    for (unsigned j = 0; j < dofsQP.shape(1); ++j) {
      for (size_t i = 0; i < points.size(); ++i) {
        dofsQP(i, j) += dofsPW(i, j);
      }
    }
  }
}

seissol::physics::TravellingWave::TravellingWave(
    const CellMaterialData& materialData, const TravellingWaveParameters& travellingWaveParameters)
    // Set phase to 0.5*M_PI, so we have a zero at the origin
    // The wave travels in direction of kVec
    // 2*pi / magnitude(kVec) is the wave length of the wave
    : Planarwave(materialData,
                 0.5 * M_PI,
                 travellingWaveParameters.kVec,
                 travellingWaveParameters.varField,
                 travellingWaveParameters.ampField),
      // origin is a point on the wavefront at time zero
      m_origin(travellingWaveParameters.origin) {
  logInfo() << "Impose a travelling wave as initial condition";
  logInfo() << "Origin = (" << m_origin[0] << ", " << m_origin[1] << ", " << m_origin[2] << ")";
  logInfo() << "kVec = (" << m_kVec[0] << ", " << m_kVec[1] << ", " << m_kVec[2] << ")";
  logInfo() << "Combine following wave modes";
  for (size_t i = 0; i < m_ampField.size(); i++) {
    logInfo() << "(" << m_varField[i] << ": " << m_ampField[i] << ")";
  }
}

seissol::physics::AcousticTravellingWaveITM::AcousticTravellingWaveITM(
    const CellMaterialData& materialData,
    const AcousticTravellingWaveParametersITM& acousticTravellingWaveParametersItm)
    : rho0(materialData.local.rho),
      c0(sqrt(materialData.local.getLambdaBar() / materialData.local.rho)),
      k(acousticTravellingWaveParametersItm.k),
      tITMMinus(acousticTravellingWaveParametersItm.itmStartingTime), tITMPlus(tITMMinus + tau),
      tau(acousticTravellingWaveParametersItm.itmDuration),
      n(acousticTravellingWaveParametersItm.itmVelocityScalingFactor) {
#ifdef USE_ANISOTROPIC
  logError() << "This has not been yet implemented for anisotropic material";
#else
  logInfo() << "Starting Test for Acoustic Travelling Wave with ITM";

  logInfo() << "rho0 = " << rho0;
  logInfo() << "c0 = " << c0;

  logInfo() << "k = " << k;

  logInfo() << "Setting up the Initial Conditions";
  init(materialData);
#endif
}

void seissol::physics::AcousticTravellingWaveITM::init(const CellMaterialData& materialData) {
#ifdef USE_ANISOTROPIC
  logError() << "This has not been yet implemented for anisotropic material";
#endif
}

void seissol::physics::AcousticTravellingWaveITM::evaluate(
    double time,
    const std::vector<std::array<double, 3>>& points,
    const CellMaterialData& materialData,
    yateto::DenseTensorView<2, real, unsigned>& dofsQP) const {
#ifdef USE_ANISOTROPIC
  logError() << "This has not been yet implemented for anisotropic material";
#else
  dofsQP.setZero();
  double pressure = 0.0;
  for (size_t i = 0; i < points.size(); ++i) {
    const auto& coordinates = points[i];
    const auto x = coordinates[0];
    const auto t = time;
    if (t <= tITMMinus) {
      pressure = c0 * rho0 * std::cos(k * x - c0 * k * t);
      dofsQP(i, 0) = -pressure;                    // sigma_xx
      dofsQP(i, 1) = -pressure;                    // sigma_yy
      dofsQP(i, 2) = -pressure;                    // sigma_zz
      dofsQP(i, 3) = 0.0;                          // sigma_xy
      dofsQP(i, 4) = 0.0;                          // sigma_yz
      dofsQP(i, 5) = 0.0;                          // sigma_xz
      dofsQP(i, 6) = std::cos(k * x - c0 * k * t); // u
      dofsQP(i, 7) = 0.0;                          // v
      dofsQP(i, 8) = 0.0;                          // w
    } else if (t <= tITMPlus) {
      pressure = -0.5 * (n - 1) * c0 * rho0 *
                     std::cos(k * x + c0 * k * n * t - (c0 * k * n + c0 * k) * tITMMinus) +
                 0.5 * (n + 1) * c0 * rho0 *
                     std::cos(k * x - c0 * k * n * t + (c0 * k * n - c0 * k) * tITMMinus);
      dofsQP(i, 0) = -pressure; // sigma_xx
      dofsQP(i, 1) = -pressure; // sigma_yy
      dofsQP(i, 2) = -pressure; // sigma_zz
      dofsQP(i, 3) = 0.0;       // sigma_xy
      dofsQP(i, 4) = 0.0;       // sigma_yz
      dofsQP(i, 5) = 0.0;       // sigma_xz
      dofsQP(i, 6) =
          0.5 * (n - 1) * std::cos(k * x + c0 * k * n * t - (c0 * k * n + c0 * k) * tITMMinus) +
          0.5 * (n + 1) * std::cos(k * x - c0 * k * n * t + (c0 * k * n - c0 * k) * tITMMinus); // u
      dofsQP(i, 7) = 0.0;                                                                       // v
      dofsQP(i, 8) = 0.0;                                                                       // w
    } else {
      pressure =
          -0.25 * (1 / n) * c0 * rho0 *
          ((-n * n + 1) * std::cos(k * x + c0 * k * t - 2.0 * c0 * k * tITMMinus -
                                   (c0 * k * n + c0 * k) * tau) +
           (n * n - 1) * std::cos(k * x + c0 * k * t - 2.0 * c0 * k * tITMMinus +
                                  (c0 * k * n - c0 * k) * tau) +
           (n * n - 2 * n + 1) * std::cos(k * x - c0 * k * t + (c0 * k * n + c0 * k) * tau) +
           (-n * n - 2 * n - 1) * std::cos(k * x - c0 * k * t - (c0 * k * n - c0 * k) * tau));
      dofsQP(i, 0) = -pressure; // sigma_xx
      dofsQP(i, 1) = -pressure; // sigma_yy
      dofsQP(i, 2) = -pressure; // sigma_zz
      dofsQP(i, 3) = 0.0;       // sigma_xy
      dofsQP(i, 4) = 0.0;       // sigma_yz
      dofsQP(i, 5) = 0.0;       // sigma_xz
      dofsQP(i, 6) =
          (-0.25 / n) *
          ((n * n - 1) *
               std::cos(k * x + c0 * k * t - 2 * c0 * k * tITMMinus - (c0 * k * n + c0 * k) * tau) +
           (-n * n + 1) *
               std::cos(k * x + c0 * k * t - 2 * c0 * k * tITMMinus + (c0 * k * n - c0 * k) * tau) +
           (n * n - 2 * n + 1) * std::cos(k * x - c0 * k * t + (c0 * k * n + c0 * k) * tau) +
           (-n * n - 2 * n - 1) * std::cos(k * x - c0 * k * t - (c0 * k * n - c0 * k) * tau)); // u
      dofsQP(i, 7) = 0.0;                                                                      // v
      dofsQP(i, 8) = 0.0;                                                                      // w
    }
  }
#endif
}

void seissol::physics::TravellingWave::evaluate(
    double time,
    const std::vector<std::array<double, 3>>& points,
    const CellMaterialData& materialData,
    yateto::DenseTensorView<2, real, unsigned>& dofsQp) const {
  dofsQp.setZero();

  auto r = yateto::DenseTensorView<2, std::complex<double>>(
      const_cast<std::complex<double>*>(m_eigenvectors.data()),
      {seissol::model::MaterialT::NumQuantities, seissol::model::MaterialT::NumQuantities});
  for (unsigned v = 0; v < m_varField.size(); ++v) {
    const auto omega = m_lambdaA[m_varField[v]];
    for (unsigned j = 0; j < dofsQp.shape(1); ++j) {
      for (size_t i = 0; i < points.size(); ++i) {
        auto arg = std::complex<double>(0.0, 1.0) *
                   (omega * time - m_kVec[0] * (points[i][0] - m_origin[0]) -
                    m_kVec[1] * (points[i][1] - m_origin[1]) -
                    m_kVec[2] * (points[i][2] - m_origin[2]) + m_phase);
        if (arg.imag() > -0.5 * M_PI && arg.imag() < 1.5 * M_PI) {
          dofsQp(i, j) += (r(j, m_varField[v]) * m_ampField[v] * std::exp(arg)).real();
        }
      }
    }
  }
}

seissol::physics::PressureInjection::PressureInjection(
    const seissol::initializer::parameters::InitializationParameters& initializationParameters)
    : m_parameters(initializationParameters) {
  const auto o1 = m_parameters.origin[0];
  const auto o2 = m_parameters.origin[1];
  const auto o3 = m_parameters.origin[2];
  const auto magnitude = m_parameters.magnitude;
  const auto width = m_parameters.width;
  logInfo() << "Prepare gaussian pressure perturbation with center at (" << o1 << ", " << o2 << ", "
            << o3 << "), magnitude = " << magnitude << ", width = " << width << ".";
}

void seissol::physics::PressureInjection::evaluate(
    double time,
    const std::vector<std::array<double, 3>>& points,
    const CellMaterialData& materialData,
    yateto::DenseTensorView<2, real, unsigned>& dofsQp) const {
  const auto o1 = m_parameters.origin[0];
  const auto o2 = m_parameters.origin[1];
  const auto o3 = m_parameters.origin[2];
  const auto magnitude = m_parameters.magnitude;
  const auto width = m_parameters.width;

  for (size_t i = 0; i < points.size(); ++i) {
    const auto& x = points[i];
    const auto x1 = x[0];
    const auto x2 = x[1];
    const auto x3 = x[2];
    const auto rSquared = std::pow(x1 - o1, 2) + std::pow(x2 - o2, 2) + std::pow(x3 - o3, 2);
    dofsQp(i, 0) = 0.0;                                     // sigma_xx
    dofsQp(i, 1) = 0.0;                                     // sigma_yy
    dofsQp(i, 2) = 0.0;                                     // sigma_yy
    dofsQp(i, 3) = 0.0;                                     // sigma_xy
    dofsQp(i, 4) = 0.0;                                     // sigma_yz
    dofsQp(i, 5) = 0.0;                                     // sigma_xz
    dofsQp(i, 6) = 0.0;                                     // u
    dofsQp(i, 7) = 0.0;                                     // v
    dofsQp(i, 8) = 0.0;                                     // w
    dofsQp(i, 9) = magnitude * std::exp(-width * rSquared); // p
    dofsQp(i, 10) = 0.0;                                    // u_f
    dofsQp(i, 11) = 0.0;                                    // v_f
    dofsQp(i, 12) = 0.0;                                    // w_f
  }
}

void seissol::physics::ScholteWave::evaluate(
    double time,
    const std::vector<std::array<double, 3>>& points,
    const CellMaterialData& materialData,
    yateto::DenseTensorView<2, real, unsigned>& dofsQp) const {
#ifndef USE_ANISOTROPIC
  const real omega = 2.0 * std::acos(-1);

  for (size_t i = 0; i < points.size(); ++i) {
    const auto& x = points[i];
    const bool isAcousticPart =
        std::abs(materialData.local.getMuBar()) < std::numeric_limits<real>::epsilon();
    const auto x1 = x[0];
    const auto x3 = x[2];
    const auto t = time;
    if (isAcousticPart) {
      dofsQp(i, 0) = 0.35944997730200889 * std::pow(omega, 2) *
                     std::exp(-0.98901344820674908 * omega * x3) *
                     std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_xx
      dofsQp(i, 1) = 0.35944997730200889 * std::pow(omega, 2) *
                     std::exp(-0.98901344820674908 * omega * x3) *
                     std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_yy
      dofsQp(i, 2) = 0.35944997730200889 * std::pow(omega, 2) *
                     std::exp(-0.98901344820674908 * omega * x3) *
                     std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_zz
      dofsQp(i, 3) = 0;                                                    // sigma_xy
      dofsQp(i, 4) = 0;                                                    // sigma_yz
      dofsQp(i, 5) = 0;                                                    // sigma_xz
      dofsQp(i, 6) = -0.50555429848461109 * std::pow(omega, 2) *
                     std::exp(-0.98901344820674908 * omega * x3) *
                     std::sin(omega * t - 1.406466352506808 * omega * x1); // u
      dofsQp(i, 7) = 0;                                                    // v
      dofsQp(i, 8) = 0.35550086150929727 * std::pow(omega, 2) *
                     std::exp(-0.98901344820674908 * omega * x3) *
                     std::cos(omega * t - 1.406466352506808 * omega * x1); // w
    } else {
      dofsQp(i, 0) =
          -2.7820282741590652 * std::pow(omega, 2) * std::exp(0.98901344820674908 * omega * x3) *
              std::sin(omega * t - 1.406466352506808 * omega * x1) +
          3.5151973269883681 * std::pow(omega, 2) * std::exp(1.2825031256883821 * omega * x3) *
              std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_xx
      dofsQp(i, 1) = -6.6613381477509402e-16 * std::pow(omega, 2) *
                         std::exp(0.98901344820674908 * omega * x3) *
                         std::sin(omega * t - 1.406466352506808 * omega * x1) +
                     0.27315475753283058 * std::pow(omega, 2) *
                         std::exp(1.2825031256883821 * omega * x3) *
                         std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_yy
      dofsQp(i, 2) =
          2.7820282741590621 * std::pow(omega, 2) * std::exp(0.98901344820674908 * omega * x3) *
              std::sin(omega * t - 1.406466352506808 * omega * x1) -
          2.4225782968570462 * std::pow(omega, 2) * std::exp(1.2825031256883821 * omega * x3) *
              std::sin(omega * t - 1.406466352506808 * omega * x1); // sigma_zz
      dofsQp(i, 3) = 0;                                             // sigma_xy
      dofsQp(i, 4) = 0;                                             // sigma_yz
      dofsQp(i, 5) =
          -2.956295201467618 * std::pow(omega, 2) * std::exp(0.98901344820674908 * omega * x3) *
              std::cos(omega * t - 1.406466352506808 * omega * x1) +
          2.9562952014676029 * std::pow(omega, 2) * std::exp(1.2825031256883821 * omega * x3) *
              std::cos(omega * t - 1.406466352506808 * omega * x1); // sigma_xz
      dofsQp(i, 6) =
          0.98901344820675241 * std::pow(omega, 2) * std::exp(0.98901344820674908 * omega * x3) *
              std::sin(omega * t - 1.406466352506808 * omega * x1) -
          1.1525489264912381 * std::pow(omega, 2) * std::exp(1.2825031256883821 * omega * x3) *
              std::sin(omega * t - 1.406466352506808 * omega * x1); // u
      dofsQp(i, 7) = 0;                                             // v
      dofsQp(i, 8) =
          1.406466352506812 * std::pow(omega, 2) * std::exp(0.98901344820674908 * omega * x3) *
              std::cos(omega * t - 1.406466352506808 * omega * x1) -
          1.050965490997515 * std::pow(omega, 2) * std::exp(1.2825031256883821 * omega * x3) *
              std::cos(omega * t - 1.406466352506808 * omega * x1); // w
    }
  }
#else
  dofsQp.setZero();
#endif
}

void seissol::physics::SnellsLaw::evaluate(
    double time,
    const std::vector<std::array<double, 3>>& points,
    const CellMaterialData& materialData,
    yateto::DenseTensorView<2, real, unsigned>& dofsQp) const {
#ifndef USE_ANISOTROPIC
  const double pi = std::acos(-1);
  const double omega = 2.0 * pi;

  for (size_t i = 0; i < points.size(); ++i) {
    const auto& x = points[i];
    const bool isAcousticPart =
        std::abs(materialData.local.getMuBar()) < std::numeric_limits<real>::epsilon();

    const auto x1 = x[0];
    const auto x3 = x[2];
    const auto t = time;
    if (isAcousticPart) {
      dofsQp(i, 0) =
          1.0 * omega *
              std::sin(omega * t - omega * (0.19866933079506119 * x1 + 0.98006657784124163 * x3)) +
          0.48055591432167399 * omega *
              std::sin(omega * t -
                       omega * (0.19866933079506149 * x1 - 0.98006657784124152 * x3)); // sigma_xx
      dofsQp(i, 1) =
          1.0 * omega *
              std::sin(omega * t - omega * (0.19866933079506119 * x1 + 0.98006657784124163 * x3)) +
          0.48055591432167399 * omega *
              std::sin(omega * t -
                       omega * (0.19866933079506149 * x1 - 0.98006657784124152 * x3)); // sigma_yy
      dofsQp(i, 2) =
          1.0 * omega *
              std::sin(omega * t - omega * (0.19866933079506119 * x1 + 0.98006657784124163 * x3)) +
          0.48055591432167399 * omega *
              std::sin(omega * t -
                       omega * (0.19866933079506149 * x1 - 0.98006657784124152 * x3)); // sigma_zz
      dofsQp(i, 3) = 0;                                                                // sigma_xy
      dofsQp(i, 4) = 0;                                                                // sigma_yz
      dofsQp(i, 5) = 0;                                                                // sigma_xz
      dofsQp(i, 6) =
          -0.19866933079506119 * omega *
              std::sin(omega * t - omega * (0.19866933079506119 * x1 + 0.98006657784124163 * x3)) -
          0.095471721907895893 * omega *
              std::sin(omega * t -
                       omega * (0.19866933079506149 * x1 - 0.98006657784124152 * x3)); // u
      dofsQp(i, 7) = 0;                                                                // v
      dofsQp(i, 8) =
          -0.98006657784124163 * omega *
              std::sin(omega * t - omega * (0.19866933079506119 * x1 + 0.98006657784124163 * x3)) +
          0.47097679041061191 * omega *
              std::sin(omega * t -
                       omega * (0.19866933079506149 * x1 - 0.98006657784124152 * x3)); // w
    } else {
      dofsQp(i, 0) =
          -0.59005639909185559 * omega *
              std::sin(omega * t -
                       1.0 / 2.0 * omega * (0.39733866159012299 * x1 + 0.91767204817721759 * x3)) +
          0.55554011463785213 * omega *
              std::sin(omega * t -
                       1.0 / 3.0 * omega *
                           (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // sigma_xx
      dofsQp(i, 1) = 0.14460396298676709 * omega *
                     std::sin(omega * t -
                              1.0 / 3.0 * omega *
                                  (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // sigma_yy
      dofsQp(i, 2) =
          0.59005639909185559 * omega *
              std::sin(omega * t -
                       1.0 / 2.0 * omega * (0.39733866159012299 * x1 + 0.91767204817721759 * x3)) +
          0.89049951522981918 * omega *
              std::sin(omega * t -
                       1.0 / 3.0 * omega *
                           (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // sigma_zz
      dofsQp(i, 3) = 0;                                                           // sigma_xy
      dofsQp(i, 4) = 0;                                                           // sigma_yz
      dofsQp(i, 5) =
          -0.55363837274201066 * omega *
              std::sin(omega * t -
                       1.0 / 2.0 * omega * (0.39733866159012299 * x1 + 0.91767204817721759 * x3)) +
          0.55363837274201 * omega *
              std::sin(omega * t -
                       1.0 / 3.0 * omega *
                           (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // sigma_xz
      dofsQp(i, 6) =
          0.37125533967075403 * omega *
              std::sin(omega * t -
                       1.0 / 2.0 * omega * (0.39733866159012299 * x1 + 0.91767204817721759 * x3)) -
          0.2585553530120539 * omega *
              std::sin(omega * t - 1.0 / 3.0 * omega *
                                       (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // u
      dofsQp(i, 7) = 0;                                                                       // v
      dofsQp(i, 8) =
          -0.16074816713222639 * omega *
              std::sin(omega * t -
                       1.0 / 2.0 * omega * (0.39733866159012299 * x1 + 0.91767204817721759 * x3)) -
          0.34834162029840349 * omega *
              std::sin(omega * t - 1.0 / 3.0 * omega *
                                       (0.59600799238518454 * x1 + 0.8029785009656123 * x3)); // w
    }
  }
#else
  dofsQp.setZero();
#endif
}

seissol::physics::Ocean::Ocean(int mode, double gravitationalAcceleration)
    : mode(mode), gravitationalAcceleration(gravitationalAcceleration) {
  if (mode < 0 || mode > 3) {
    throw std::runtime_error("Wave mode " + std::to_string(mode) + " is not supported.");
  }
}
void seissol::physics::Ocean::evaluate(double time,
                                       const std::vector<std::array<double, 3>>& points,
                                       const CellMaterialData& materialData,
                                       yateto::DenseTensorView<2, real, unsigned>& dofsQp) const {
#ifndef USE_ANISOTROPIC
  for (size_t i = 0; i < points.size(); ++i) {
    const auto x = points[i][0];
    const auto y = points[i][1];
    const auto z = points[i][2];
    const auto t = time;

    const auto g = gravitationalAcceleration;
    if (std::abs(g - 9.81e-3) > 10e-15) {
      logError() << "Ocean scenario only supports g=9.81e-3 currently!";
    }
    if (materialData.local.getMuBar() > 10e-15) {
      logError() << "Ocean scenario only works for acoustic material (mu = 0.0)!";
    }
    const double pi = std::acos(-1);
    const double rho = materialData.local.rho;

    const double lx = 10.0;    // km
    const double ly = 10.0;    // km
    const double kX = pi / lx; // 1/km
    const double kY = pi / ly; // 1/km

    constexpr auto KStars =
        std::array<double, 3>{0.4433813748841239, 1.5733628061766445, 4.713305873881573};

    // Note: Could be computed on the fly but it's better to pre-compute them with higher precision!
    constexpr auto Omegas =
        std::array<double, 3>{0.0425599572628432, 2.4523337594491745, 7.1012991617572165};

    const auto kStar = KStars[mode];
    const auto omega = Omegas[mode];

    const auto b = g * kStar / (omega * omega);
    constexpr auto ScalingFactor = 1;

    // Shear stresses are zero for elastic
    dofsQp(i, 3) = 0.0;
    dofsQp(i, 4) = 0.0;
    dofsQp(i, 5) = 0.0;

    if (mode == 0) {
      // Gravity mode
      const auto pressure = -std::sin(kX * x) * std::sin(kY * y) * std::sin(omega * t) *
                            (std::sinh(kStar * z) + b * std::cosh(kStar * z));
      dofsQp(i, 0) = ScalingFactor * pressure;
      dofsQp(i, 1) = ScalingFactor * pressure;
      dofsQp(i, 2) = ScalingFactor * pressure;

      dofsQp(i, 6) = ScalingFactor * (kX / (omega * rho)) * std::cos(kX * x) * std::sin(kY * y) *
                     std::cos(omega * t) * (std::sinh(kStar * z) + b * std::cosh(kStar * z));
      dofsQp(i, 7) = ScalingFactor * (kY / (omega * rho)) * std::sin(kX * x) * std::cos(kY * y) *
                     std::cos(omega * t) * (std::sinh(kStar * z) + b * std::cosh(kStar * z));
      dofsQp(i, 8) = ScalingFactor * (kStar / (omega * rho)) * std::sin(kX * x) * std::sin(kY * y) *
                     std::cos(omega * t) * (std::cosh(kStar * z) + b * std::sinh(kStar * z));
    } else {
      // Elastic-acoustic mode
      const auto pressure = -std::sin(kX * x) * std::sin(kY * y) * std::sin(omega * t) *
                            (std::sin(kStar * z) + b * std::cos(kStar * z));
      dofsQp(i, 0) = ScalingFactor * pressure;
      dofsQp(i, 1) = ScalingFactor * pressure;
      dofsQp(i, 2) = ScalingFactor * pressure;
      dofsQp(i, 6) = ScalingFactor * (kX / (omega * rho)) * std::cos(kX * x) * std::sin(kY * y) *
                     std::cos(omega * t) * (std::sin(kStar * z) + b * std::cos(kStar * z));
      dofsQp(i, 7) = ScalingFactor * (kY / (omega * rho)) * std::sin(kX * x) * std::cos(kY * y) *
                     std::cos(omega * t) * (std::sin(kStar * z) + b * std::cos(kStar * z));
      dofsQp(i, 8) = ScalingFactor * (kStar / (omega * rho)) * std::sin(kX * x) * std::sin(kY * y) *
                     std::cos(omega * t) * (std::cos(kStar * z) - b * std::sin(kStar * z));
    }
  }
#else
  dofsQp.setZero();
#endif
}
