// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_INITIALFIELD_H_
#define SEISSOL_SRC_PHYSICS_INITIALFIELD_H_

#include <Numerical/Eigenvalues.h>
#include <array>
#include <complex>
#include <vector>

#include <Eigen/Dense>

#include "Equations/Setup.h"
#include "Model/Common.h"

#include "GeneratedCode/init.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"

namespace seissol::physics {
class InitialField {
  public:
  virtual ~InitialField() = default;
  virtual void evaluate(double time,
                        const std::array<double, 3>* points,
                        std::size_t count,
                        const CellMaterialData& materialData,
                        yateto::DenseTensorView<2, real, unsigned>& dofsQP) const = 0;
};

class ZeroField : public InitialField {
  public:
  void evaluate(double time,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override {
    dofsQP.setZero();
  }
};

class PressureInjection : public InitialField {
  public:
  PressureInjection(
      const seissol::initializer::parameters::InitializationParameters& initializationParameters);

  void evaluate(double time,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override;

  private:
  seissol::initializer::parameters::InitializationParameters m_parameters;
};

// A planar wave travelling in direction kVec
template <typename MaterialT>
class Planarwave : public InitialField {
  public:
  // Choose phase in [0, 2*pi]
  Planarwave(const CellMaterialData& materialData,
             double phase,
             Eigen::Vector3d kVec,
             std::vector<int> varField,
             std::vector<std::complex<double>> ampField)
      : m_varField(std::move(varField)), m_ampField(std::move(ampField)), m_phase(phase),
        m_kVec(std::move(kVec)) {
    init(materialData);
  }
  explicit Planarwave(const CellMaterialData& materialData,
                      double phase = 0.0,
                      Eigen::Vector3d kVec = {M_PI, M_PI, M_PI})
      : m_phase(phase), m_kVec(std::move(kVec)) {

    if constexpr (MaterialT::Type == model::MaterialType::Acoustic) {
      // Acoustic materials has the following wave modes:
      // P, N, N, -P
      // Here we impose the P mode
      m_varField = {0};
      m_ampField = {1.0};
    } else if constexpr (MaterialT::Type == model::MaterialType::Poroelastic) {
      // Poroelastic materials have the following wave modes:
      //-P, -S2, -S1, -Ps, N, N, N, N, N, Ps, S1, S2, P
      // Here we impose -S1, -Ps and P
      m_varField = {2, 3, 12};
      m_ampField = {1.0, 1.0, 1.0};
    } else {
      const auto isAcoustic = materialData.local->getMuBar() <= 1e-15;
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
    }
    init(materialData);
  }

  void evaluate(double time,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override {
    dofsQP.setZero();

    auto r = yateto::DenseTensorView<2, std::complex<double>>(
        const_cast<std::complex<double>*>(m_eigenvectors.data()),
        {MaterialT::NumQuantities, MaterialT::NumQuantities});
    for (unsigned v = 0; v < m_varField.size(); ++v) {
      const auto omega = m_lambdaA[m_varField[v]];
      for (unsigned j = 0; j < dofsQP.shape(1); ++j) {
        for (size_t i = 0; i < count; ++i) {
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

  protected:
  std::vector<int> m_varField;
  std::vector<std::complex<double>> m_ampField;
  const double m_phase;
  const Eigen::Vector3d m_kVec;
  std::array<std::complex<double>, MaterialT::NumQuantities> m_lambdaA;
  std::array<std::complex<double>, MaterialT::NumQuantities * MaterialT::NumQuantities>
      m_eigenvectors;

  private:
  void init(const CellMaterialData& materialData) {
    assert(m_varField.size() == m_ampField.size());

    std::array<std::complex<double>, MaterialT::NumQuantities * MaterialT::NumQuantities>
        planeWaveOperator{};
    seissol::model::getPlaneWaveOperator(
        *dynamic_cast<MaterialT*>(materialData.local), m_kVec.data(), planeWaveOperator.data());
    seissol::eigenvalues::Eigenpair<std::complex<double>, MaterialT::NumQuantities>
        eigendecomposition;
    computeEigenvalues(planeWaveOperator, eigendecomposition);
    m_lambdaA = eigendecomposition.values;
    m_eigenvectors = eigendecomposition.vectors;
  }
};

// superimpose three planar waves travelling into different directions
template <typename MaterialT>
class SuperimposedPlanarwave : public InitialField {
  public:
  //! Choose phase in [0, 2*pi]
  SuperimposedPlanarwave(const CellMaterialData& materialData, real phase = 0.0)
      : m_kVec({{{M_PI, 0.0, 0.0}, {0.0, M_PI, 0.0}, {0.0, 0.0, M_PI}}}), m_phase(phase),
        m_pw({Planarwave<MaterialT>(materialData, phase, m_kVec.at(0)),
              Planarwave<MaterialT>(materialData, phase, m_kVec.at(1)),
              Planarwave<MaterialT>(materialData, phase, m_kVec.at(2))}) {}

  void evaluate(double time,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override {
    dofsQP.setZero();

    std::vector<real> dofsPwVector(dofsQP.size());
    auto dofsPW = yateto::DenseTensorView<2, real, unsigned>(dofsPwVector.data(),
                                                             {count, MaterialT::NumQuantities},
                                                             {0, 0},
                                                             {count, MaterialT::NumQuantities});

    for (int pw = 0; pw < 3; pw++) {
      // evaluate each planarwave
      m_pw.at(pw).evaluate(time, points, count, materialData, dofsPW);
      // and add results together
      for (unsigned j = 0; j < dofsQP.shape(1); ++j) {
        for (size_t i = 0; i < count; ++i) {
          dofsQP(i, j) += dofsPW(i, j);
        }
      }
    }
  }

  private:
  const std::array<Eigen::Vector3d, 3> m_kVec;
  const double m_phase;
  std::array<Planarwave<MaterialT>, 3> m_pw;
};

// A part of a planar wave travelling in one direction
template <typename MaterialT>
class TravellingWave : public Planarwave<MaterialT> {
  public:
  TravellingWave(
      const CellMaterialData& materialData,
      const TravellingWaveParameters&
          travellingWaveParameters) // Set phase to 0.5*M_PI, so we have a zero at the origin
      // The wave travels in direction of kVec
      // 2*pi / magnitude(kVec) is the wave length of the wave
      : Planarwave<MaterialT>(materialData,
                              0.5 * M_PI,
                              travellingWaveParameters.kVec,
                              travellingWaveParameters.varField,
                              travellingWaveParameters.ampField),
        // origin is a point on the wavefront at time zero
        m_origin(travellingWaveParameters.origin) {
    logInfo() << "Impose a travelling wave as initial condition";
    logInfo() << "Origin = (" << m_origin[0] << ", " << m_origin[1] << ", " << m_origin[2] << ")";
    logInfo() << "kVec = (" << this->m_kVec[0] << ", " << this->m_kVec[1] << ", " << this->m_kVec[2]
              << ")";
    logInfo() << "Combine following wave modes";
    for (size_t i = 0; i < this->m_ampField.size(); i++) {
      logInfo() << "(" << this->m_varField[i] << ": " << this->m_ampField[i] << ")";
    }
  }

  void evaluate(double time,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQp) const override {
    dofsQp.setZero();

    auto r = yateto::DenseTensorView<2, std::complex<double>>(
        const_cast<std::complex<double>*>(this->m_eigenvectors.data()),
        {MaterialT::NumQuantities, MaterialT::NumQuantities});
    for (unsigned v = 0; v < this->m_varField.size(); ++v) {
      const auto omega = this->m_lambdaA[this->m_varField[v]];
      for (unsigned j = 0; j < dofsQp.shape(1); ++j) {
        for (size_t i = 0; i < count; ++i) {
          auto arg = std::complex<double>(0.0, 1.0) *
                     (omega * time - this->m_kVec[0] * (points[i][0] - m_origin[0]) -
                      this->m_kVec[1] * (points[i][1] - m_origin[1]) -
                      this->m_kVec[2] * (points[i][2] - m_origin[2]) + this->m_phase);
          if (arg.imag() > -0.5 * M_PI && arg.imag() < 1.5 * M_PI) {
            dofsQp(i, j) +=
                (r(j, this->m_varField[v]) * this->m_ampField[v] * std::exp(arg)).real();
          }
        }
      }
    }
  }

  private:
  Eigen::Vector3d m_origin;
};

class AcousticTravellingWaveITM : public InitialField {
  public:
  AcousticTravellingWaveITM(
      const CellMaterialData& materialData,
      const AcousticTravellingWaveParametersITM& acousticTravellingWaveParametersITM);
  void evaluate(double time,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override;

  private:
  void init(const CellMaterialData& materialData);
  double rho0;
  double c0;
  double k;
  double tITMMinus;
  double tau;
  double tITMPlus;
  double n;
};

class ScholteWave : public InitialField {
  public:
  ScholteWave() = default;
  void evaluate(double time,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override;
};
class SnellsLaw : public InitialField {
  public:
  SnellsLaw() = default;
  void evaluate(double time,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override;
};
/*
 * From
 * Abrahams, L. S., Krenz, L., Dunham, E. M., & Gabriel, A. A. (2019, December).
 * Verification of a 3D fully-coupled earthquake and tsunami model.
 * In AGU Fall Meeting Abstracts (Vol. 2019, pp. NH43F-1000).
 * A 3D extension of the 2D scenario in
 * Lotto, G. C., & Dunham, E. M. (2015).
 * High-order finite difference modeling of tsunami generation in a compressible ocean from offshore
 * earthquakes. Computational Geosciences, 19(2), 327-340.
 */
class Ocean : public InitialField {
  private:
  const int mode;
  const double gravitationalAcceleration;

  public:
  Ocean(int mode, double gravitationalAcceleration);
  void evaluate(double time,
                const std::array<double, 3>* points,
                std::size_t count,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override;
};
} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_INITIALFIELD_H_
