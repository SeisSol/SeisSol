// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_INITIALFIELD_H_
#define SEISSOL_SRC_PHYSICS_INITIALFIELD_H_

#include <array>
#include <complex>
#include <vector>

#include <Eigen/Dense>

#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Precision.h"
#include "generated_code/init.h"

namespace seissol::physics {
class InitialField {
  public:
  virtual ~InitialField() = default;
  virtual void evaluate(double time,
                        const std::vector<std::array<double, 3>>& points,
                        const CellMaterialData& materialData,
                        yateto::DenseTensorView<2, real, unsigned>& dofsQP) const = 0;
};

class ZeroField : public InitialField {
  public:
  void evaluate(double,
                const std::vector<std::array<double, 3>>&,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override {
    dofsQP.setZero();
  }
};

class PressureInjection : public InitialField {
  public:
  PressureInjection(
      const seissol::initializer::parameters::InitializationParameters& initializationParameters);

  void evaluate(double,
                const std::vector<std::array<double, 3>>&,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override;

  private:
  seissol::initializer::parameters::InitializationParameters m_parameters;
};

// A planar wave travelling in direction kVec
class Planarwave : public InitialField {
  public:
  // Choose phase in [0, 2*pi]
  Planarwave(const CellMaterialData& materialData,
             double phase,
             Eigen::Vector3d kVec,
             std::vector<int> varField,
             std::vector<std::complex<double>> ampField);
  explicit Planarwave(const CellMaterialData& materialData,
                      double phase = 0.0,
                      Eigen::Vector3d kVec = {M_PI, M_PI, M_PI});

  void evaluate(double time,
                const std::vector<std::array<double, 3>>& points,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override;

  protected:
  std::vector<int> m_varField;
  std::vector<std::complex<double>> m_ampField;
  const double m_phase;
  const Eigen::Vector3d m_kVec;
  std::array<std::complex<double>, seissol::model::MaterialT::NumQuantities> m_lambdaA;
  std::array<std::complex<double>,
             seissol::model::MaterialT::NumQuantities * seissol::model::MaterialT::NumQuantities>
      m_eigenvectors;

  private:
  void init(const CellMaterialData& materialData);
};

// superimpose three planar waves travelling into different directions
class SuperimposedPlanarwave : public InitialField {
  public:
  //! Choose phase in [0, 2*pi]
  SuperimposedPlanarwave(const CellMaterialData& materialData, real phase = 0.0);

  void evaluate(double time,
                const std::vector<std::array<double, 3>>& points,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override;

  private:
  const std::array<Eigen::Vector3d, 3> m_kVec;
  const double m_phase;
  std::array<Planarwave, 3> m_pw;
};

// A part of a planar wave travelling in one direction
class TravellingWave : public Planarwave {
  public:
  TravellingWave(const CellMaterialData& materialData,
                 const TravellingWaveParameters& travellingWaveParameters);

  void evaluate(double time,
                const std::vector<std::array<double, 3>>& points,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override;

  private:
  Eigen::Vector3d m_origin;
};

class AcousticTravellingWaveITM : public InitialField {
  public:
  AcousticTravellingWaveITM(
      const CellMaterialData& materialData,
      const AcousticTravellingWaveParametersITM& acousticTravellingWaveParametersITM);
  void evaluate(double time,
                const std::vector<std::array<double, 3>>& points,
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
                const std::vector<std::array<double, 3>>& points,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override;
};
class SnellsLaw : public InitialField {
  public:
  SnellsLaw() = default;
  void evaluate(double time,
                const std::vector<std::array<double, 3>>& points,
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
                const std::vector<std::array<double, 3>>& points,
                const CellMaterialData& materialData,
                yateto::DenseTensorView<2, real, unsigned>& dofsQP) const override;
};
} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_INITIALFIELD_H_
