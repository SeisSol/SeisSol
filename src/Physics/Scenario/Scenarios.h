// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_SCENARIO_SCENARIOS_H_
#define SEISSOL_SRC_PHYSICS_SCENARIO_SCENARIOS_H_

#include "Common/ConfigRegistry.h"
#include "GeneratedCode/init.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include "Initializer/Typedefs.h"
#include "Physics/InitialField.h"

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

namespace seissol::physics {

class PressureInjection : public InitialFieldOf<PressureInjection> {
  public:
  explicit PressureInjection(
      const seissol::initializer::parameters::InitializationParameters& initializationParameters);

  template <typename RealT>
  void evaluateIn(double time,
                  const std::array<double, 3>* points,
                  std::size_t count,
                  const CellMaterialData& materialData,
                  yateto::DenseTensorView<2, RealT, unsigned>& dofsQP) const;

  private:
  seissol::initializer::parameters::InitializationParameters parameters_;
};

// A planar wave travelling in direction kVec
class Planarwave : public InitialFieldOf<Planarwave> {
  public:
  // Choose phase in [0, 2*pi]; the wave is set up for the material of the configuration `config`
  Planarwave(const CellMaterialData& materialData,
             ConfigId config,
             double phase,
             Eigen::Vector3d kVec,
             std::vector<int> varField,
             std::vector<std::complex<double>> ampField);
  Planarwave(const CellMaterialData& materialData,
             ConfigId config,
             double phase = 0.0,
             Eigen::Vector3d kVec = {M_PI, M_PI, M_PI});

  template <typename RealT>
  void evaluateIn(double time,
                  const std::array<double, 3>* points,
                  std::size_t count,
                  const CellMaterialData& materialData,
                  yateto::DenseTensorView<2, RealT, unsigned>& dofsQP) const;

  protected:
  std::vector<int> varField_;
  std::vector<std::complex<double>> ampField_;
  double phase_;
  Eigen::Vector3d kVec_;
  // the quantities of the material, and the eigenvalues and eigenvectors of the plane wave
  // operator for them
  std::size_t numQuantities_{};
  std::vector<std::complex<double>> lambdaA_;
  std::vector<std::complex<double>> eigenvectors_;

  private:
  void init(const CellMaterialData& materialData, ConfigId config);
};

// superimpose three planar waves travelling into different directions
class SuperimposedPlanarwave : public InitialFieldOf<SuperimposedPlanarwave> {
  public:
  //! Choose phase in [0, 2*pi]
  SuperimposedPlanarwave(const CellMaterialData& materialData, ConfigId config, double phase = 0.0);

  template <typename RealT>
  void evaluateIn(double time,
                  const std::array<double, 3>* points,
                  std::size_t count,
                  const CellMaterialData& materialData,
                  yateto::DenseTensorView<2, RealT, unsigned>& dofsQP) const;

  private:
  std::array<Eigen::Vector3d, 3> kVec_;
  std::array<Planarwave, 3> pw_;
};

// A part of a planar wave travelling in one direction
class TravellingWave : public InitialFieldOf<TravellingWave, Planarwave> {
  public:
  TravellingWave(const CellMaterialData& materialData,
                 ConfigId config,
                 const TravellingWaveParameters& travellingWaveParameters);

  template <typename RealT>
  void evaluateIn(double time,
                  const std::array<double, 3>* points,
                  std::size_t count,
                  const CellMaterialData& materialData,
                  yateto::DenseTensorView<2, RealT, unsigned>& dofsQP) const;

  private:
  Eigen::Vector3d origin_;
};

class AcousticTravellingWaveITM : public InitialFieldOf<AcousticTravellingWaveITM> {
  public:
  AcousticTravellingWaveITM(
      const CellMaterialData& materialData,
      const AcousticTravellingWaveParametersITM& acousticTravellingWaveParametersITM);
  template <typename RealT>
  void evaluateIn(double time,
                  const std::array<double, 3>* points,
                  std::size_t count,
                  const CellMaterialData& materialData,
                  yateto::DenseTensorView<2, RealT, unsigned>& dofsQP) const;

  private:
  void init(const CellMaterialData& materialData);
  double rho0_;
  double c0_;
  double k_;
  double tITMMinus_;
  double tau_;
  double tITMPlus_;
  double n_;
};

class ScholteWave : public InitialFieldOf<ScholteWave> {
  public:
  ScholteWave() = default;
  template <typename RealT>
  void evaluateIn(double time,
                  const std::array<double, 3>* points,
                  std::size_t count,
                  const CellMaterialData& materialData,
                  yateto::DenseTensorView<2, RealT, unsigned>& dofsQP) const;
};
class SnellsLaw : public InitialFieldOf<SnellsLaw> {
  public:
  SnellsLaw() = default;
  template <typename RealT>
  void evaluateIn(double time,
                  const std::array<double, 3>* points,
                  std::size_t count,
                  const CellMaterialData& materialData,
                  yateto::DenseTensorView<2, RealT, unsigned>& dofsQP) const;
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
class Ocean : public InitialFieldOf<Ocean> {
  private:
  int mode_;
  double gravitationalAcceleration_;
  // whether the material of the configuration is elastic (and not acoustic), and the index of its
  // first velocity
  bool elastic_;
  std::size_t velocityOffset_;

  public:
  Ocean(int mode, double gravitationalAcceleration, ConfigId config);
  template <typename RealT>
  void evaluateIn(double time,
                  const std::array<double, 3>* points,
                  std::size_t count,
                  const CellMaterialData& materialData,
                  yateto::DenseTensorView<2, RealT, unsigned>& dofsQP) const;
};

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_SCENARIO_SCENARIOS_H_
