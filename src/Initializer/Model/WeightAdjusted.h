// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_WEIGHTADJUSTED_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_WEIGHTADJUSTED_H_

#include "Common/Constants.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Model/CommonDatastructures.h"
#include "Model/OperatorLayout.h"
#include "Numerical/BasisFunction.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>

namespace seissol::initializer {

/**
 * How much longer than 1 / r times the one its wave speeds permit the time step of a
 * weight-adjusted cell may be at an order of the scheme, r the spread of its impedances (see
 * ImpedanceSpread).
 *
 * On a periodic cube of 1296 cells at CFL 0.5, an acoustic slab whose jumps cut the cells, its
 * impedances by ratios r of 3, 10 and 20, ran through to t = 1 at time steps of up to c / r times
 * the one its wave speeds permit, the smallest c over those ratios 3.0 at order 3, 1.5 at order 4,
 * 1.8 at order 5 and 1.5 at order 6 (r of 3 and 10 there); the smallest c it grew with was 3.4,
 * 1.8, 2.0 and 1.65. An elastic slab by 10 at order 4 ran through up to 2.0. The scale keeps a
 * margin of 1.5 on the c that ran through; an order that was not measured takes the smallest.
 */
constexpr double weightAdjustedTimestepScale(std::size_t order) {
  switch (order) {
  case 3:
    return 2.0;
  case 5:
    return 1.2;
  default:
    // orders 4 and 6, and the ones that were not measured
    return 1.0;
  }
}

/**
 * Which cells are weight-adjusted (see WeightAdjustable), and how much shorter their time step is.
 *
 * A cell is where its impedances spread by more than the ratio the run asks for
 * (ModelParameters::weightAdjustedRatio) over the points its material is asked at -- the samples
 * the operator is formed at and the points of its faces the flux is integrated at -- and they are
 * not close to a material that varies smoothly (ImpedanceJump). Lifting the fluxes with the
 * material at the face grows where the material jumps inside a cell, the more the higher the
 * contrast. In an acoustic prototype it did so, for a jump at one place against the points, at
 * order 4 from an impedance ratio of 4 on, at order 5 at 10, at orders 3 and 6 not up to 10; at
 * order 3 from 5 on where the jump lies elsewhere, and nowhere for a transition over a tenth of a
 * cell.
 *
 * A weight-adjusted cell pays with its time step. Its faces lift the dissipation of the Riemann
 * problem without material, and the inverse weight at the samples scales it with the stiffest
 * material in the cell where the face may see the softest. Its time step is shortened by c / r, r
 * the spread of its impedances and c what the order of the scheme permits
 * (weightAdjustedTimestepScale). An elastic medium carries the same for its shear waves, which are
 * slower: what they ask for counts with the ratio of the speeds.
 */
class ImpedanceSpread {
  public:
  /// One more point of the cell.
  void add(const seissol::model::Material& material) {
    const double density = material.getDensity();
    const double pWave = material.getPWaveSpeed();
    const double sWave = material.getSWaveSpeed();
    include(density * pWave, minimumP_, maximumP_);
    pWaveSpeed_ = std::max(pWaveSpeed_, pWave);
    // a point without shear has no shear impedance to compare
    if (sWave > 0.0) {
      include(density * sWave, minimumS_, maximumS_);
      sWaveSpeed_ = std::max(sWaveSpeed_, sWave);
    }
  }

  /// The largest ratio between two P impedances of the cell, or between two S impedances of its
  /// points that carry shear. One where the material does not vary.
  [[nodiscard]] double ratio() const { return std::max(ratioP(), ratioS()); }

  /// What the time step of a weight-adjusted cell is shortened by against the one its largest wave
  /// speed permits, at the order of the scheme: c / r for the P waves, c as the order permits
  /// (weightAdjustedTimestepScale), and the same for the S waves weighed with their speed against
  /// the largest P wave speed; one at most.
  [[nodiscard]] double timestepFactor(std::size_t order = ConvergenceOrder) const {
    const double pWaves = ratioP();
    const double sWaves = pWaveSpeed_ > 0.0 ? sWaveSpeed_ / pWaveSpeed_ * ratioS() : 0.0;
    return std::min(1.0, weightAdjustedTimestepScale(order) / std::max(pWaves, sWaves));
  }

  private:
  static void include(double impedance, double& minimum, double& maximum) {
    minimum = std::min(minimum, impedance);
    maximum = std::max(maximum, impedance);
  }

  [[nodiscard]] double ratioP() const { return maximumP_ > 0.0 ? maximumP_ / minimumP_ : 1.0; }
  [[nodiscard]] double ratioS() const { return maximumS_ > 0.0 ? maximumS_ / minimumS_ : 1.0; }

  double minimumP_{std::numeric_limits<double>::max()};
  double maximumP_{0.0};
  double minimumS_{std::numeric_limits<double>::max()};
  double maximumS_{0.0};
  double pWaveSpeed_{0.0};
  double sWaveSpeed_{0.0};
};

/**
 * What of the spread of the impedances of a cell a smooth material does not have: the logarithms
 * of its P and of its S impedances at the samples the operator is formed at, their projection onto
 * the polynomials of one degree below the cell's with the weights of the rule the samples are of,
 * and the spread of what the projection leaves, at the samples and, where the flux reads the
 * material at the points of the faces (NodalFlux), at those points as well, as a ratio: one where
 * it leaves nothing.
 *
 * A material that varies smoothly inside the cell is close to such a polynomial, a jump is not,
 * and the plain lift grows where a jump cuts a cell (see ImpedanceSpread). A jump between a face
 * and the samples nearest to it is not among the samples at all, only at the points of the face,
 * where the polynomial the samples give misses it. In the acoustic prototype of the scheme at
 * orders 3 to 5, weight-adjusting the cells where what the projection leaves at the samples and at
 * the points of the faces spreads by more than 2 kept every layer with a jump stable, also one next
 * to a face, which the samples alone missed; it took no cell of a smooth layer that spread by up to
 * 5 inside a cell (docs/weight-adjusted-cells.rst). The projection of the cell's own degree would
 * not tell them apart: it follows the samples of a jump so closely that the cells it leaves out
 * make the plain lift grow.
 */
class ImpedanceJump {
  public:
  /// how many samples the material has in a cell (init::materialNodes)
  static constexpr std::size_t Samples = tensor::materialEval::Shape[0];

  /// how many points of the faces of a cell the material is read at where the flux varies along a
  /// face: those it is integrated at, face by face (FacePointGenerator::fluxPoints)
  static constexpr std::size_t FacePoints = NodalFlux ? Cell::NumFaces * FluxFaceNodes : 0;

  /// The modes of the cell the polynomials up to one degree below the cell's are given in: the
  /// projection onto them with the weights of the rule of the samples (init::materialFit) has a row
  /// per mode, which in a hierarchical and orthogonal basis is zero past the modes of that degree.
  static constexpr std::size_t FitModes = tensor::materialFit::Shape[0];

  /// those modes at the points of the faces
  using FaceBasis = std::array<std::array<double, FitModes>, FacePoints>;

  /// The modes at the points of the faces, from where point `point` lies in the reference cell,
  /// `referencePoint(point)`.
  static FaceBasis
      faceBasis(const std::function<std::array<double, Cell::Dim>(std::size_t)>& referencePoint) {
    FaceBasis basis{};
    for (std::size_t point = 0; point < FacePoints; ++point) {
      const auto xi = referencePoint(point);
      const basisFunction::SampledBasisFunctions<double> sampled(
          ConvergenceOrder, xi[0], xi[1], xi[2]);
      for (std::size_t mode = 0; mode < FitModes; ++mode) {
        basis[point][mode] = sampled.data()[mode];
      }
    }
    return basis;
  }

  explicit ImpedanceJump(const FaceBasis& faceBasis) : faceBasis_(faceBasis) {}

  /// The material at sample `sample` of the cell.
  void add(std::size_t sample, const seissol::model::Material& material) {
    addAt(logP_[sample], logS_[sample], material);
  }

  /// The material at point `point` of the faces of the cell.
  void addFace(std::size_t point, const seissol::model::Material& material) {
    addAt(faceLogP_[point], faceLogS_[point], material);
  }

  /// The ratio the jump accounts for, the larger one of the P and of the S impedances. Where some
  /// points carry shear and others do not, the shear modulus itself jumps, and so does the ratio:
  /// it is unbounded, as it is where an impedance is not a positive number.
  [[nodiscard]] double ratio() const {
    if (withShear_ > 0 && withShear_ < Samples + FacePoints) {
      return std::numeric_limits<double>::infinity();
    }
    const double pRatio = residualRatio(logP_, faceLogP_);
    return withShear_ == Samples + FacePoints ? std::max(pRatio, residualRatio(logS_, faceLogS_))
                                              : pRatio;
  }

  private:
  void addAt(double& logP, double& logS, const seissol::model::Material& material) {
    const double density = material.getDensity();
    logP = std::log(density * material.getPWaveSpeed());
    const double sWave = material.getSWaveSpeed();
    // a point without shear has no shear impedance to fit
    if (sWave > 0.0) {
      logS = std::log(density * sWave);
      ++withShear_;
    }
  }

  [[nodiscard]] double residualRatio(const std::array<double, Samples>& values,
                                     const std::array<double, FacePoints>& faceValues) const {
    const auto project = init::materialFit::view::create(init::materialFit::Values);
    const auto evaluate = init::materialEval::view::create(init::materialEval::Values);
    std::array<double, FitModes> modes{};
    for (std::size_t mode = 0; mode < FitModes; ++mode) {
      for (std::size_t sample = 0; sample < Samples; ++sample) {
        if (project.isInRange(mode, sample)) {
          modes[mode] += project(mode, sample) * values[sample];
        }
      }
    }
    double lowest = std::numeric_limits<double>::max();
    double highest = std::numeric_limits<double>::lowest();
    const auto include = [&](double residual) {
      lowest = std::min(lowest, residual);
      highest = std::max(highest, residual);
    };
    for (std::size_t sample = 0; sample < Samples; ++sample) {
      double residual = values[sample];
      for (std::size_t mode = 0; mode < FitModes; ++mode) {
        if (evaluate.isInRange(sample, mode)) {
          residual -= evaluate(sample, mode) * modes[mode];
        }
      }
      include(residual);
    }
    for (std::size_t point = 0; point < FacePoints; ++point) {
      double residual = faceValues[point];
      for (std::size_t mode = 0; mode < FitModes; ++mode) {
        residual -= faceBasis_[point][mode] * modes[mode];
      }
      include(residual);
    }
    // an impedance that is no positive number makes a logarithm that is no finite one, which the
    // comparisons above pass over
    const double spread = highest - lowest;
    const bool finite =
        std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); }) &&
        std::all_of(
            faceValues.begin(), faceValues.end(), [](double v) { return std::isfinite(v); });
    return finite ? std::exp(spread) : std::numeric_limits<double>::infinity();
  }

  const FaceBasis& faceBasis_;
  std::array<double, Samples> logP_{};
  std::array<double, Samples> logS_{};
  std::array<double, FacePoints> faceLogP_{};
  std::array<double, FacePoints> faceLogS_{};
  std::size_t withShear_{0};
};

/// Whether a cell is weight-adjusted, for the ratios the run asks for (ModelParameters): where its
/// impedances spread by more than `ratioThreshold`, and less their projection still by more than
/// `jumpThreshold` (ImpedanceJump). None is where `ratioThreshold` is zero; with a `jumpThreshold`
/// of one at most, every cell whose impedances spread by more than the ratio is.
inline bool
    isWeightAdjusted(double ratio, double jump, double ratioThreshold, double jumpThreshold) {
  return ratioThreshold > 0.0 && ratio > ratioThreshold &&
         (jumpThreshold <= 1.0 || jump > jumpThreshold);
}

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_WEIGHTADJUSTED_H_
