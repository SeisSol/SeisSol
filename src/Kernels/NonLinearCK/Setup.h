// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_NONLINEARCK_SETUP_H_
#define SEISSOL_SRC_KERNELS_NONLINEARCK_SETUP_H_

// IWYU pragma: begin_exports

#include "Kernels/NonLinearCK/Local.h"
#include "Kernels/NonLinearCK/Neighbor.h"
#include "Kernels/NonLinearCK/Time.h"

// IWYU pragma: end_exports

#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/quantities.h"
#include "Initializer/BasicTypedefs.h"
#include "Kernels/NonLinearCK/Solver.h"
#include "Model/Common.h"

#include <cstddef>
#include <string>
#include <string_view>
#include <utils/logger.h>

namespace seissol::model {

/// What a dissipation matrix is written through.
// (note: the nullptr cast is needed to differentiate const vs non-const view)
using ProjectorView = decltype(init::fluxDissipation::view::create(static_cast<real*>(nullptr)));

/**
 * The initial strain is a material parameter that the kernels read as a
 * tensor, so it is converted once, here, rather than at every timestep: the
 * material holds it in double and the kernels work in the solver's precision.
 */
template <typename MaterialT>
struct SolverSetup<kernels::solver::nonlinearck::Solver, MaterialT>
    : public SolverSetupDefaults<kernels::solver::nonlinearck::Solver, MaterialT> {
  static void initializeSpecificLocalData(const MaterialT& material,
                                          double /*timeStepWidth*/,
                                          typename MaterialT::Solver::LocalData* localData) {
    auto epsInit = init::epsInit::view::create(localData->epsInit);
    epsInit.setZero();
    epsInit(0) = material.epsInitXX;
    epsInit(1) = material.epsInitYY;
    epsInit(2) = material.epsInitZZ;
    epsInit(3) = material.epsInitXY;
    epsInit(4) = material.epsInitYZ;
    epsInit(5) = material.epsInitXZ;

    localData->maxWaveSpeedBound = static_cast<real>(material.getMaxWaveSpeed());

    // Filled by name against the order the codegen chose, so the two cannot
    // drift apart: a parameter in the wrong slot would be a different
    // material, silently.
    const auto set = [&](std::string_view name, double value) {
      for (std::size_t i = 0; i < generated::MaterialParameterNames.size(); ++i) {
        if (generated::MaterialParameterNames[i] == name) {
          localData->parameters[i] = static_cast<real>(value);
          return;
        }
      }
      logError() << "The kernels do not read a material parameter called" << name.data();
    };

    set("rhoInv", 1.0 / material.rho);
    set("lambda0", material.lambda0);
    set("mu0", material.mu0);
    set("gammaR", material.gammaR);
    set("xi0", material.xi0);
    set("damageRate", material.damageRate);
    set("breakageRate", material.breakageRate);
    set("healingRate", material.healingRate);
    set("betaAlpha", material.betaAlpha);
    for (std::size_t i = 0; i < material.aB.size(); ++i) {
      set("aB" + std::to_string(i), material.aB.at(i));
    }
  }

  /// The pair of matrices a face applies, built from the flux of the face
  /// normal rather than from a Godunov state. A Godunov state is a state, so
  /// the solver it builds maps the state onto itself; what a cell transports
  /// here is wider than its state, so the flux of the face normal is
  /// assembled instead.
  /// One entry of a wave family in the face frame: which state row it feeds,
  /// which transported quantity feeds it, and with what weight.
  struct FamilyEntry {
    std::size_t row;
    std::size_t column;
    double weight;
  };

  /// The two matrices a face frame splits its coupled quantities into, rotated
  /// into the frame the flux is applied in.
  ///
  /// Together they are the absolute value of the flux Jacobian of the face
  /// normal, which is what an upwind flux dissipates with. For a system whose
  /// state is a strain that is
  ///
  ///   |A| = diag( S Q^(-1/2) T , Q^(1/2) ),
  ///
  /// with Q the acoustic tensor of the face normal, T the map from a strain to
  /// the traction it carries over the density, and S the lift of a vector back
  /// onto a strain. The point of writing it this way is that T of the state is
  /// the traction, and the traction is transported: the strain rows read the
  /// stress a cell hands over, not its own strain, so no material crosses a
  /// face here either. In the face frame it reads
  ///
  ///   eps_nn  <- sigma_nn  / (rho c_p),   v_n <- c_p v_n,
  ///   eps_nt  <- sigma_nt  / (2 rho c_s), v_t <- c_s v_t,
  ///
  /// and the three strain components a face normal does not transport are fed
  /// by nothing, which is what their zero eigenvalue says.
  ///
  /// Reading the strain instead would not merely be a different amount of
  /// dissipation. Dissipativity is M |A| being symmetric and positive
  /// semi-definite for the energy metric M = diag(W C, rho), and a family
  /// written as an index projector on the strain fails both: it has
  /// eigenvalues of either sign, so a face built that way feeds the modes it
  /// is there to damp. The form above is symmetric and semi-definite by
  /// construction -- the strain part of the rate is rho tau . Q^(-1/2) tau
  /// with tau the traction jump over the density -- and it stays so for any
  /// positive pair of speeds, which is what keeps it safe once the damage has
  /// made the tangent anisotropic and the two families only a surrogate.
  ///
  /// The weights carry the undamaged moduli, so that a scalar per family is
  /// still all a step has to supply: the speeds the kernel multiplies with are
  /// the damaged ones, and rho c^2 of the undamaged material is what turns
  /// them back into the reciprocal the strain rows want. That is exact in the
  /// elastic limit and errs towards more dissipation away from it.
  static void writeUpwindDissipation(const double* normal,
                                     const double* tangent1,
                                     const double* tangent2,
                                     double scale,
                                     const MaterialT& material,
                                     ProjectorView& pressure,
                                     ProjectorView& shear) {
    real rotationData[tensor::ghostMap::size()]{};
    real inverseData[tensor::ghostMap::size()]{};
    auto rotation = init::ghostMap::view::create(rotationData);
    auto inverse = init::ghostMap::view::create(inverseData);
    rotation.setZero();
    inverse.setZero();
    model::detail::writeRotationBlocks<false>(
        MaterialT::TransportGroups, normal, tangent1, tangent2, rotation);
    model::detail::writeRotationBlocks<true>(
        MaterialT::TransportGroups, normal, tangent1, tangent2, inverse);

    const double lambda = material.getLambdaBar();
    const double mu = material.getMuBar();

    // Where each group sits: the strain and the velocity of the state, and the
    // stress and the velocity of what is transported. Asked for by role, so a
    // layout that moves a group does not silently move a dissipation with it.
    constexpr std::size_t StrainRow = roleOffset(MaterialT::PrimaryGroups, FaceRole::Traction);
    constexpr std::size_t StrainExtent = roleExtent(MaterialT::PrimaryGroups, FaceRole::Traction);
    constexpr std::size_t VelocityRow = roleOffset(MaterialT::PrimaryGroups, FaceRole::Velocity);
    constexpr std::size_t VelocityExtent = roleExtent(MaterialT::PrimaryGroups, FaceRole::Velocity);
    constexpr std::size_t StressColumn = roleOffset(MaterialT::TransportGroups, FaceRole::Traction);
    constexpr std::size_t VelocityColumn =
        roleOffset(MaterialT::TransportGroups, FaceRole::Velocity);

    // Voigt order, so the normal component of a face frame is 0 and the two
    // shear components naming it are 3 and 5.
    const std::array<FamilyEntry, 2> Pressure{
        {{StrainRow + 0, StressColumn + 0, 1.0 / (lambda + 2.0 * mu)},
         {VelocityRow + 0, VelocityColumn + 0, 1.0}}};
    const std::array<FamilyEntry, 4> Shear{{{StrainRow + 3, StressColumn + 3, 0.5 / mu},
                                            {StrainRow + 5, StressColumn + 5, 0.5 / mu},
                                            {VelocityRow + 1, VelocityColumn + 1, 1.0},
                                            {VelocityRow + 2, VelocityColumn + 2, 1.0}}};

    // Per pair of groups rather than over the block: a rotation mixes the rows
    // of a group among themselves and no further, so the conjugation is zero
    // between any two that this does not name. The matrix is indexed by what
    // feeds it first and what it feeds second, which is the order the flux
    // kernel contracts it in.
    const auto conjugate = [&](const auto& family,
                               std::size_t rowBegin,
                               std::size_t rowExtent,
                               std::size_t columnBegin,
                               std::size_t columnExtent,
                               auto& target) {
      for (std::size_t row = rowBegin; row < rowBegin + rowExtent; ++row) {
        for (std::size_t column = columnBegin; column < columnBegin + columnExtent; ++column) {
          real sum = 0.0;
          for (const auto& entry : family) {
            const bool inBlock = entry.row >= rowBegin && entry.row < rowBegin + rowExtent &&
                                 entry.column >= columnBegin &&
                                 entry.column < columnBegin + columnExtent;
            if (inBlock) {
              sum += rotation(row, entry.row) * entry.weight * inverse(entry.column, column);
            }
          }
          target(column, row) = scale * sum;
        }
      }
    };
    conjugate(Pressure, StrainRow, StrainExtent, StressColumn, StrainExtent, pressure);
    conjugate(Pressure, VelocityRow, VelocityExtent, VelocityColumn, VelocityExtent, pressure);
    conjugate(Shear, StrainRow, StrainExtent, StressColumn, StrainExtent, shear);
    conjugate(Shear, VelocityRow, VelocityExtent, VelocityColumn, VelocityExtent, shear);
  }

  static void assembleTabulatedFaceFlux(FaceType faceType,
                                        std::size_t /*side*/,
                                        double surface,
                                        double volume,
                                        const double* normal,
                                        const double* tangent1,
                                        const double* tangent2,
                                        const MaterialT& materialLocal,
                                        real* aPlusT,
                                        real* aMinusT,
                                        real* aMinusTShear,
                                        bool upwind) {
    static_assert(tensor::fluxConstant::size() == tensor::AplusT::size(),
                  "The constant half of the flux solver is stored in the slot of the "
                  "solver it is half of.");

    real normalData[3];
    for (std::size_t i = 0; i < 3; ++i) {
      normalData[i] = static_cast<real>(normal[i]);
    }

    // A face without a neighbour has a ghost rule, and the rule folds
    // into the pair: outflow is the local state on both sides, so its
    // average is the local flux and its jump is nothing. Everything
    // else that has no neighbour needs a mirror, and the mirror needs
    // the rotation of the transported quantities.
    const bool outflow = faceType == FaceType::Outflow;
    const bool freeSurface = faceType == FaceType::FreeSurface;

    // A fault carries its own flux, imposed by the friction it is
    // under, so the pair of a rupture face is zero on both matrices
    // and the face contributes nothing from here. The linear solver
    // says the same thing by scaling its flux solver with zero.
    const bool rupture = faceType == FaceType::DynamicRupture;

    if (faceType != FaceType::Regular && !outflow && !freeSurface && !rupture) {
      logError() << "The nonlinear solver has no ghost rule for face type"
                 << static_cast<int>(faceType) << "yet.";
    }

    kernel::damageFluxSolver fluxSolver;
    // Scale with |S_side|/|J|, negated because the flux matrices are
    // subtracted -- as the linear path does a few lines further down.
    // An outflow face applies its own half twice, since the ghost
    // state is its own and there is no second half to come.
    const double faceScale = -2.0 * surface / (6.0 * volume);
    fluxSolver.fluxScale = rupture ? 0.0 : (outflow ? 2.0 : 1.0) * faceScale;
    fluxSolver.rhoInv = 1.0 / materialLocal.rho;
    fluxSolver.faceNormal = normalData;
    fluxSolver.fluxConstant = aPlusT;
    fluxSolver.bindGlobals(Pool::host());
    fluxSolver.execute();

    // The dissipation of a face with a neighbour is the identity on
    // the quantities the two cells couple through, scaled the way the
    // flux is; an outflow face dissipates nothing, because there is no
    // jump to dissipate.
    // Which flux this is, is what stands in these two and nothing else: the
    // kernel adds a bound times each of them, and never learns which it has.
    // The identity in the first and nothing in the second means the first
    // bound scales every mode -- and the first bound is the larger of the
    // two, so that is Rusanov. The two projectors mean each family is scaled
    // with its own speed, and the modes that do not propagate with neither.
    auto dissipation = init::fluxDissipation::view::create(aMinusT);
    auto shear = init::fluxDissipation::view::create(aMinusTShear);
    dissipation.setZero();
    shear.setZero();
    if (!outflow && !rupture) {
      if (upwind) {
        writeUpwindDissipation(
            normal, tangent1, tangent2, 0.5 * faceScale, materialLocal, dissipation, shear);
      } else {
        for (std::size_t row = 0; row < generated::CoupledQuantities; ++row) {
          dissipation(row, row) = 0.5 * faceScale;
        }
      }
    }

    if (freeSurface) {
      // The traction of the face has to vanish, and with the stress
      // transported the condition is that and nothing else: a ghost
      // state whose traction is the negated local one averages to
      // zero traction, and its velocity is the local one, so the
      // jump the dissipation sees is the traction alone.
      //
      // Negating a traction is a reflection in the face-local frame,
      // so the map is the rotation there, the signs, and the rotation
      // back. It is folded into the pair here, which is why a free
      // surface costs a regular face's arithmetic afterwards.
      real toGlobalData[tensor::ghostMap::size()]{};
      real toFaceData[tensor::ghostMap::size()]{};
      auto toGlobal = init::ghostMap::view::create(toGlobalData);
      auto toFace = init::ghostMap::view::create(toFaceData);
      toGlobal.setZero();
      toFace.setZero();
      model::detail::writeRotationBlocks<false>(
          model::MaterialT::TransportGroups, normal, tangent1, tangent2, toGlobal);
      model::detail::writeRotationBlocks<true>(
          model::MaterialT::TransportGroups, normal, tangent1, tangent2, toFace);

      std::array<real, tensor::ghostMap::Shape[0]> mirror{};
      mirror.fill(1.0);
      std::size_t offset = 0;
      for (const auto& group : model::MaterialT::TransportGroups) {
        if (group.kind == model::QuantityKind::SymTensor2) {
          for (const auto component : model::SymTensor2Traction) {
            mirror[offset + component] = -1.0;
          }
        }
        offset += group.extent();
      }

      real ghostData[tensor::ghostMap::size()]{};
      auto ghost = init::ghostMap::view::create(ghostData);
      ghost.setZero();
      for (std::size_t row = 0; row < tensor::ghostMap::Shape[0]; ++row) {
        for (std::size_t column = 0; column < tensor::ghostMap::Shape[1]; ++column) {
          real sum = 0.0;
          for (std::size_t k = 0; k < tensor::ghostMap::Shape[0]; ++k) {
            sum += toGlobal(row, k) * mirror[k] * toFace(k, column);
          }
          ghost(row, column) = sum;
        }
      }

      real foldedData[tensor::fluxFolded::size()]{};
      kernel::damageFluxGhost fold;
      fold.ghostMap = ghostData;
      fold.fluxFolded = foldedData;
      // Every matrix the two halves are assembled from, with the sign the
      // ghost state enters that half with: the constant part adds its mirror,
      // each dissipation subtracts it, and the halves then read
      // C + G^T C +- lambda (D - G^T D). One matrix per wave family, so the
      // families are folded one by one -- a family left out of this would
      // scale the local trace of the face against nothing.
      for (const auto& [sign, slot] : {std::pair<double, real*>{1.0, aPlusT},
                                       std::pair<double, real*>{-1.0, aMinusT},
                                       std::pair<double, real*>{-1.0, aMinusTShear}}) {
        fold.ghostSign = sign;
        fold.fluxSource = slot;
        fold.execute();
        std::copy_n(foldedData, tensor::fluxFolded::size(), slot);
      }
    }
  }

  /// The flux of what the friction imposed, rather than the flux of a state:
  /// the tables the kernel is built from are constants of the pool, and the
  /// density is the only thing per cell it needs.
  template <typename KernelT>
  static void bindFaultFluxSolver(KernelT& krnl, const MaterialT& material, const real* /*star*/) {
    krnl.rhoInv = 1.0 / material.rho;
  }

  /// The second expansion, in its own basis: the fused projection sums both,
  /// because what a face reads is both.
  template <typename KernelT>
  static void bindFaultTimeCoefficient(KernelT& krnl,
                                       std::size_t index,
                                       const kernels::TimeCoefficients& coeffs,
                                       std::size_t power) {
    krnl.coeffDR(index) = coeffs.state[power];
    krnl.extraCoeffDR(index) = coeffs.extra[power];
  }

  static void
      initializeSpecificNeighborData(const MaterialT& material,
                                     typename MaterialT::Solver::NeighborData* neighborData) {
    neighborData->maxWaveSpeedBound.fill(static_cast<real>(material.getMaxWaveSpeed()));
  }
};

} // namespace seissol::model

#endif // SEISSOL_SRC_KERNELS_NONLINEARCK_SETUP_H_
