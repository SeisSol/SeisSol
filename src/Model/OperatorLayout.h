// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_MODEL_OPERATORLAYOUT_H_
#define SEISSOL_SRC_MODEL_OPERATORLAYOUT_H_

#include "Common/Typedefs.h"
#include "Config.h"
#include "GeneratedCode/coefficients.h"
#include "Initializer/BasicTypedefs.h"

#include <cstddef>

/// How the operator of a cell is stored. The flags sit here rather than beside
/// the cell data, because the fault reads them too and its own declarations
/// cannot include the cell data without a cycle.
namespace seissol {

/// Whether a cell carries the coefficients its operator is linear in together
/// with the rows of its Jacobian, instead of the star matrices the two fold
/// into. The build decides, and a solver that does not declare the
/// decomposition keeps the matrices whatever the build asks for.
constexpr bool FactoredStar = Config::FactoredStar && generated::SolverNumCoefficients > 0;

/// How many scalars the operator of a cell is linear in. The generator states
/// it, since this is needed where a solver's declaration cannot be
/// instantiated; StarOperands.h checks the two against each other.
constexpr std::size_t StarCoefficientCount = generated::SolverNumCoefficients;

/// Whether the material varies inside a cell, so that a cell carries one
/// coefficient per sample point rather than one for itself. It needs the
/// factored star, since the coefficients are what varies.
constexpr bool NodalMaterial = Config::MaterialNodal && FactoredStar;

/// How many samples of the material a cell carries. One, where it does not
/// vary inside the cell.
constexpr std::size_t MaterialSampleCount = generated::MaterialSampleCount;

/// How many points the operator is formed at, one where it is one per cell.
constexpr std::size_t OperatorPointCount = generated::OperatorPointCount;

/// How many scalars the source term of a cell is linear in, and whether it is
/// formed where the material is sampled. A solver without a source term, or
/// one that does not declare its decomposition, carries none.
constexpr std::size_t SourceCoefficientCount = generated::SolverNumSourceCoefficients;
constexpr bool NodalSource = NodalMaterial && SourceCoefficientCount > 0;

/// How many of those a cell carries a second time, as what a sample point
/// deviates from the cell. A solver that factorises the source term into a
/// solve it does once per cell needs them; one that applies it as a product
/// reads the samples themselves and carries none.
constexpr std::size_t SourceDeviationCount = generated::SolverNumSourceDeviations;
constexpr bool NodalSourceDeviation = NodalSource && SourceDeviationCount > 0;

/// How many scalars the flux operator of a face is linear in, in the
/// coordinates of that face, and at how many points of a face -- those of the
/// rule its flux is integrated with -- it carries them. None, where the operator
/// of a face is not those scalars.
constexpr std::size_t FluxCoefficientCount = generated::FluxNumCoefficients;
constexpr std::size_t FluxFaceNodes = generated::FaceNodes;

/// Whether the flux operator varies along a face. A material that varies
/// inside a cell varies along its faces too, but only where the operator there
/// is the handful of scalars above can a face carry it that way; otherwise the
/// face keeps the one operator per side that is built from the cell's material.
constexpr bool NodalFlux = NodalMaterial && FluxCoefficientCount > 0;

/// How many scalars the lift of a fault face reads at each quadrature point of
/// it, and whether it does. The lift is the coefficient matrix of the fault
/// normal, which is a handful of scalars of the face in its own coordinates
/// for every material, so where the material varies inside a cell a fault face
/// carries those at its points; otherwise it carries one matrix per side.
constexpr std::size_t FaultFluxCoefficientCount = generated::FaultFluxNumCoefficients;
constexpr bool NodalFaultFlux = NodalMaterial && FaultFluxCoefficientCount > 0;

/// Whether a cell may be curved. It then carries the rows of its inverse Jacobian at the points
/// the operator is formed at, and a face carries its rotation and its scale at its nodes; which
/// cells are curved, the mesh decides (see geometry::MeshReader::geometryOrder).
constexpr bool Curvilinear = Config::Curvilinear && NodalFlux;
static_assert(Config::Curvilinear == generated::Curvilinear,
              "The build and the generated code disagree about whether a cell may be curved.");
static_assert(!Config::Curvilinear || NodalFlux,
              "A curved face applies its flux at its nodes, which this build does not.");
static_assert(!Config::Curvilinear || NodalFaultFlux,
              "A curved fault face applies its lift at its points, which this build does not.");

/// Whether a cell can lift its fluxes without material and apply the inverse weight of its medium
/// at the sample points afterwards: weight-adjusted DG (Chan, Hewett and Warburton 2017). The
/// operator of a medium is the inverse of its weight times one without material, A = W^-1 B; the
/// faces of such a cell lift their fluxes with B into a buffer of its own (LTS::FluxBuffer), and
/// once all of them are in, the cell adds the buffer read at the sample points, times W^-1 there,
/// projected back (kernel::weightAdjust). Where the material jumps inside a cell, that keeps the
/// scheme stable where lifting with the material at the face does not. It needs the material at
/// the sample points and a medium that states its weight. A medium with memory variables takes the
/// weight of its base material; the rows its faces lift for the memory variables carry no material
/// and pass through. A face whose flux is one matrix per side forms it without material just as a
/// face with its scalars per point does. Which cells take it, the run decides (see
/// seissol::initializer::ImpedanceSpread).
constexpr bool WeightAdjustable = NodalMaterial && generated::WeightAdjustable;

/// Whether the corrector applies its volume term in the strong form. It does wherever the operator
/// varies inside a cell: the weak form lifts the cell's own trace with the operator inside the
/// cell, while the flux states it at the face, and the two only cancel where they are the same
/// operator. In the strong form the corrector of every face subtracts the normal flux of the
/// cell's own trace instead (see initializer::toCorrectorForm, and nodalVolumeDerivatives in the
/// generator).
constexpr bool StrongCorrector = NodalMaterial;

/// Whether the fault takes the subtraction of the strong form on itself: its lift then applies the
/// flux of the fault normal to the imposed state less the side's own trace, at the points of the
/// fault, where it otherwise is the corrector of the face, at the points of the face. Where a cell
/// may be curved, the scale of a point carries the inverse of the Jacobian determinant there, which
/// no quadrature integrates exactly; the lift and the subtraction cancel for a state that is
/// continuous across a locked fault only where they are taken at the same points.
constexpr bool FaultSubtractsOwnTrace = StrongCorrector && Curvilinear;

/// Whether the corrector of a cell applies the flux of a face of this type to the traces of the
/// face. A fault face supplies its flux through the fault, so in the weak form it has none; in the
/// strong form what the corrector applies there is the subtraction of the normal flux of the cell's
/// own trace alone, unless the fault takes that on itself.
constexpr bool appliesFaceCorrector(FaceType faceType) {
  return (StrongCorrector && !FaultSubtractsOwnTrace) || faceType != FaceType::DynamicRupture;
}

} // namespace seissol

#endif // SEISSOL_SRC_MODEL_OPERATORLAYOUT_H_
