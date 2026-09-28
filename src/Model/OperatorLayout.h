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
/// coordinates of that face, and how many nodes of a face carry them. None,
/// where the operator of a face is not those scalars.
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

} // namespace seissol

#endif // SEISSOL_SRC_MODEL_OPERATORLAYOUT_H_
