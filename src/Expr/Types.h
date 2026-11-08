// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_TYPES_H_
#define SEISSOL_SRC_EXPR_TYPES_H_

// The types of the values of a program (see ValueType).
//
// THE RULES. A value has the compute type of its program (Real), unless the script converts it
// (Kind::Cast: f64(), f32() and real() in sderiv) or it is a state of another type. The operands of
// an operation have one type, which is the type of its value; a literal takes the type of the other
// operands (of the consumer, if it has no others: `f64(2.0)`), and a value of literals alone is
// Real. There is no implicit conversion, also not for a condition: `select(c, a, b)` needs c, a and
// b of one type, as `a * b` needs a and b. So a program is valid or not whatever its compute type
// is, and a value in binary64 does not silently become one in binary32, nor the other way round.
// The coordinates of a lookup have its type; a contraction is Real; a state and the next value of
// it have the type of the state.
//
// WHY EXPLICIT. The point of a value in binary64 in a program in binary32 is precision where the
// script needs it: a time late in a run (binary32 resolves 1.9e-6 s at 20 s), a sum over many
// steps, a difference that cancels. Promotion by the operands would spread binary64 over everything
// that meets such a value, which is the cost the program in binary32 exists to avoid (on many GPUs,
// binary64 runs at 1/32 to 1/64 of the rate); demotion would lose what the conversion was for.
//
// THE DATA. The backends move the data of a program in one type: its inputs and outputs, its
// states and its hoisted values (dataType()). That is binary64 where a value of the program is in
// binary64, and the compute type otherwise. A value of another type converts where it is loaded and
// stored, which is exact where it is in binary32, so the values in the compute type are the same
// bits as in a program without types.

#include "Expr/Ir.h"
#include "Expr/Program.h"

namespace seissol::expr {

/// The type the backends compute a value of type `type` in, in a program whose compute type is
/// `computeType`.
[[nodiscard]] ComputeType resolve(ValueType type, ComputeType computeType);

/// The program, with the type of every value in Node::type, every literal in the type it is used
/// in (one node per type), and the conversions that change nothing gone: a conversion of an input
/// becomes a load of the input in the type converted to. The inputs, outputs, states, grids,
/// matrices and blocks stay as they are. Throws std::invalid_argument for operands of different
/// types, and for a state whose next value has another type than the state.
[[nodiscard]] Program typed(const Program& program);

/// The type in which the backends move the data of `program` (see above).
[[nodiscard]] ComputeType dataType(const Program& program);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_TYPES_H_
