// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EXPR_REWRITE_H_
#define SEISSOL_SRC_EXPR_REWRITE_H_

// Rewrites of a Program by its consumer.
//
// A model is written pointwise: it reads `v1`, `s_xx`, `x` or `t` as ordinary channels and knows
// nothing about basis functions. Which of those channels a consumer supplies as a column and which
// it obtains from the degrees of freedom of a cell is the consumer's decision, not the model's --
// so the same program runs against already projected columns (one buffer more) or against the
// modal coefficients themselves (fused, no buffer), and only the binding differs.

#include "Expr/Ir.h"
#include "Expr/Program.h"

#include <functional>
#include <map>
#include <string>
#include <unordered_map>

namespace seissol::expr {

/// Builds the replacement of a channel in the arena of the rewritten program.
using ChannelBuilder = std::function<NodeId(Arena& arena)>;

/// Replaces every read of an input channel named in `builders` by the node its builder makes in
/// the rewritten arena -- an expression over other channels, a contraction, a constant. The
/// replacement is not rewritten again, so a builder may read the channel it replaces.
///
/// State channels are never replaced, and a name the program does not read is not an error.
/// Inputs nothing reads any more are dropped from the signature, and the channels the builders
/// introduce are added to it, in the order they are first created; the result is validated again.
///
/// Throws std::invalid_argument when the rewritten program does not validate.
void substituteChannels(Program& program, const std::map<std::string, ChannelBuilder>& builders);

/// Replaces every read of an input channel named in `blocks` by a contraction of the block mapped
/// to it against `matrix` (cf. Kind::Contract). A mapping, not a list: the assignment is by name.
///
/// A name the program does not read is not an error, and a channel that is not mapped stays an
/// ordinary column; state channels are never replaced. Inputs nothing reads any more are dropped
/// from the signature, which therefore changes along with the fingerprint, and the result is
/// validated again. The matrix and the blocks have to be declared in `program` beforehand
/// (Program::internMatrix, Program::internBlock).
///
/// Throws std::invalid_argument when the rewritten program does not validate, e.g. for a block
/// whose length does not match the columns of `matrix`.
void substituteByContraction(Program& program,
                             MatrixId matrix,
                             const std::map<std::string, BlockId>& blocks);

/// `root` with every read of a channel in `replacements` replaced by the node mapped to it, in the
/// same arena. A subtree that reads none of them stays the node it is, so a program that holds both
/// the original and the replaced expression computes the shared part once -- e.g. one instance of a
/// script per time, in which only what depends on `t` is computed for each. The replacements are
/// not rewritten again. `memo` keeps the rewritten nodes; one memo serves one set of replacements,
/// across any number of roots.
///
/// Throws std::invalid_argument for a node of a kind it does not rewrite (the derived-output
/// kinds, which only the derived-output lowering has).
NodeId replaceChannels(Arena& arena,
                       NodeId root,
                       const std::unordered_map<int, NodeId>& replacements,
                       std::unordered_map<NodeId, NodeId>& memo);

/// The outputs and the states of `second` added to those of `first`, as one program on one DAG:
/// a channel of the same name is one channel, so an input both read is read once, and an
/// expression both compute is computed once. The grids and the contractions of `second` are
/// declared in the result as well. It computes in F64 if either does.
///
/// Throws std::invalid_argument for an output or a state both have, for a state of one that the
/// other reads as an input, and for an input both read with different types.
[[nodiscard]] Program mergePrograms(const Program& first, const Program& second);

} // namespace seissol::expr

#endif // SEISSOL_SRC_EXPR_REWRITE_H_
