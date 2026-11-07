// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "Expr/Ir.h"
#include "Expr/Lower.h"
#include "Expr/Program.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::expr::test {

namespace {

using DT = reader::scripting::DataType;

std::uint64_t bitsOf(double v) {
  std::uint64_t u = 0;
  std::memcpy(&u, &v, sizeof(u));
  return u;
}

int countOpcode(const StageCode& stage, Opcode op) {
  return static_cast<int>(std::count_if(
      stage.code.begin(), stage.code.end(), [op](const Instruction& i) { return i.op == op; }));
}

} // namespace

TEST_SUITE("Expr::Ir") {

  TEST_CASE("interning collapses structurally equal subtrees") {
    Arena arena;
    const NodeId x = arena.field("x");
    const NodeId y = arena.field("y");
    const NodeId a = arena.pw(Fn::Mul, arena.pw(Fn::Add, x, y), arena.konst(2.0));
    const NodeId b = arena.pw(Fn::Mul, arena.pw(Fn::Add, x, y), arena.konst(2.0));
    CHECK(a == b);
    // x, y, add, 2.0, mul
    CHECK(arena.size() == 5);
  }

  // The sderiv NodeEq compares Const by `==`, which is not reflexive for NaN, so
  // a NaN key is inserted and then never found again: every konst(NaN) allocates.
  TEST_CASE("NaN constants intern to one node") {
    Arena arena;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const NodeId first = arena.konst(nan);
    for (int i = 0; i < 16; ++i) {
      CHECK(arena.konst(nan) == first);
    }
    CHECK(arena.size() == 1);
  }

  // The other direction of the same bug: `+0.0 == -0.0` is true, so the two
  // collapse and whichever lost the race silently changes the sign of 1/x.
  TEST_CASE("signed zeros stay distinct") {
    Arena arena;
    const NodeId plus = arena.konst(0.0);
    const NodeId minus = arena.konst(-0.0);
    CHECK(plus != minus);
    CHECK(arena.size() == 2);
    CHECK(bitsOf(arena[plus].value) != bitsOf(arena[minus].value));
    CHECK(std::isinf(1.0 / arena[plus].value));
    CHECK(1.0 / arena[minus].value < 0.0);
  }

  TEST_CASE("lookup interning covers the coordinate span") {
    Arena arena;
    const NodeId x = arena.field("x");
    const NodeId y = arena.field("y");
    const NodeId z = arena.field("z");

    const NodeId l0 = arena.lookup(0, 0, {x, y, z});
    const NodeId l1 = arena.lookup(0, 0, {x, y, z});
    CHECK(l0 == l1); // same span -> same node

    const NodeId l2 = arena.lookup(0, 0, {x, z, y}); // permuted coordinates
    const NodeId l3 = arena.lookup(0, 1, {x, y, z}); // different component
    const NodeId l4 = arena.lookup(0, 0, {x, y});    // different arity
    CHECK(l2 != l0);
    CHECK(l3 != l0);
    CHECK(l4 != l0);

    const std::vector<NodeId> kids = arena.children(l0);
    CHECK(kids == std::vector<NodeId>{x, y, z});
  }

  // NodeHash/NodeEq hold a back-pointer to the owning Arena so they can reach
  // args_ for Lookup nodes. A defaulted move would carry the pool over with the
  // pointer still aimed at the moved-from arena.
  TEST_CASE("arena survives move and copy") {
    Arena source;
    const NodeId x = source.field("x");
    const NodeId y = source.field("y");
    const NodeId look = source.lookup(0, 0, {x, y});
    const std::size_t sizeBefore = source.size();

    Arena moved(std::move(source));
    CHECK(moved.size() == sizeBefore);
    CHECK(moved.lookup(0, 0, {x, y}) == look); // must still dedup
    CHECK(moved.size() == sizeBefore);
    CHECK(moved.channelName(moved[x].ch) == "x");

    Arena copied(moved);
    CHECK(copied.lookup(0, 0, {x, y}) == look);
    CHECK(copied.size() == sizeBefore);

    Arena assigned;
    assigned = copied;
    CHECK(assigned.lookup(0, 0, {x, y}) == look);
    CHECK(assigned.size() == sizeBefore);
  }

  TEST_CASE("property: independent construction of the same tree yields one node") {
    std::mt19937 rng(20260821);
    for (int trial = 0; trial < 200; ++trial) {
      Arena arena;
      const std::vector<NodeId> leaves{arena.field("x"),
                                       arena.field("y"),
                                       arena.konst(0.0),
                                       arena.konst(-0.0),
                                       arena.konst(std::numeric_limits<double>::quiet_NaN())};
      const std::size_t leafCount = arena.size();

      // Two independent builds of the same random tree must return the same id
      // and must not grow the arena the second time round.
      std::vector<std::uint32_t> script(24);
      for (auto& s : script) {
        s = rng();
      }
      std::size_t cursor = 0;
      const auto build = [&](auto&& self, int depth) -> NodeId {
        const std::uint32_t token = script[cursor++ % script.size()];
        if (depth == 0) {
          return leaves[token % leaves.size()];
        }
        const NodeId a = self(self, depth - 1);
        const NodeId b = self(self, depth - 1);
        static constexpr Fn Ops[] = {Fn::Add, Fn::Mul, Fn::Min, Fn::Sub};
        return arena.pw(Ops[token % 4], a, b);
      };

      cursor = 0;
      const NodeId first = build(build, 3);
      const std::size_t afterFirst = arena.size();
      cursor = 0;
      const NodeId second = build(build, 3);

      CHECK(first == second);
      CHECK(arena.size() == afterFirst);
      CHECK(afterFirst >= leafCount);
    }
  }

  TEST_CASE("name tables round-trip") {
    for (int i = 0; i <= static_cast<int>(Fn::Select); ++i) {
      const auto f = static_cast<Fn>(i);
      Fn back{};
      REQUIRE(fnFromName(name(f), back));
      CHECK(back == f);
    }
    for (int i = 0; i <= static_cast<int>(Red::Last); ++i) {
      const auto r = static_cast<Red>(i);
      Red back{};
      REQUIRE(redFromName(name(r), back));
      CHECK(back == r);
    }
    Fn unknown{};
    CHECK_FALSE(fnFromName("definitely-not-an-op", unknown));
  }
}

TEST_SUITE("Expr::Validate") {

  TEST_CASE("a well-formed pointwise program validates") {
    Program p;
    p.addInput("x", DT::F64);
    p.addInput("t", DT::F64);
    Arena& a = p.arena();
    p.addOutput("u", DT::F64, a.pw(Fn::Mul, a.field("x"), a.pw(Fn::Sin, a.field("t"))));
    CHECK_NOTHROW(validate(p));
  }

  TEST_CASE("stateful and element-context kinds are rejected") {
    SUBCASE("cumint") {
      Program p;
      p.addInput("x", DT::F64);
      Arena& a = p.arena();
      p.addOutput("u", DT::F64, a.cumint(a.field("x")));
      CHECK_THROWS_AS(validate(p), std::invalid_argument);
    }
    SUBCASE("fold") {
      Program p;
      p.addInput("x", DT::F64);
      Arena& a = p.arena();
      p.addOutput("u", DT::F64, a.fold(Red::Max, a.field("x")));
      CHECK_THROWS_AS(validate(p), std::invalid_argument);
    }
    SUBCASE("dx") {
      Program p;
      p.addInput("x", DT::F64);
      Arena& a = p.arena();
      p.addOutput("u", DT::F64, a.dx(0, a.field("x")));
      CHECK_THROWS_AS(validate(p), std::invalid_argument);
    }
    SUBCASE("sample") {
      Program p;
      p.addInput("x", DT::F64);
      Arena& a = p.arena();
      p.addOutput("u", DT::F64, a.sample(a.field("x")));
      CHECK_THROWS_AS(validate(p), std::invalid_argument);
    }
  }

  TEST_CASE("an undeclared channel is rejected") {
    Program p;
    p.addInput("x", DT::F64);
    Arena& a = p.arena();
    p.addOutput("u", DT::F64, a.pw(Fn::Add, a.field("x"), a.field("undeclared")));
    CHECK_THROWS_AS(validate(p), std::invalid_argument);
  }

  TEST_CASE("state is an accepted channel source") {
    Program p;
    p.addInput("f", DT::F64);
    Arena& a = p.arena();
    const NodeId next = a.pw(Fn::Add, a.field("acc"), a.field("f"));
    p.addState("acc", 0.0, next);
    p.addOutput("u", DT::F64, a.field("acc"));
    CHECK_NOTHROW(validate(p));
  }

  TEST_CASE("a name cannot be both input and state") {
    Program p;
    p.addInput("acc", DT::F64);
    Arena& a = p.arena();
    p.addState("acc", 0.0, a.field("acc"));
    p.addOutput("u", DT::F64, a.field("acc"));
    CHECK_THROWS_AS(validate(p), std::invalid_argument);
  }

  TEST_CASE("an out-of-range grid id is rejected") {
    Program p;
    p.addInput("x", DT::F64);
    Arena& a = p.arena();
    p.addOutput("u", DT::F64, a.lookup(3, 0, {a.field("x")}));
    CHECK_THROWS_AS(validate(p), std::invalid_argument);
  }

  TEST_CASE("duplicate names are rejected") {
    Program p;
    p.addInput("x", DT::F64);
    p.addInput("x", DT::F64);
    Arena& a = p.arena();
    p.addOutput("u", DT::F64, a.field("x"));
    CHECK_THROWS_AS(validate(p), std::invalid_argument);
  }
}

TEST_SUITE("Expr::Fingerprint") {

  TEST_CASE("the fingerprint follows DAG shape, not construction order") {
    const auto make = [](bool reversed) {
      Program p;
      p.addInput("x", DT::F64);
      p.addInput("y", DT::F64);
      Arena& a = p.arena();
      NodeId x = NoNode;
      NodeId y = NoNode;
      if (reversed) {
        y = a.field("y");
        x = a.field("x");
      } else {
        x = a.field("x");
        y = a.field("y");
      }
      const NodeId mul = a.pw(Fn::Mul, x, y);
      const NodeId add = a.pw(Fn::Add, x, y);
      p.addOutput("u", DT::F64, a.pw(Fn::Sub, mul, add));
      return p;
    };
    CHECK(make(false).canonicalForm() == make(true).canonicalForm());
    CHECK(make(false).fingerprint() == make(true).fingerprint());
  }

  TEST_CASE("channel names do not reach the fingerprint, slot order does") {
    // Renaming both channels leaves the gather layout alone: sub(in0, in1) in
    // both cases, so it is one kernel.
    const auto renamed = [](const std::string& first, const std::string& second) {
      Program p;
      p.addInput(first, DT::F64);
      p.addInput(second, DT::F64);
      Arena& a = p.arena();
      p.addOutput("u", DT::F64, a.pw(Fn::Sub, a.field(first), a.field(second)));
      return p;
    };
    CHECK(renamed("x", "y").fingerprint() == renamed("alpha", "beta").fingerprint());

    // Declaring the same two inputs the other way round while keeping the
    // expression fixed does change the layout: sub(in0, in1) vs sub(in1, in0).
    const auto reordered = [](bool swap) {
      Program p;
      if (swap) {
        p.addInput("y", DT::F64);
        p.addInput("x", DT::F64);
      } else {
        p.addInput("x", DT::F64);
        p.addInput("y", DT::F64);
      }
      Arena& a = p.arena();
      p.addOutput("u", DT::F64, a.pw(Fn::Sub, a.field("x"), a.field("y")));
      return p;
    };
    CHECK(reordered(false).fingerprint() != reordered(true).fingerprint());
  }

  TEST_CASE("compute type and constant bit patterns are part of the identity") {
    const auto make = [](ComputeType ct, double c) {
      Program p;
      p.addInput("x", DT::F64);
      Arena& a = p.arena();
      p.setComputeType(ct);
      p.addOutput("u", DT::F64, a.pw(Fn::Add, a.field("x"), a.konst(c)));
      return p;
    };
    CHECK(make(ComputeType::F64, 1.0).fingerprint() != make(ComputeType::F32, 1.0).fingerprint());
    CHECK(make(ComputeType::F64, 0.0).fingerprint() != make(ComputeType::F64, -0.0).fingerprint());
    CHECK(make(ComputeType::F64, 1.0).fingerprint() == make(ComputeType::F64, 1.0).fingerprint());
  }

  TEST_CASE("output order is part of the identity") {
    const auto make = [](bool swapped) {
      Program p;
      p.addInput("x", DT::F64);
      Arena& a = p.arena();
      const NodeId s = a.pw(Fn::Sin, a.field("x"));
      const NodeId c = a.pw(Fn::Cos, a.field("x"));
      p.addOutput("u", DT::F64, swapped ? c : s);
      p.addOutput("v", DT::F64, swapped ? s : c);
      return p;
    };
    CHECK(make(false).fingerprint() != make(true).fingerprint());
  }

  TEST_CASE("lowering options are not in the program fingerprint") {
    LowerOptions bare;
    LowerOptions hoisting;
    hoisting.invariantInputs = {"x"};
    CHECK(bare.fingerprint() != hoisting.fingerprint());

    LowerOptions reordered;
    reordered.invariantInputs = {"y", "x"};
    LowerOptions sorted;
    sorted.invariantInputs = {"x", "y"};
    CHECK(reordered.fingerprint() == sorted.fingerprint());
  }
}

TEST_SUITE("Expr::Lower") {

  // The whole reason the pass exists: without reuse a 22-node chain needs 22
  // slots, which at f64 x 256 lanes is 44 KiB and no longer an L1 structure.
  TEST_CASE("liveness reuses slots along a chain") {
    Program p;
    p.addInput("x", DT::F64);
    Arena& a = p.arena();
    NodeId n = a.field("x");
    const NodeId one = a.konst(1.0);
    for (int i = 0; i < 20; ++i) {
      n = a.pw(Fn::Add, n, one);
    }
    p.addOutput("u", DT::F64, n);
    validate(p);

    const LoweredProgram lowered = lower(p);
    CHECK(lowered.run().code.size() == 22); // load x, const 1, 20 adds
    // x dies at the first add; the shared constant lives to the last one.
    CHECK(lowered.peakSlotCount() == 2);
    CHECK_FALSE(lowered.hasPrecompute());
  }

  TEST_CASE("a diamond is emitted once, not exponentially") {
    Program p;
    p.addInput("x", DT::F64);
    Arena& a = p.arena();
    NodeId n = a.field("x");
    for (int i = 0; i < 20; ++i) {
      n = a.pw(Fn::Mul, n, n);
    }
    p.addOutput("u", DT::F64, n);
    validate(p);

    const LoweredProgram lowered = lower(p);
    CHECK(lowered.run().code.size() == 21); // 1 load + 20 squarings, not 2^20
    CHECK(lowered.peakSlotCount() == 1);
  }

  TEST_CASE("nothing is hoisted by default") {
    Program p;
    p.addInput("x", DT::F64);
    p.addInput("t", DT::F64);
    Arena& a = p.arena();
    const GridId g = p.internGrid(reader::datafield::GridDesc{});
    const NodeId look = a.lookup(g, 0, {a.field("x")});
    p.addOutput("u", DT::F64, a.pw(Fn::Mul, look, a.pw(Fn::Sin, a.field("t"))));
    validate(p);

    const LoweredProgram lowered = lower(p);
    CHECK_FALSE(lowered.hasPrecompute());
    CHECK(lowered.persistentSlotCount() == 0);
    CHECK(countOpcode(lowered.run(), Opcode::Lookup) == 1);
  }

  // The analytic boundary condition: x is fixed for the run, t is not, so the
  // grid gather moves out of the timestep loop.
  TEST_CASE("an invariant grid lookup is hoisted out of the call") {
    Program p;
    p.addInput("x", DT::F64);
    p.addInput("t", DT::F64);
    Arena& a = p.arena();
    const GridId g = p.internGrid(reader::datafield::GridDesc{});
    const NodeId look = a.lookup(g, 0, {a.field("x")});
    p.addOutput("u", DT::F64, a.pw(Fn::Mul, look, a.pw(Fn::Sin, a.field("t"))));
    validate(p);

    LowerOptions options;
    options.invariantInputs = {"x"};
    options.invariantGrids = {g};
    const LoweredProgram lowered = lower(p, options);

    CHECK(lowered.hasPrecompute());
    CHECK(lowered.persistentSlotCount() == 1);
    CHECK(countOpcode(lowered.precompute(), Opcode::Lookup) == 1);
    CHECK(countOpcode(lowered.run(), Opcode::Lookup) == 0);
    CHECK(countOpcode(lowered.run(), Opcode::LoadPersistent) == 1);
    CHECK(lowered.precompute().persistent.size() == 1);
  }

  TEST_CASE("an updating grid is not hoisted even when its coordinates are") {
    Program p;
    p.addInput("x", DT::F64);
    Arena& a = p.arena();
    const GridId g = p.internGrid(reader::datafield::GridDesc{});
    p.addOutput("u", DT::F64, a.lookup(g, 0, {a.field("x")}));
    validate(p);

    LowerOptions options;
    options.invariantInputs = {"x"}; // but invariantGrids stays empty
    const LoweredProgram lowered = lower(p, options);
    CHECK_FALSE(lowered.hasPrecompute());
    CHECK(countOpcode(lowered.run(), Opcode::Lookup) == 1);
  }

  TEST_CASE("cheap invariants are recomputed rather than materialised") {
    Program p;
    p.addInput("x", DT::F64);
    p.addInput("t", DT::F64);
    Arena& a = p.arena();
    // x + x is invariant but costs one add; a persistent buffer plus a streaming
    // read is not a good trade for that.
    const NodeId cheap = a.pw(Fn::Add, a.field("x"), a.field("x"));
    p.addOutput("u", DT::F64, a.pw(Fn::Mul, cheap, a.field("t")));
    validate(p);

    LowerOptions options;
    options.invariantInputs = {"x"};
    const LoweredProgram lowered = lower(p, options);
    CHECK_FALSE(lowered.hasPrecompute());

    // Raising the price of the subtree flips the decision.
    Program q;
    q.addInput("x", DT::F64);
    q.addInput("t", DT::F64);
    Arena& b = q.arena();
    const NodeId pricey = b.pw(Fn::Exp, b.field("x"));
    q.addOutput("u", DT::F64, b.pw(Fn::Mul, pricey, b.field("t")));
    validate(q);
    const LoweredProgram lowered2 = lower(q, options);
    CHECK(lowered2.hasPrecompute());
    CHECK(lowered2.persistentSlotCount() == 1);
  }

  TEST_CASE("state reads and writes lower to persistent slots") {
    Program p;
    p.addInput("f", DT::F64);
    Arena& a = p.arena();
    const NodeId next = a.pw(Fn::Add, a.field("acc"), a.field("f"));
    p.addState("acc", 2.5, next);
    p.addOutput("u", DT::F64, a.field("acc"));
    validate(p);

    const LoweredProgram lowered = lower(p);
    CHECK(lowered.stateSlotCount() == 1);
    CHECK(lowered.persistentSlotCount() == 1);
    CHECK(countOpcode(lowered.run(), Opcode::LoadPersistent) == 1);
    REQUIRE(lowered.run().persistent.size() == 1);
    CHECK(lowered.run().persistent[0].target == 0);
    REQUIRE(lowered.run().outputs.size() == 1);

    // Parallel assignment: the output reads `acc` from before the update, so the
    // output store and the state store must not share a source slot with a value
    // that has already been overwritten.
    CHECK(lowered.run().outputs[0].source != lowered.run().persistent[0].source);
  }

  TEST_CASE("a value feeding two outputs is computed once") {
    Program p;
    p.addInput("x", DT::F64);
    Arena& a = p.arena();
    const NodeId shared = a.pw(Fn::Exp, a.field("x"));
    p.addOutput("u", DT::F64, a.pw(Fn::Add, shared, a.konst(1.0)));
    p.addOutput("v", DT::F64, a.pw(Fn::Mul, shared, a.konst(2.0)));
    validate(p);

    const LoweredProgram lowered = lower(p);
    CHECK(countOpcode(lowered.run(), Opcode::Pw) == 3); // exp, add, mul
    CHECK(lowered.run().outputs.size() == 2);
  }

  TEST_CASE("operand slots are live where they are read") {
    Program p;
    p.addInput("x", DT::F64);
    p.addInput("y", DT::F64);
    Arena& a = p.arena();
    const NodeId u = a.pw(Fn::Select,
                          a.pw(Fn::Lt, a.field("x"), a.field("y")),
                          a.pw(Fn::Sqrt, a.field("x")),
                          a.pw(Fn::Neg, a.field("y")));
    p.addOutput("u", DT::F64, u);
    validate(p);

    const LoweredProgram lowered = lower(p);
    // Every operand must name a slot below the stage's high-water mark, and the
    // ternary must carry exactly three of them.
    bool sawSelect = false;
    for (const Instruction& i : lowered.run().code) {
      for (std::int32_t k = 0; k < i.operandCount; ++k) {
        const std::int32_t slot = lowered.operands()[i.operandBegin + k];
        CHECK(slot >= 0);
        CHECK(slot < lowered.run().slotCount);
      }
      if (i.op == Opcode::Pw && i.fn == Fn::Select) {
        sawSelect = true;
        CHECK(i.operandCount == 3);
      }
    }
    CHECK(sawSelect);
  }
}

} // namespace seissol::expr::test
