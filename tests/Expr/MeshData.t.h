// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

// Data on the mesh (datasets of the kind cells): a value per cell of the mesh, read as `d.c`
// without coordinates -- an input channel of the program, after all others, which a consumer that
// knows the cell of every point binds; in the core format and in the full language alike.

#include <doctest.h>

#include "Expr/Backend.h"
#include "Expr/Binding.h"
#include "Expr/Core.h"
#include "Expr/Program.h"
#include "Expr/SderivFrontend.h"
#include "Reader/Datafield/Grid.h"
#include "Reader/Scripting/DataTable.h"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::unit_test::meshdata {

using namespace seissol::expr;

/// A core file that reads the data per cell `d` of the file material.h5.
constexpr const char* Core = R"(#!sderiv 1
in x
dataset d (kind = "cells", path = "material.h5", variable = "material",
           components = ("rho", "vs"))
def t0 = d.vs * d.vs
def t1 = d.rho * t0
def t2 = t1 + x
out mu = t2
)";

/// The core text of `program` without its `meta line`s, the positions the full language gives.
std::string withoutPositions(const Program& program) {
  std::string out;
  std::string text = printCore(program);
  std::size_t from = 0;
  while (from < text.size()) {
    const std::size_t end = text.find('\n', from);
    const std::string line = text.substr(from, end - from + 1);
    if (line.rfind("meta line", 0) != 0) {
      out += line;
    }
    from = end == std::string::npos ? text.size() : end + 1;
  }
  return out;
}

std::vector<std::string> inputNames(const Program& program) {
  std::vector<std::string> names;
  for (const auto& input : program.inputs()) {
    names.push_back(input.name);
  }
  return names;
}

/// The message of the error `compile` throws, empty where it throws none.
template <typename F>
std::string failureOf(F&& compile) {
  try {
    compile();
  } catch (const std::exception& error) {
    return error.what();
  }
  return "";
}

TEST_SUITE("expr data on the mesh") {
  TEST_CASE("a core file reads data per cell as input channels, after the others") {
    const Program program = loadCore(Core, "/data/run/material.sderiv");
    // the components after `in`, in the order of the declaration
    CHECK(inputNames(program) == std::vector<std::string>{"x", "d.rho", "d.vs"});
    REQUIRE(program.grids().size() == 1);
    const auto& desc = program.grids()[0];
    CHECK(desc.kind == reader::datafield::GridKind::Cells);
    CHECK(!reader::datafield::sampledAtCoordinates(desc.kind));
    CHECK(desc.name == "d");
    CHECK(desc.path == "/data/run/material.h5");
    CHECK(desc.variable == "material");
    CHECK(desc.componentNames == std::vector<std::string>{"rho", "vs"});
    // the normal form is what the file is, and reads back into the same program
    const Program local = loadCore(Core, "");
    CHECK(printCore(local) == Core);
    CHECK(loadCore(printCore(local), "").fingerprint() == local.fingerprint());
  }

  TEST_CASE("the full language reads data per cell as the core format does") {
    const Program program = compileSderivModule(
        "dataset d (kind = \"cells\", path = \"material.h5\", variable = \"material\",\n"
        "           components = (\"rho\", \"vs\"))\n"
        "out def mu = d.rho * (d.vs * d.vs) + x\n");
    CHECK(inputNames(program) == std::vector<std::string>{"x", "d.rho", "d.vs"});
    CHECK(withoutPositions(program) == Core);
    CHECK(program.fingerprint() == loadCore(Core, "").fingerprint());
    // a component it does not read is an input all the same, as the declaration gives it
    const Program some = compileSderivModule(
        "dataset d (kind = \"cells\", path = \"m.h5\", variable = \"v\", components = (\"a\", "
        "\"b\"))\nout def u = d.b\n");
    CHECK(inputNames(some) == std::vector<std::string>{"d.a", "d.b"});
  }

  TEST_CASE("two names of one file are two datasets") {
    const Program program = loadCore("#!sderiv 1\n"
                                     "dataset a (kind = \"cells\", path = \"m.h5\", variable = "
                                     "\"v\", components = (\"rho\"))\n"
                                     "dataset b (kind = \"cells\", path = \"m.h5\", variable = "
                                     "\"v\", components = (\"rho\"))\n"
                                     "def t0 = a.rho + b.rho\n"
                                     "out u = t0\n",
                                     "");
    CHECK(program.grids().size() == 2);
    CHECK(inputNames(program) == std::vector<std::string>{"a.rho", "b.rho"});
    CHECK(loadCore(printCore(program), "").fingerprint() == program.fingerprint());
  }

  TEST_CASE("data per cell has no coordinates, no interpolation and no boundary") {
    const std::string declaration =
        "dataset d (kind = \"cells\", path = \"m.h5\", variable = \"v\", components = (\"rho\"))\n";
    const auto core = [&](const std::string& rest) {
      return failureOf(
          [&] { static_cast<void>(loadCore("#!sderiv 1\nin x\n" + declaration + rest, "")); });
    };
    const auto full = [&](const std::string& text) {
      return failureOf([&] { static_cast<void>(compileSderivModule(text)); });
    };
    CHECK(core("out u = d.rho(x)\n").find("read without coordinates, as `d.rho`") !=
          std::string::npos);
    CHECK(core("out u = d\n").find("holds a value per cell of the mesh, read with a component") !=
          std::string::npos);
    CHECK(core("out u = d.vp\n").find("has no component `vp`") != std::string::npos);
    CHECK(full(declaration + "out def u = d.rho(x)\n").find("read without coordinates") !=
          std::string::npos);
    CHECK(full(declaration + "out def u = d\n").find("a value per cell of the mesh") !=
          std::string::npos);
    CHECK(full(declaration + "out def u = d(x)\n").find("read without coordinates") !=
          std::string::npos);
    CHECK(full(declaration + "out def u = d.vp\n").find("has no component `vp`") !=
          std::string::npos);
    // its attributes: those of its kind, without the interpolation and the boundary of a grid
    CHECK(core("out u = d.rho\n").empty());
    CHECK(failureOf([] {
            static_cast<void>(loadCore("#!sderiv 1\ndataset d (kind = \"cells\", path = \"m.h5\", "
                                       "variable = \"v\", components = (\"rho\"), interpolation = "
                                       "\"linear\")\nout u = d.rho\n",
                                       ""));
          }).find("a dataset of the kind \"cells\" has no attribute `interpolation`") !=
          std::string::npos);
    CHECK(failureOf([] {
            static_cast<void>(loadCore("#!sderiv 1\ndataset d (kind = \"cells\", path = \"m.h5\", "
                                       "components = (\"rho\"))\nout u = d.rho\n",
                                       ""));
          }).find("gives no variable") != std::string::npos);
    CHECK(full("dataset d (kind = \"cells\", path = \"m.h5\", components = (\"rho\"))\n"
               "out def u = d.rho\n")
              .find("gives no variable") != std::string::npos);
    // and a grid is sampled at coordinates, as before
    CHECK(failureOf([] {
            static_cast<void>(loadCore(
                "#!sderiv 1\ndataset m (kind = \"inline\", components = (\"v\"), knots = ((0, "
                "1)), values = ((0), (1)), interpolation = \"linear\", boundary = \"clamp\")\n"
                "out u = m.v\n",
                ""));
          }).find("is sampled at coordinates") != std::string::npos);
    // a channel of data on the mesh is no other input
    CHECK(failureOf([] {
            static_cast<void>(loadCore("#!sderiv 1\nin `d.rho`\ndataset d (kind = \"cells\", path "
                                       "= \"m.h5\", variable = \"v\", components = (\"rho\"))\n"
                                       "out u = d.rho\n",
                                       ""));
          }).find("has the name of a component of the data on the mesh") != std::string::npos);
  }

  TEST_CASE("a consumer binds data per cell as inputs, and one that does not names it") {
    const Program program = loadCore(Core, "");
    std::vector<double> x{1.0, 2.0};
    std::vector<double> rho{2000.0, 3000.0};
    std::vector<double> vs{10.0, 20.0};
    std::vector<double> mu(2, 0.0);
    reader::scripting::DataTable table(2);
    table.bindViewConst<double>("x", reader::scripting::Direction::In, x.data());
    table.bindViewConst<double>("d.rho", reader::scripting::Direction::In, rho.data());
    table.bindViewConst<double>("d.vs", reader::scripting::Direction::In, vs.data());
    table.bindView<double>("mu", reader::scripting::Direction::Out, mu.data());
    auto binding = Binding::bind(program, table);
    // the store has no grid to open for it
    reader::datafield::GridStore store;
    CHECK(internGrids(program, store) == std::vector<std::size_t>{NoStoreGrid});
    CHECK(store.size() == 0);
    BackendOptions options;
    options.quiet = true;
    const auto kernel = makeKernel(program, binding, store, options);
    kernel->precompute(table);
    kernel->run(table);
    CHECK(mu[0] == 2000.0 * 100.0 + 1.0);
    CHECK(mu[1] == 3000.0 * 400.0 + 2.0);

    reader::scripting::DataTable without(2);
    without.bindViewConst<double>("x", reader::scripting::Direction::In, x.data());
    without.bindView<double>("mu", reader::scripting::Direction::Out, mu.data());
    CHECK(failureOf([&] { static_cast<void>(Binding::bind(program, without)); })
              .find("a value per cell of the mesh (the dataset `d` of material.h5), which this "
                    "consumer does not bind") != std::string::npos);
  }
}

} // namespace seissol::unit_test::meshdata
