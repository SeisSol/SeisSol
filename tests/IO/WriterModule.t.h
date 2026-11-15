// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Writer/Module/WriterModule.h"

namespace seissol::unit_test {

TEST_CASE("IO/WriterModule: a restart resumes the output numbering" * doctest::test_suite("io")) {
  using seissol::io::writer::module::outputCountBefore;

  // a run that starts over writes output 0 first
  CHECK(outputCountBefore(0.0, 0.5) == 0);

  // one resuming from a checkpoint continues behind what the run it continues already wrote,
  // the output at the checkpoint time included
  CHECK(outputCountBefore(0.5, 0.5) == 2);
  CHECK(outputCountBefore(2.0, 0.5) == 5);

  // a checkpoint that does not sit on an output time belongs behind the last one before it
  CHECK(outputCountBefore(2.4, 0.5) == 5);
  CHECK(outputCountBefore(2.6, 0.5) == 6);

  // a checkpoint time that should sit on an output may arrive a rounding error below it, and
  // still belongs behind that output rather than in front of it
  CHECK(outputCountBefore(2.0 - 1e-13, 0.5) == 5);

  // an interval that never fires leaves the numbering alone
  CHECK(outputCountBefore(2.0, 0.0) == 0);
}

} // namespace seissol::unit_test
