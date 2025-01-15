// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "tests/TestHelper.h"
#include <cstdlib>

#include "SourceTerm/NRF.h"
#include "SourceTerm/NRFReader.h"

#include "SlipRatesData.h"

namespace seissol::unit_test {
TEST_CASE("NRF Reader") {
  seissol::sourceterm::NRF nrf;
  seissol::sourceterm::readNRF("Testing/source_loh.nrf", nrf);

  REQUIRE(nrf.centres[0](0) == AbsApprox(0.0));
  REQUIRE(nrf.centres[0](1) == AbsApprox(0.0));
  REQUIRE(nrf.centres[0](2) == AbsApprox(2000.0));

  REQUIRE(nrf.subfaults[0].tan1(0) == AbsApprox(0.0));
  REQUIRE(nrf.subfaults[0].tan1(1) == AbsApprox(1.0));
  REQUIRE(nrf.subfaults[0].tan1(2) == AbsApprox(0.0));

  REQUIRE(nrf.subfaults[0].tan2(0) == AbsApprox(0.0));
  REQUIRE(nrf.subfaults[0].tan2(1) == AbsApprox(0.0));
  REQUIRE(nrf.subfaults[0].tan2(2) == AbsApprox(1.0));

  REQUIRE(nrf.subfaults[0].normal(0) == AbsApprox(1.0));
  REQUIRE(nrf.subfaults[0].normal(1) == AbsApprox(0.0));
  REQUIRE(nrf.subfaults[0].normal(2) == AbsApprox(0.0));

  REQUIRE(nrf.subfaults[0].area == AbsApprox(3.0866008336686616479e+07));
  REQUIRE(nrf.subfaults[0].tinit == AbsApprox(0.0));
  REQUIRE(nrf.subfaults[0].timestep == AbsApprox(0.0002));
  REQUIRE(nrf.subfaults[0].mu == AbsApprox(0.0));

  REQUIRE(nrf.size() == 1);
  REQUIRE(nrf.subfaults.size() == nrf.centres.size());
  REQUIRE(nrf.centres.size() + 1 == nrf.sroffsets.size());

  for (size_t dim = 0; dim < 3; dim++) {
    for (unsigned i = 0;
         i < std::min((size_t)nrf.sroffsets[1][dim] - nrf.sroffsets[0][dim], SlipRates[dim].size());
         i++) {
      // NRF File is in cm
      REQUIRE(nrf.sliprates[dim][i] == AbsApprox(SlipRates[dim][i] / 100));
    }
  }
}
} // namespace seissol::unit_test
