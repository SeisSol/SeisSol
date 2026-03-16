// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SlipRatesData.h"
#include "SourceTerm/NRF.h"
#include "SourceTerm/NRFReader.h"
#include "TestHelper.h"

#include <cstdlib>

namespace seissol::unit_test {
TEST_CASE("NRF Reader" * doctest::test_suite("reader")) {
  seissol::sourceterm::NRF nrf;
  seissol::sourceterm::readNRF(tpath("Testing/source_loh.nrf").c_str(), nrf);

  CHECK(nrf.centres[0](0) == AbsApprox(0.0));
  CHECK(nrf.centres[0](1) == AbsApprox(0.0));
  CHECK(nrf.centres[0](2) == AbsApprox(2000.0));

  CHECK(nrf.subfaults[0].tan1(0) == AbsApprox(0.0));
  CHECK(nrf.subfaults[0].tan1(1) == AbsApprox(1.0));
  CHECK(nrf.subfaults[0].tan1(2) == AbsApprox(0.0));

  CHECK(nrf.subfaults[0].tan2(0) == AbsApprox(0.0));
  CHECK(nrf.subfaults[0].tan2(1) == AbsApprox(0.0));
  CHECK(nrf.subfaults[0].tan2(2) == AbsApprox(1.0));

  CHECK(nrf.subfaults[0].normal(0) == AbsApprox(1.0));
  CHECK(nrf.subfaults[0].normal(1) == AbsApprox(0.0));
  CHECK(nrf.subfaults[0].normal(2) == AbsApprox(0.0));

  CHECK(nrf.subfaults[0].area == AbsApprox(3.0866008336686616479e+07));
  CHECK(nrf.subfaults[0].tinit == AbsApprox(0.0));
  CHECK(nrf.subfaults[0].timestep == AbsApprox(0.0002));
  CHECK(nrf.subfaults[0].mu == AbsApprox(0.0));

  CHECK(nrf.size() == 1);
  CHECK(nrf.subfaults.size() == nrf.centres.size());
  CHECK(nrf.centres.size() + 1 == nrf.sroffsets.size());

  for (size_t dim = 0; dim < 3; dim++) {
    for (unsigned i = 0;
         i < std::min(static_cast<size_t>(nrf.sroffsets[1][dim]) - nrf.sroffsets[0][dim],
                      SlipRates[dim].size());
         i++) {
      // NRF File is in cm
      CHECK(nrf.sliprates[dim][i] == AbsApprox(SlipRates[dim][i] / 100));
    }
  }
}
} // namespace seissol::unit_test
