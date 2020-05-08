#include <cxxtest/TestSuite.h>

#include <SourceTerm/NRFReader.h>
#include <SourceTerm/NRF.h>

#include "slipRatesData.h"

namespace seissol {
  namespace unit_test {
    class NRFReaderTestSuite;
  }
}

class seissol::unit_test::NRFReaderTestSuite : public CxxTest::TestSuite
{
public:
  void testNRFReader()
  {
    //NRF stores all values always in double, so use eps = 1e-16
    constexpr double eps = 1e-16;
    seissol::sourceterm::NRF nrf;
    seissol::sourceterm::readNRF("Testing/source_loh.nrf", nrf);

    TS_ASSERT_DELTA(nrf.centres[0].x,    0.0, eps);
    TS_ASSERT_DELTA(nrf.centres[0].y,    0.0, eps);
    TS_ASSERT_DELTA(nrf.centres[0].z, 2000.0, eps);

    TS_ASSERT_DELTA(nrf.subfaults[0].tan1.x, 0.0, eps); 
    TS_ASSERT_DELTA(nrf.subfaults[0].tan1.y, 1.0, eps); 
    TS_ASSERT_DELTA(nrf.subfaults[0].tan1.z, 0.0, eps); 

    TS_ASSERT_DELTA(nrf.subfaults[0].tan2.x, 0.0, eps); 
    TS_ASSERT_DELTA(nrf.subfaults[0].tan2.y, 0.0, eps); 
    TS_ASSERT_DELTA(nrf.subfaults[0].tan2.z, 1.0, eps); 

    TS_ASSERT_DELTA(nrf.subfaults[0].normal.x, 1.0, eps); 
    TS_ASSERT_DELTA(nrf.subfaults[0].normal.y, 0.0, eps); 
    TS_ASSERT_DELTA(nrf.subfaults[0].normal.z, 0.0, eps);

    TS_ASSERT_DELTA(nrf.subfaults[0].area, 3.0866008336686616479e+07, eps);
    TS_ASSERT_DELTA(nrf.subfaults[0].tinit, 0.0, eps);
    TS_ASSERT_DELTA(nrf.subfaults[0].timestep, 0.0002, eps);
    TS_ASSERT_DELTA(nrf.subfaults[0].mu, 0.0, eps);
    TS_ASSERT_EQUALS(nrf.source, 1);

    for (int i = 0; i < 80; i++) {
      //NRF File is in cm
      TS_ASSERT_DELTA(nrf.sliprates[0][i], slipRates[i]/100, eps);
      TS_ASSERT_DELTA(nrf.sliprates[1][i], 0, eps);
      TS_ASSERT_DELTA(nrf.sliprates[2][i], 0, eps);
    }


  }
};
