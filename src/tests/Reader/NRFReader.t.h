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
    //NRF stores all values always in double, so use epsilon = 1e-16
    const double epsilon = std::numeric_limits<double>::epsilon();
    seissol::sourceterm::NRF nrf;
    seissol::sourceterm::readNRF("Testing/source_loh.nrf", nrf);

    TS_ASSERT_DELTA(nrf.centres[0](0),    0.0, epsilon);
    TS_ASSERT_DELTA(nrf.centres[0](1),    0.0, epsilon);
    TS_ASSERT_DELTA(nrf.centres[0](2), 2000.0, epsilon);

    TS_ASSERT_DELTA(nrf.subfaults[0].tan1(0), 0.0, epsilon); 
    TS_ASSERT_DELTA(nrf.subfaults[0].tan1(1), 1.0, epsilon); 
    TS_ASSERT_DELTA(nrf.subfaults[0].tan1(2), 0.0, epsilon); 

    TS_ASSERT_DELTA(nrf.subfaults[0].tan2(0), 0.0, epsilon); 
    TS_ASSERT_DELTA(nrf.subfaults[0].tan2(1), 0.0, epsilon); 
    TS_ASSERT_DELTA(nrf.subfaults[0].tan2(2), 1.0, epsilon); 

    TS_ASSERT_DELTA(nrf.subfaults[0].normal(0), 1.0, epsilon); 
    TS_ASSERT_DELTA(nrf.subfaults[0].normal(1), 0.0, epsilon); 
    TS_ASSERT_DELTA(nrf.subfaults[0].normal(2), 0.0, epsilon);

    TS_ASSERT_DELTA(nrf.subfaults[0].area, 3.0866008336686616479e+07, epsilon);
    TS_ASSERT_DELTA(nrf.subfaults[0].tinit, 0.0, epsilon);
    TS_ASSERT_DELTA(nrf.subfaults[0].timestep, 0.0002, epsilon);
    TS_ASSERT_DELTA(nrf.subfaults[0].mu, 0.0, epsilon);
    TS_ASSERT_EQUALS(nrf.source, 1);

    for (unsigned i = 0; i < slipRates.size(); i++) {
      //NRF File is in cm
      TS_ASSERT_DELTA(nrf.sliprates[0][i], slipRates[i]/100, epsilon);
      TS_ASSERT_DELTA(nrf.sliprates[1][i], 0, epsilon);
      TS_ASSERT_DELTA(nrf.sliprates[2][i], 0, epsilon);
    }


  }
};
