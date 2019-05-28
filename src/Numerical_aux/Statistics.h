/**
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 **/

#ifndef NUMERICAL_AUX_STATISTICS_H_
#define NUMERICAL_AUX_STATISTICS_H_

#include <vector>

namespace seissol {
  namespace statistics {
    struct Summary {
      Summary(double value = 0.0);
      Summary(std::vector<double> const& values);

      double mean;
      double std;
      double min;
      double median;
      double max;
    };

    auto parallelSummary(double value) -> Summary;
  }
}

#endif
