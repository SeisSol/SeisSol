#ifndef SEISSOL_FILTER_H
#define SEISSOL_FILTER_H

#include <vector>
#include <array>

#include "Kernels/common.hpp"
#include "Kernels/Interface.hpp"
#include "Kernels/NeighborBase.h"
namespace seissol::kernels {

class ExponentialFilter {
  public:
  struct Configuration {
    real alpha = -std::log(std::numeric_limits<real>::epsilon());
    unsigned int filterOrder = 32;
    unsigned int filterCutoff = 0;
  };

  ExponentialFilter() : conf() {};
  ExponentialFilter(Configuration conf) : conf(conf) {
    assert((conf.filterOrder % 2) == 0);
    unsigned int i = 0;
    for (unsigned int ord = 0; ord < CONVERGENCE_ORDER; ord++) {
      for (unsigned int k = 0; k <= ord; k++) {
        for (unsigned int j = 0; j <= ord - k; j++) {
          assert(i < filterMatrix.size());
          filterMatrix[i++] = evaluateFilter(ord - j - k, j, k);
          //std::cout << ord << "\t" <<ord - j - k << "\t" << j << "\t" << k << "\t" << filterMatrix[i-1] << "\n"<< std::endl;
        }
      }
    }
  }

  real inline getFilterCoeff(unsigned idx) const {
    assert(idx < filterMatrix.size());
    return filterMatrix[idx];
  }

  private:
  real inline evaluateFilter(unsigned i, unsigned j, unsigned k) const {
    const auto n = i + j + k;
    const auto N = CONVERGENCE_ORDER;
    // Hesthaven p. 129
    // Hesthaven eq. 5.16
    const auto eta = (1.0 * n) / N;
    //std::cout << n << "\t" << eta << std::endl;
    const auto etaCutoff = (1.0 * conf.filterCutoff) / N;

    if (0 <= eta && eta <= etaCutoff) {
      return 1.0;
    }

    return std::exp(-conf.alpha *
                    std::pow(((eta - etaCutoff) / (1 - etaCutoff)), conf.filterOrder));
  }

  std::array<real, getNumberOfBasisFunctions()> filterMatrix{};
  Configuration conf{};
};

void inline applyFilter(NeighborData& data, const ExponentialFilter& filter) {
  // TODO(Lukas) Create yateto kernel
  auto qView = init::Q::view::create(data.dofs);
  for (unsigned basisIdx = 0; basisIdx < tensor::Q::Shape[0]; ++basisIdx) {
    for (unsigned dofIdx = 0; dofIdx < tensor::Q::Shape[1]; ++dofIdx) {
      qView(basisIdx, dofIdx) *= filter.getFilterCoeff(basisIdx);
    }
  }
}

} // namespace seissol::kernels
#endif // SEISSOL_FILTER_H
