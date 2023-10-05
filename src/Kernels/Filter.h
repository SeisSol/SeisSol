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
      /*
      // TODO(Lukas) Remove
  struct Configuration {
    real alpha = -std::log(std::numeric_limits<real>::epsilon());
    unsigned int filterOrder = 32;
    unsigned int filterCutoff = 0;
  };
       */

  ExponentialFilter() : conf() {};
  ExponentialFilter(initializer::parameters::FilterParameters conf) : conf(conf) {
    if (conf.type != initializer::parameters::FilterTypes::Exponential) {
      std::fill(filterMatrix.begin(), filterMatrix.end(), 1.0);
      return;
    }

    assert((conf.order % 2) == 0);
    unsigned int i = 0;
    for (unsigned int ord = 0; ord < CONVERGENCE_ORDER; ord++) {
      for (unsigned int k = 0; k <= ord; k++) {
        for (unsigned int j = 0; j <= ord - k; j++) {
          assert(i < filterMatrix.size());
          filterMatrix[i++] = evaluateFilter(ord - j - k, j, k);
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
    // Hesthaven Nodal DG p. 129, eq. 5.16
    const auto eta = (1.0 * n) / N;
    const auto etaCutoff = (1.0 * conf.cutoff) / N;

    if (0 <= eta && eta <= etaCutoff) {
      return 1.0;
    }

    return std::exp(-conf.alpha *
                    std::pow(((eta - etaCutoff) / (1 - etaCutoff)), conf.order));
  }

  std::array<real, getNumberOfAlignedBasisFunctions()> filterMatrix{ 0.0 };
  initializer::parameters::FilterParameters conf{};
};

void inline applyFilter(NeighborData& data, const ExponentialFilter& filter) {
  auto qView = init::Q::view::create(data.dofs);
  for (unsigned basisIdx = 0; basisIdx < tensor::Q::Shape[0]; ++basisIdx) {
    for (unsigned dofIdx = 0; dofIdx < tensor::Q::Shape[1]; ++dofIdx) {
      qView(basisIdx, dofIdx) *= filter.getFilterCoeff(basisIdx);
    }
  }
}

} // namespace seissol::kernels
#endif // SEISSOL_FILTER_H
