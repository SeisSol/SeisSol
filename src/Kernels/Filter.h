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
  ExponentialFilter() : conf() {};
  ExponentialFilter(initializer::parameters::FilterParameters conf, unsigned int dimensions)
      : conf(conf), dimensions(dimensions) {
    // TODO Better error handling
    assert((conf.order % 2) == 0);
    assert((dimensions == 2) || (dimensions == 3));
    if (dimensions == 2) {
      for (unsigned int ord = 0; ord <= CONVERGENCE_ORDER; ++ord) {
        for (unsigned int j = 0; j <= ord; ++j) {
          filterMatrix.emplace_back(evaluateFilter(ord - j, j, 0));
        }
      }
    } else {
        for (unsigned int ord = 0; ord < CONVERGENCE_ORDER; ord++) {
          for (unsigned int k = 0; k <= ord; k++) {
            for (unsigned int j = 0; j <= ord - k; j++) {
              filterMatrix.emplace_back(evaluateFilter(ord - j - k, j, k));
            }
          }
        }
      }
  }

  real inline getFilterCoeff(unsigned idx) const {
    // TODO Depending on the implementation, it might be faster to resize filterMatrix to
    // include an alignment
    if (conf.type != initializer::parameters::FilterTypes::Exponential) return 1;

    if(idx < filterMatrix.size()) {
      return filterMatrix[idx];
    } else {
      return 1.0; // Don't touch other modes
    }
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

  std::vector<real> filterMatrix{};
  initializer::parameters::FilterParameters conf{};
  unsigned int dimensions = 3;
};

void inline applyFilter(NeighborData& data, const ExponentialFilter& filter) {
  auto qView = init::Q::view::create(data.dofs);
  for (unsigned dofIdx = 0; dofIdx < tensor::Q::Shape[1]; ++dofIdx) {
    for (unsigned basisIdx = 0; basisIdx < tensor::Q::Shape[0]; ++basisIdx) {
      qView(basisIdx, dofIdx) *= filter.getFilterCoeff(basisIdx);
    }
  }
}

auto inline computeDRFilterMatrix(const ExponentialFilter& filter) {
  auto filterMatrixData = std::array<real, tensor::filter::Size>{0.0};
  auto filterWeightsData = std::array<real, tensor::filter::Size>{0.0};
  auto filterWeights = init::filterWeights::view::create(filterWeightsData.data());
  for (unsigned i = 0; i < filterWeights.shape(0); ++i) {
    filterWeights(i,i)  = filter.getFilterCoeff(i);
  }

  auto compFilterKrnl = dynamicRupture::kernel::computeFilterMatrix{};
  compFilterKrnl.filter = filterMatrixData.data();
  compFilterKrnl.filterWeights = filterWeightsData.data();
  compFilterKrnl.V2QuadTo2m = init::V2QuadTo2m::Values;
  compFilterKrnl.V2mTo2Quad = init::V2mTo2Quad::Values;
  compFilterKrnl.execute();

  return filterMatrixData;
}

} // namespace seissol::kernels
#endif // SEISSOL_FILTER_H
