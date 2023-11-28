#include "Filter.h"

namespace seissol::kernels {
using namespace initializer::parameters;

Filter::Filter(FilterParameters conf, unsigned int dimensions) {}

IdentityFilter::IdentityFilter(FilterParameters conf, unsigned int dimensions)
    : Filter(conf, dimensions){};

real IdentityFilter::getFilterCoeff(unsigned) const { return 1.0; }
ExponentialFilter::ExponentialFilter(FilterParameters conf, unsigned int dimensions)
    : Filter(conf, dimensions) {
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
ExponentialFilter::ExponentialFilter() : Filter() {}

real ExponentialFilter::getFilterCoeff(unsigned int idx) const {
  if (idx < filterMatrix.size()) {
    return filterMatrix[idx];
  } else {
    return 1.0; // Don't touch other modes
  }
}

real ExponentialFilter::evaluateFilter(unsigned int i, unsigned int j, unsigned int k) const {
  const auto n = i + j + k;
  const auto N = CONVERGENCE_ORDER;
  // Hesthaven Nodal DG p. 129, eq. 5.16
  const auto eta = (1.0 * n) / N;
  const auto etaCutoff = (1.0 * conf.cutoff) / N;

  if (0 <= eta && eta <= etaCutoff) {
    return 1.0;
  }

  return std::exp(-conf.alpha * std::pow(((eta - etaCutoff) / (1 - etaCutoff)), conf.order));
}

std::array<real, tensor::drFilter::Size> computeDRFilterMatrix(const Filter& filter) {
  auto filterMatrixData = std::array<real, tensor::drFilter::Size>{0.0};
  auto filterWeightsData = std::array<real, tensor::drFilter::Size>{0.0};
  auto filterWeights = init::filterWeights::view::create(filterWeightsData.data());
  for (unsigned i = 0; i < filterWeights.shape(0); ++i) {
    filterWeights(i, i) = filter.getFilterCoeff(i);
  }

  auto compFilterKrnl = dynamicRupture::kernel::computeFilterMatrix{};
  compFilterKrnl.drFilter = filterMatrixData.data();
  compFilterKrnl.filterWeights = filterWeightsData.data();
  compFilterKrnl.V2QuadTo2m = init::V2QuadTo2m::Values;
  compFilterKrnl.V2mTo2Quad = init::V2mTo2Quad::Values;
  compFilterKrnl.execute();

  return filterMatrixData;
}

std::array<real, tensor::volumeFilter::Size> computeFilterMatrix(const Filter& filter) {
  auto filterMatrixData = std::array<real, tensor::volumeFilter::Size>{0.0};
  auto filterMatrix = init::volumeFilter::view::create(filterMatrixData.data());
  for (unsigned i = 0; i < filterMatrix.shape(0); ++i) {
    filterMatrix(i, i) = filter.getFilterCoeff(i);
  }
  return filterMatrixData;
}

std::unique_ptr<Filter> makeFilter(const initializer::parameters::FilterParameters& conf,
                                   unsigned int dimensions) {
  using namespace initializer::parameters;
  switch (conf.type) {
  case FilterTypes::Identity:
    return std::make_unique<IdentityFilter>(conf, dimensions);
  case FilterTypes::Exponential:
    logInfo() << "Creating exponential filter with order " << conf.order;
    return std::make_unique<ExponentialFilter>(conf, dimensions);
  default:
    logError() << "FilterType" << static_cast<int>(conf.type) << "is not supported by makeFilter()";
    return nullptr;
  }
}


} // namespace seissol::kernels
