#ifndef SEISSOL_FILTER_H
#define SEISSOL_FILTER_H

#include <vector>
#include <array>

#include "Initializer/InputParameters.hpp"
#include "Kernels/common.hpp"
#include "Kernels/Interface.hpp"
#include "Kernels/NeighborBase.h"
namespace seissol::kernels {

class Filter {
  public:
  virtual real getFilterCoeff(unsigned idx) const = 0;

  protected:
  Filter() = default;
  Filter(initializer::parameters::FilterParameters conf, unsigned int dimensions);

  initializer::parameters::FilterParameters conf{};
  unsigned int dimensions = 3;
};

class IdentityFilter : public Filter {
  public:
  IdentityFilter(initializer::parameters::FilterParameters conf, unsigned int dimensions);
  real getFilterCoeff(unsigned idx) const;
};

class ExponentialFilter : public Filter {
  public:
  ExponentialFilter();
  ExponentialFilter(initializer::parameters::FilterParameters conf, unsigned int dimensions);

  real getFilterCoeff(unsigned idx) const;

  private:
  real evaluateFilter(unsigned i, unsigned j, unsigned k) const;
  std::vector<real> filterMatrix{};
};

std::array<real, tensor::drFilter::Size> computeDRFilterMatrix(const Filter& filter);

std::array<real, tensor::volumeFilter::Size> computeFilterMatrix(const Filter& filter);

std::unique_ptr<Filter> makeFilter(const initializer::parameters::FilterParameters& conf,
                                   unsigned int dimensions);

} // namespace seissol::kernels
#endif // SEISSOL_FILTER_H
