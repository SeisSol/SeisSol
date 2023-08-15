#ifndef COMMON_CONFIGTENSOR_HPP
#define COMMON_CONFIGTENSOR_HPP

#include "generated_code/kernel.h"
#include "generated_code/init.h"
#include "generated_code/tensor.h"

namespace seissol {

namespace model {
class ElasticMaterial;
template <std::size_t>
class ViscoElasticMaterial;
class PoroElasticMaterial;
class AnisotropicMaterial;
} // namespace model

template <typename MaterialT>
struct MaterialToYATeTo {};

template <>
struct MaterialToYATeTo<seissol::model::ElasticMaterial> {
  constexpr static std::size_t Position = 0;
};

template <>
struct MaterialToYATeTo<seissol::model::ViscoElasticMaterial<3>> {
  constexpr static std::size_t Position = 1;
};

template <>
struct MaterialToYATeTo<seissol::model::PoroElasticMaterial> {
  constexpr static std::size_t Position = 2;
};

template <>
struct MaterialToYATeTo<seissol::model::AnisotropicMaterial> {
  constexpr static std::size_t Position = 3;
};

template <typename Config>
struct Yateto {
  using Kernel = seissol::
      kernel<MaterialToYATeTo<typename Config::MaterialT>::Position, Config::ConvergenceOrder, 0>;
  using Init = seissol::
      init<MaterialToYATeTo<typename Config::MaterialT>::Position, Config::ConvergenceOrder, 0>;
  using Tensor = seissol::
      tensor<MaterialToYATeTo<typename Config::MaterialT>::Position, Config::ConvergenceOrder, 0>;
};

} // namespace seissol

#endif
