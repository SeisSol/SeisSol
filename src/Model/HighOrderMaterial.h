#pragma once

#include <Kernels/Common.h>
#include <Model/CommonDatastructures.h>
#include <array>
#include <string>
#include <type_traits>

namespace seissol::model {

template <typename BaseMaterialT, std::size_t Order>
class HighOrderMaterial : public Material {
  public:
  static_assert(std::is_base_of_v<Material, BaseMaterialT>,
                "BaseMaterialT must inherit from Material.");

  static constexpr std::size_t MaterialOrder = Order;
  static constexpr std::size_t Functions3D = kernels::getNumberOfBasisFunctions(Order);
  static constexpr std::size_t Functions2D = Order * (Order + 1) / 2;
  std::array<BaseMaterialT, Functions3D> materials;

  static constexpr std::size_t NumQuantities = BaseMaterialT::NumQuantities;
  static constexpr std::size_t NumElasticQuantities = BaseMaterialT::NumElasticQuantities;
  static constexpr std::size_t NumberPerMechanism = BaseMaterialT::NumberPerMechanism;
  static constexpr std::size_t Mechanisms = BaseMaterialT::Mechanisms;
  static constexpr MaterialType Type = BaseMaterialT::Type;
  static constexpr LocalSolver Solver = BaseMaterialT::Solver;
  static inline const std::string Text = BaseMaterialT::Text + "-h" + std::to_string(Order);
  static inline const std::array<std::string, NumQuantities> Quantities = BaseMaterialT::Quantities;

  using LocalSpecificData = typename BaseMaterialT::LocalSpecificData;
  using NeighborSpecificData = typename BaseMaterialT::NeighborSpecificData;

  ~HighOrderMaterial() override = default;

  HighOrderMaterial() = default;
  HighOrderMaterial(const std::vector<double>& data) {
    const auto perMaterial = data.size() / materials.size();
    for (std::size_t i = 0; i < materials.size(); ++i) {
      std::vector<double> perData(data.begin() + perMaterial * i,
                                  data.begin() + perMaterial * (i + 1));
      materials[i] = BaseMaterialT(perData);
    }
    rho = materials[0].rho;
  }

  // TODO: refine all of those
  [[nodiscard]] double getMaxWaveSpeed() const override { return materials[0].getMaxWaveSpeed(); }
  [[nodiscard]] double getPWaveSpeed() const override { return materials[0].getPWaveSpeed(); }
  [[nodiscard]] double getSWaveSpeed() const override { return materials[0].getSWaveSpeed(); }
  [[nodiscard]] double getMuBar() const override { return materials[0].getMuBar(); }
  [[nodiscard]] double getLambdaBar() const override { return materials[0].getLambdaBar(); }
  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override {
    // for now, take the base material here
    materials[0].getFullStiffnessTensor(fullTensor);
  }
  [[nodiscard]] MaterialType getMaterialType() const override { return Type; }
};

} // namespace seissol::model
