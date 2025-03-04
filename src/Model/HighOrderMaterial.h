#pragma once

#include <Kernels/Common.h>
#include <Model/CommonDatastructures.h>
#include <array>
#include <equation-elastic-6-double/init.h>
#include <limits>
#include <numeric>
#include <string>
#include <type_traits>

namespace seissol::model {

template <typename BaseMaterialT, std::size_t Order>
class HighOrderMaterial : public Material {
  public:
  static_assert(std::is_base_of_v<Material, BaseMaterialT>,
                "BaseMaterialT must inherit from Material.");

  static constexpr std::size_t MaterialOrder = Order;
  static constexpr std::size_t Samples3D = seissol::init::hompoints::Shape[0];

  // for now, interpolate
  std::array<BaseMaterialT, Samples3D> materials;

  static constexpr std::size_t NumQuantities = BaseMaterialT::NumQuantities;
  static constexpr std::size_t NumElasticQuantities = BaseMaterialT::NumElasticQuantities;
  static constexpr std::size_t NumberPerMechanism = BaseMaterialT::NumberPerMechanism;
  static constexpr std::size_t TractionQuantities = BaseMaterialT::TractionQuantities;
  static constexpr std::size_t Mechanisms = BaseMaterialT::Mechanisms;
  static constexpr MaterialType Type = BaseMaterialT::Type;
  static constexpr LocalSolver Solver = BaseMaterialT::Solver;
  static inline const std::string Text = BaseMaterialT::Text + "-h" + std::to_string(Order);
  static inline const std::array<std::string, NumQuantities> Quantities = BaseMaterialT::Quantities;
  static constexpr std::size_t Parameters = Samples3D * BaseMaterialT::Parameters;

  static constexpr bool SupportsDR = BaseMaterialT::SupportsDR;
  static constexpr bool SupportsLTS = BaseMaterialT::SupportsLTS;

  using LocalSpecificData = typename BaseMaterialT::LocalSpecificData;
  using NeighborSpecificData = typename BaseMaterialT::NeighborSpecificData;

  ~HighOrderMaterial() override = default;

  HighOrderMaterial() = default;
  HighOrderMaterial(const std::vector<double>& data) {
    const auto perMaterial = data.size() / materials.size();
    assert(perMaterial == BaseMaterialT::Parameters);
    for (std::size_t i = 0; i < materials.size(); ++i) {
      std::vector<double> perData(data.begin() + perMaterial * i,
                                  data.begin() + perMaterial * (i + 1));
      materials[i] = BaseMaterialT(perData);
    }
    rho = materials[0].rho;
  }

  // TODO: refine all of those
  [[nodiscard]] double getMaxWaveSpeed() const override {
    // TODO: search maximum of polynomial
    double maxWavespeed = 0;
    for (const auto& material : materials) {
      maxWavespeed = std::max(maxWavespeed, material.getMaxWaveSpeed());
    }
    return maxWavespeed;
  }
  [[nodiscard]] double getPWaveSpeed() const override {
    double maxWavespeed = 0;
    for (const auto& material : materials) {
      maxWavespeed = std::max(maxWavespeed, material.getPWaveSpeed());
    }
    return maxWavespeed;
  }
  [[nodiscard]] double getSWaveSpeed() const override {
    double maxWavespeed = 0;
    for (const auto& material : materials) {
      maxWavespeed = std::max(maxWavespeed, material.getSWaveSpeed());
    }
    return maxWavespeed;
  }
  [[nodiscard]] double getMuBar() const override {
    double mubar = 0;
    for (const auto& material : materials) {
      mubar += material.getMuBar();
    }
    return mubar / materials.size();
  }
  [[nodiscard]] double getLambdaBar() const override {
    double lambdabar = 0;
    for (const auto& material : materials) {
      lambdabar += material.getLambdaBar();
    }
    return lambdabar / materials.size();
  }
  void getFullStiffnessTensor(std::array<double, 81>& fullTensor) const override {
    // for now, take the base material here
    materials[0].getFullStiffnessTensor(fullTensor);
  }
  [[nodiscard]] MaterialType getMaterialType() const override { return Type; }
  [[nodiscard]] double maximumTimestep() const override {
    double timestepbound = std::numeric_limits<double>::infinity();
    for (const auto& material : materials) {
      timestepbound = std::min(timestepbound, material.getMaxWaveSpeed());
    }
    return timestepbound;
  }
};

} // namespace seissol::model
