#ifndef SEISSOL_ITMCONFIG_H
#define SEISSOL_ITMCONFIG_H

#include <memory>
#include <yaml-cpp/yaml.h>

namespace seissol::initializers::ITM {
class ITMParameters {
  private:
  double ITMTime;
  double ITMVelocityScalingFactor;
  double ITMStartingTime;
  bool ITMToggle;

  public:
  [[nodiscard]] double getITMTime() const;
  [[nodiscard]] double getITMVelocityScalingFactor() const;
  [[nodiscard]] double getITMStartingTime() const;
  [[nodiscard]] bool getITMToggle() const;
  ITMParameters(double ITMTime, double ITMVelocityScalingFactor, double ITMStartingTime, bool ITMToggle);
};

ITMParameters readITMParametersFromYaml(std::shared_ptr<YAML::Node>& params);
} // namespace seissol::initializers::ITM

#endif // SEISSOL_ITMCONFIG_H