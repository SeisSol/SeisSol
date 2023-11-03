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
  int reflectionType;

  public:
  [[nodiscard]] double getITMTime() const;
  [[nodiscard]] double getITMVelocityScalingFactor() const;
  [[nodiscard]] double getITMStartingTime() const;
  [[nodiscard]] bool getITMToggle() const;
  [[nodiscard]] int getReflectionType() const;
  ITMParameters(double ITMTime, double ITMVelocityScalingFactor, double ITMStartingTime, bool ITMToggle, int reflectionType);
};

ITMParameters readITMParametersFromYaml(std::shared_ptr<YAML::Node>& params);
} // namespace seissol::initializers::ITM

#endif // SEISSOL_ITMCONFIG_H