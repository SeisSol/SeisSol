#ifndef SEISSOL_ITMCONFIG_H
#define SEISSOL_ITMCONFIG_H

#include <memory>
#include <yaml-cpp/yaml.h>

namespace seissol::initializers::ITM {
class ITMParameters {
  private:
  double ITMTime;
  public:
  [[nodiscard]] double getITMTime() const;
  ITMParameters(double ITMTime);
};

ITMParameters readITMParametersFromYaml(std::shared_ptr<YAML::Node>& params);
} // namespace seissol::initializers::ITM

#endif //SEISSOL_ITMCONFIG_H