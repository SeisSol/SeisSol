
#ifndef INPUT_PARAMETERS_HPP_
#define INPUT_PARAMETERS_HPP_

#include <string>
#include <yaml-cpp/yaml.h>

#include <xdmfwriter/XdmfWriter.h>


namespace seissol::initializer::parameters {

enum class OutputRefinement : int { NoRefine = 0, Refine4 = 1, Refine8 = 2, Refine32 = 3 };


struct EndParameters {
  double endTime;
};
} // namespace seissol::initializer::parameters

#endif
