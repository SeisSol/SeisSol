#include "dispatcher.hpp"
#include "DispatcherImplementation/implcpu.hpp"
#include "DispatcherImplementation/implgpu.hpp"
#include "Common/configs.hpp"

#include <variant>

namespace seissol::waveprop {
std::unique_ptr<WavePropDispatcherBase> getDispatcher(const seissol::initializers::LTS& lts,
                                                      seissol::initializers::Layer& layer,
                                                      const SupportedConfigs& config) {
#ifdef ACL_DEVICE
  return std::visit(
      [&](auto&& elem) {
        using ConfigT = std::decay_t<decltype(elem)>;
        auto ptr = std::make_unique<WavePropDispatcherGPU<ConfigT>>(lts, layer);
        return static_cast<std::unique_ptr<WavePropDispatcherBase>>(std::move(ptr));
      },
      config);
#else
  return std::visit(
      [&](auto&& elem) {
        using ConfigT = std::decay_t<decltype(elem)>;
        auto ptr = std::make_unique<WavePropDispatcherCPU<ConfigT>>(lts, layer);
        return static_cast<std::unique_ptr<WavePropDispatcherBase>>(std::move(ptr));
      },
      config);
#endif
}
} // namespace seissol::waveprop
