#include "dispatcher.hpp"
#include "DispatcherImplementation/implcpu.hpp"
#include "DispatcherImplementation/implgpu.hpp"

namespace seissol::waveprop {
    std::unique_ptr<WavePropDispatcherBase> getDispatcher(const seissol::initializers::LTS& lts, seissol::initializers::Layer& layer, bool plasticity) {
#ifdef ACL_DEVICE
        if (plasticity) {
            return std::make_unique<WavePropDispatcherGPU<true>>(lts, layer);
        }
        else {
            return std::make_unique<WavePropDispatcherGPU<false>>(lts, layer);
        }
#else
        if (plasticity) {
            return std::make_unique<WavePropDispatcherCPU<true>>(lts, layer);
        }
        else {
            return std::make_unique<WavePropDispatcherCPU<false>>(lts, layer);
        }
#endif
    }
}
