#pragma once

#include <Common/configs.hpp>
#include <Common/templating.hpp>
#include <functional>
#include <vector>
#include <utility>
#include <variant>

namespace seissol::initializers {

// a wrapper type for both enums, but also static variants, like we use them
template <typename T>
class EnumLayer {
  public:
  EnumLayer(const std::vector<T>& supportedValues) : supportedValues(supportedValues) {}

  int size() const { return supportedValues.size(); }

  int argument(int index) const { return supportedValues.at(index); }

  using Type = T;

  private:
  std::vector<T> supportedValues;
};

class RangeLayer {
  public:
  RangeLayer(int lower, int upper) : lower(lower), upper(upper) {}

  int size() const { return upper - lower + 1; }

  int argument(int index) const { return index + lower; }

  using Type = void;

  private:
  int lower, upper;
};

class StopLayerSet {
  public:
  template <typename F, typename... Args>
  void call(int color, F&& func, Args... args) {
    std::invoke(std::forward<F>(func), args...);
  }

  int color() { return 0; }

  int size() { return 1; }

  using Type = VariadicContainer<>;
};

template <typename Definition, typename SubLayerSet>
class LayerSet {
  public:
  LayerSet(Definition&& definition, SubLayerSet&& subLayerSet)
      : definition(definition), subLayerSet(subLayerSet) {}

  template <typename F, typename... Args>
  void call(int color, F&& func, Args... args) {
    int index = color % definition.size();
    int subColor = color / definition.size();
    subLayerSet.call(subColor, std::forward<F>(func), args..., definition.argument(index));
  }

  template <typename... RestArgs>
  int color(int index, RestArgs... rest) {
    return index + definition.size() * subLayerSet.color(rest...);
  }

  int size() { return definition.size() * subLayerSet.size(); }

  using Type =
      std::conditional_t<std::is_same_v<typename Definition::Type, void>,
                         typename SubLayerSet::Type,
                         PrependVariadicT<typename Definition::Type, typename SubLayerSet::Type>>;

  private:
  Definition definition;
  SubLayerSet subLayerSet;
};

template <typename... Definitions>
class ColorMap {
  private:
  template <typename Head, typename... Rest>
  struct NestTypes {
    using Type = LayerSet<Head, typename NestTypes<Rest...>::Type>;
    static void create(Head&& head, Rest&&... rest) {
      return Type(head, NestTypes<Rest...>::create(rest...));
    }
  };

  template <typename Head>
  struct NestTypes<Head> {
    using Type = LayerSet<Head, StopLayerSet>;
    static void create(Head&& head) { return Type(head, StopLayerSet()); }
  };
  using NestedLayerSets = typename NestTypes<Definitions...>::Type;
  NestedLayerSets LayerSets;

  public:
  ColorMap(Definitions&&... definitions) : LayerSets(NestedLayerSets::create(definitions...)) {}

  template <typename F, typename... Args>
  void call(int color, F&& func) {
    LayerSets.call(color, std::forward<F>(func));
  }

  template <typename... Args>
  int color(Args... args) {
    return LayerSets.color(args...);
  }

  int size() { return LayerSets.size(); }

  using Type = typename NestedLayerSets::Type;
};

// to replicate SeisSol v1.0.1 behavior, select ColorMap<RangeLayer> instead

// will yield Type = VariadicContainer<SupportedConfigs>
using ClusterColorMap = ColorMap<EnumLayer<SupportedConfigs>, RangeLayer>;

// to get a full layer representation (including Copy, Ghost, Interior), do
// ColorMap<EnumLayer<SupportedConfigs>, RangeLayer, EnumLayer<LayerType>> to also take into account
// different devices etc. (e.g. to get a single-process-per-node representation), define a variant
// std::variant<CpuType, GpuType>, and use an EnumLayer.

} // namespace seissol::initializers
