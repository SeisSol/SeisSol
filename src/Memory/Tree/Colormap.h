// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_MEMORY_TREE_COLORMAP_H_
#define SEISSOL_SRC_MEMORY_TREE_COLORMAP_H_

#include <Common/Templating.h>
#include <Config.h>
#include <Initializer/BasicTypedefs.h>
#include <cstddef>
#include <functional>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace seissol::initializer {

// a wrapper type for both enums, but also static variants, like we use them
template <typename T>
class EnumLayer {
  public:
  EnumLayer(const std::vector<T>& supportedValues) : supportedValues(supportedValues) {
    for (std::size_t i = 0; i < supportedValues.size(); ++i) {
      reverse[supportedValues[i]] = i;
    }
  }

  [[nodiscard]] std::size_t size() const { return supportedValues.size(); }

  [[nodiscard]] T argument(std::size_t index) const { return supportedValues.at(index); }

  [[nodiscard]] std::size_t color(const T& value) const { return reverse.at(value); }

  using Type = T;

  using VariantType = void;

  private:
  std::vector<T> supportedValues;
  std::unordered_map<T, std::size_t> reverse;
};

template <typename T>
class TraitLayer {
  public:
  TraitLayer(const std::vector<T>& supportedValues) : supportedValues(supportedValues) {
    for (std::size_t i = 0; i < supportedValues.size(); ++i) {
      reverse[supportedValues[i].index()] = i;
    }
  }

  [[nodiscard]] std::size_t size() const { return supportedValues.size(); }

  [[nodiscard]] T argument(std::size_t index) const { return supportedValues.at(index); }

  [[nodiscard]] std::size_t color(const T& value) const { return reverse.at(value.index()); }

  using Type = T;

  using VariantType = T;

  private:
  std::vector<T> supportedValues;
  std::unordered_map<std::size_t, std::size_t> reverse;
};

class StopLayerSet {
  public:
  template <typename F, typename... Args>
  void call(std::size_t color, F&& func, Args... args) {
    std::invoke(std::forward<F>(func), args...);
  }

  [[nodiscard]] std::size_t color() const { return 0; }

  [[nodiscard]] std::size_t size() const { return 1; }

  using Type = std::tuple<>;

  using VariantType = std::variant<>;

  [[nodiscard]] Type argument(std::size_t color) const { return {}; }
};

template <typename Definition, typename SubLayerSet>
class LayerSet {
  public:
  LayerSet(Definition&& definition, SubLayerSet&& subLayerSet)
      : definition(std::move(definition)), subLayerSet(std::move(subLayerSet)) {}

  template <typename F, typename... Args>
  void call(std::size_t color, F&& func, Args... args) {
    const std::size_t index = color % definition.size();
    const std::size_t subColor = color / definition.size();
    subLayerSet.call(subColor, std::forward<F>(func), args..., definition.argument(index));
  }

  template <typename... RestArgs>
  [[nodiscard]] std::size_t color(const typename Definition::Type& value,
                                  const RestArgs&... rest) const {
    return definition.color(value) + definition.size() * subLayerSet.color(rest...);
  }

  [[nodiscard]] std::size_t size() const { return definition.size() * subLayerSet.size(); }

  using Type = PrependVariadicT<typename Definition::Type, typename SubLayerSet::Type>;

  using VariantType = std::conditional_t<
      std::is_same_v<typename Definition::VariantType, void>,
      typename SubLayerSet::VariantType,
      PrependVariadicT<typename Definition::VariantType, typename SubLayerSet::VariantType>>;

  [[nodiscard]] Type argument(std::size_t color) const {
    const std::size_t index = color % definition.size();
    const std::size_t subColor = color / definition.size();
    return std::tuple_cat(std::make_tuple(definition.argument(index)),
                          subLayerSet.argument(subColor));
  }

  private:
  Definition definition;
  SubLayerSet subLayerSet;
};

template <typename ConvenienceType, typename... Definitions>
class ColorMap {
  private:
  template <typename Head, typename... Rest>
  struct NestTypes {
    using Type = LayerSet<Head, typename NestTypes<Rest...>::Type>;
    static Type create(Head&& head, Rest&&... rest) {
      return Type(std::move(head), NestTypes<Rest...>::create(std::forward<Rest>(rest)...));
    }
  };

  template <typename Head>
  struct NestTypes<Head> {
    using Type = LayerSet<Head, StopLayerSet>;
    static Type create(Head&& head) { return Type(std::move(head), StopLayerSet()); }
  };
  using NestedLayerSets = typename NestTypes<Definitions...>::Type;
  NestedLayerSets layerSets;

  public:
  ColorMap(Definitions&&... definitions)
      : layerSets(NestTypes<Definitions...>::create(std::forward<Definitions>(definitions)...)) {}

  template <typename F, typename... Args>
  void call(int color, F&& func) const {
    layerSets.call(color, std::forward<F>(func));
  }

  template <typename... Args>
  [[nodiscard]] std::size_t color(const Args&... args) const {
    return layerSets.color(args...);
  }

  template <std::size_t... Idx>
  [[nodiscard]] std::size_t colorId(const ConvenienceType& data,
                                    std::index_sequence<Idx...> /*...*/) const {
    return color(std::get<Idx>(data.toTuple())...);
  }

  [[nodiscard]] std::size_t colorId(const ConvenienceType& data) const {
    return colorId(data, std::make_index_sequence<std::tuple_size_v<Type>>());
  }

  [[nodiscard]] std::size_t size() const { return layerSets.size(); }

  [[nodiscard]] std::vector<std::size_t> sizes() const { return layerSets.sizes(); }

  ConvenienceType argument(std::size_t color) const {
    const auto res = layerSets.argument(color);
    return ConvenienceType::fromTuple(res);
  }

  using Type = typename NestedLayerSets::Type;
};

struct LayerIdentifier {
  HaloType halo;
  ConfigVariant config;
  std::size_t lts;

  LayerIdentifier() = default;

  LayerIdentifier(HaloType halo, ConfigVariant config, std::size_t lts)
      : halo(halo), config(config), lts(lts) {}

  [[nodiscard]] std::tuple<ConfigVariant, std::size_t, HaloType> toTuple() const {
    return {config, lts, halo};
  }

  static LayerIdentifier fromTuple(const std::tuple<ConfigVariant, std::size_t, HaloType>& tuple) {
    return LayerIdentifier(std::get<2>(tuple), std::get<0>(tuple), std::get<1>(tuple));
  }
};

using LTSColorMap = ColorMap<LayerIdentifier,
                             TraitLayer<ConfigVariant>,
                             EnumLayer<std::size_t>,
                             EnumLayer<HaloType>>;

// to get a full layer representation (including Copy, Ghost, Interior), do
// ColorMap<EnumLayer<SupportedConfigs>, RangeLayer, EnumLayer<LayerType>> to also take into account
// different devices etc. (e.g. to get a single-process-per-node representation), define a variant
// std::variant<CpuType, GpuType>, and use an EnumLayer.

} // namespace seissol::initializer
#endif // SEISSOL_SRC_MEMORY_TREE_COLORMAP_H_
