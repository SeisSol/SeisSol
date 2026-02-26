// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_MEMORY_TREE_COLORMAP_H_
#define SEISSOL_SRC_MEMORY_TREE_COLORMAP_H_

#include "Common/Templating.h"
#include "Config.h"
#include "Initializer/BasicTypedefs.h"

#include <cstddef>
#include <functional>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace seissol::initializer {

/**
  A wrapper type for enums and integer lists, to be used in `ColorMap`.
 */
template <typename T>
class EnumLayer {
  public:
  explicit EnumLayer(const std::vector<T>& supportedValues) : supportedValues_(supportedValues) {
    for (std::size_t i = 0; i < supportedValues.size(); ++i) {
      reverse_[supportedValues[i]] = i;
    }
  }

  [[nodiscard]] std::size_t size() const { return supportedValues_.size(); }

  [[nodiscard]] T argument(std::size_t index) const { return supportedValues_.at(index); }

  [[nodiscard]] std::size_t color(const T& value) const { return reverse_.at(value); }

  using Type = T;

  using VariantType = void;

  private:
  std::vector<T> supportedValues_;
  std::unordered_map<T, std::size_t> reverse_;
};

/**
  A wrapper type around a `std::variant`, to be used in `ColorMap`.
 */
template <typename T>
class TraitLayer {
  public:
  explicit TraitLayer(const std::vector<T>& supportedValues) : supportedValues_(supportedValues) {
    for (std::size_t i = 0; i < supportedValues.size(); ++i) {
      reverse_[supportedValues[i].index()] = i;
    }
  }

  [[nodiscard]] std::size_t size() const { return supportedValues_.size(); }

  [[nodiscard]] T argument(std::size_t index) const { return supportedValues_.at(index); }

  [[nodiscard]] std::size_t color(const T& value) const { return reverse_.at(value.index()); }

  using Type = T;

  using VariantType = T;

  private:
  std::vector<T> supportedValues_;
  std::unordered_map<std::size_t, std::size_t> reverse_;
};

/**
  Helper class for `ColorMap`; indicates that we have reached the end of the `LayerSet`.
  Not intended for outside use.
 */
class StopLayerSet {
  public:
  template <typename F, typename... Args>
  void call(std::size_t /*color*/, F&& func, Args... args) {
    std::invoke(std::forward<F>(func), args...);
  }

  [[nodiscard]] std::size_t color() const { return 0; }

  [[nodiscard]] std::size_t size() const { return 1; }

  using Type = std::tuple<>;

  using VariantType = std::variant<>;

  [[nodiscard]] Type argument(std::size_t /*color*/) const { return {}; }
};

/**
  Helper class for `ColorMap`; contains either an `EnumLayer` or a `TraitLayer`.
  Not intended for outside use.
 */
template <typename Definition, typename SubLayerSet>
class LayerSet {
  public:
  LayerSet(Definition&& definition, SubLayerSet&& subLayerSet)
      : definition_(std::move(definition)), subLayerSet_(std::move(subLayerSet)) {}

  template <typename F, typename... Args>
  void call(std::size_t color, F&& func, Args... args) {
    const std::size_t index = color % definition_.size();
    const std::size_t subColor = color / definition_.size();
    subLayerSet_.call(subColor, std::forward<F>(func), args..., definition_.argument(index));
  }

  template <typename... RestArgs>
  [[nodiscard]] std::size_t color(const typename Definition::Type& value,
                                  const RestArgs&... rest) const {
    return definition_.color(value) + definition_.size() * subLayerSet_.color(rest...);
  }

  [[nodiscard]] std::size_t size() const { return definition_.size() * subLayerSet_.size(); }

  using Type = PrependVariadicT<typename Definition::Type, typename SubLayerSet::Type>;

  using VariantType = std::conditional_t<
      std::is_same_v<typename Definition::VariantType, void>,
      typename SubLayerSet::VariantType,
      PrependVariadicT<typename Definition::VariantType, typename SubLayerSet::VariantType>>;

  [[nodiscard]] Type argument(std::size_t color) const {
    const std::size_t index = color % definition_.size();
    const std::size_t subColor = color / definition_.size();
    return std::tuple_cat(std::make_tuple(definition_.argument(index)),
                          subLayerSet_.argument(subColor));
  }

  private:
  Definition definition_;
  SubLayerSet subLayerSet_;
};

/**
  Defines an ordering for different storage layers in a Storage structure,
  in a hierarchitcal fashion; with a convenience type (instead of bare tuples).
 */
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
  NestedLayerSets layerSets_;

  public:
  explicit ColorMap(Definitions&&... definitions)
      : layerSets_(NestTypes<Definitions...>::create(std::forward<Definitions>(definitions)...)) {}

  template <typename F, typename... Args>
  void call(int color, F&& func) const {
    layerSets_.call(color, std::forward<F>(func));
  }

  template <typename... Args>
  [[nodiscard]] std::size_t color(const Args&... args) const {
    return layerSets_.color(args...);
  }

  template <std::size_t... Idx>
  [[nodiscard]] std::size_t colorId(const ConvenienceType& data,
                                    std::index_sequence<Idx...> /*...*/) const {
    return color(std::get<Idx>(data.toTuple())...);
  }

  [[nodiscard]] std::size_t colorId(const ConvenienceType& data) const {
    return colorId(data, std::make_index_sequence<std::tuple_size_v<Type>>());
  }

  [[nodiscard]] std::size_t size() const { return layerSets_.size(); }

  [[nodiscard]] std::vector<std::size_t> sizes() const { return layerSets_.sizes(); }

  ConvenienceType argument(std::size_t color) const {
    const auto res = layerSets_.argument(color);
    return ConvenienceType::fromTuple(res);
  }

  using Type = typename NestedLayerSets::Type;
};

using ConfigVariant = std::variant<Config>;

/**
  A convenience data structure for the Layer identifier type we use.
  In essence, it provides basic infos about the layer, like e.g. if it is a ghost, copy or interior
  layer or which LTS cluster it belongs to; and maps identifier to a contiguous ID (a.k.a. color),
  and back.
  */
struct LayerIdentifier {
  HaloType halo{HaloType::Interior};
  ConfigVariant config;
  std::size_t lts{};

  LayerIdentifier() = default;

  LayerIdentifier(HaloType halo, ConfigVariant config, std::size_t lts)
      : halo(halo), config(config), lts(lts) {}

  [[nodiscard]] std::tuple<HaloType, std::size_t, ConfigVariant> toTuple() const {
    return {halo, lts, config};
  }

  static LayerIdentifier fromTuple(const std::tuple<HaloType, std::size_t, ConfigVariant>& tuple) {
    return LayerIdentifier(std::get<0>(tuple), std::get<2>(tuple), std::get<1>(tuple));
  }
};

// NOTE: keep the ordering like this until we merge #1411. Otherwise, there will be a data layout
// mismatch with the LtsLayout class.
using LTSColorMap = ColorMap<LayerIdentifier,
                             EnumLayer<HaloType>,
                             EnumLayer<std::size_t>,
                             TraitLayer<ConfigVariant>>;

// to get a full layer representation (including Copy, Ghost, Interior), do
// ColorMap<EnumLayer<SupportedConfigs>, RangeLayer, EnumLayer<LayerType>> to also take into account
// different devices etc. (e.g. to get a single-process-per-node representation), define a variant
// std::variant<CpuType, GpuType>, and use an EnumLayer.

} // namespace seissol::initializer
#endif // SEISSOL_SRC_MEMORY_TREE_COLORMAP_H_
