#pragma once

#include <Common/configs.hpp>
#include <Initializer/Boundary.h>
#include <Initializer/DynamicRupture.h>
#include <variant>
#include <memory>

#include "LTSTree.hpp"
#include "Initializer/LTS.h"
#include "utils/logger.h"

namespace seissol {
namespace initializers {
template <typename LtsT>
struct LTSRef {
  LTSTree& tree;
  const LtsT& lts;
  std::size_t id;

  LTSRef(std::size_t id, LTSTree& tree)
      : id(id), tree(tree), lts(dynamic_cast<const LtsT&>(tree.lts())) {}

  template <typename T>
  T* var(const Variable<T>& variable) {
    return tree.var(variable);
  }
};

template <typename LtsT>
struct LTSLayerRef {
  LTSTree& tree;
  Layer& layer;
  const LtsT& lts;
  std::size_t config;
  std::size_t cluster;
  LayerType icg;

  LTSLayerRef(std::size_t config, std::size_t cluster, LayerType icg, LTSTree& tree, Layer& layer)
      : config(config), cluster(cluster), icg(icg), tree(tree), layer(layer),
        lts(dynamic_cast<const LtsT&>(tree.lts())) {}

  template <typename T>
  T* var(const Variable<T>& variable) {
    return layer.var(variable);
  }
};

template <typename T, template <typename> typename LtsT>
class LTSForest {
  public:
  void initialize(unsigned numberTimeClusters) {
    initializeInternal<0>();
    visit([&](LTSTree& tree, const auto& lts) {
      tree.setNumberOfTimeClusters(numberTimeClusters);
      tree.fixate();
    });
  }

  // F == void(LtsRef<Config>), where Config is from SupportedConfigs
  template <typename F>
  constexpr void visit(F&& visitor) {
    visitInternal<0>(std::forward<F>(visitor));
  }

  // F == void(LtsRef<Config>), where Config is from SupportedConfigs
  template <typename F>
  constexpr void visitIdx(std::size_t idx, F&& visitor) {
    visitIdxInternal<0>(idx, std::forward<F>(visitor));
  }

  // F == void(int, int, LayerType, LtsRef<Config>), where Config is from SupportedConfigs
  template <typename F>
  constexpr void visitLayers(F&& visitor) {
    visitLayersInternal<0>(std::forward<F>(visitor));
  }

  // F == void(LtsRef<Lts<Config>>, LtsRef<Lts2<Config>>), where Config is from SupportedConfigs
  template <typename F, template <typename> typename Lts2T>
  constexpr void visitTwo(F&& visitor, LTSForest<T, Lts2T>& two) {
    visitTwoInternal<0, F, Lts2T>(std::forward<F>(visitor), two);
  }

  template <typename F, template <typename> typename Lts2T>
  constexpr void visitTwoLayers(F&& visitor, LTSForest<T, Lts2T>& two) {
    visitTwoLayersInternal<0, F, Lts2T>(std::forward<F>(visitor), two);
  }

  void allocateTouchVariables() {
    visit([](auto&& ltsview) {
      ltsview.tree.allocateVariables();
      ltsview.tree.touchVariables();
    });
  }

  void allocateBuckets() {
    visit([](auto&& ltsview) { ltsview.tree.allocateBuckets(); });
  }

  std::array<LTSTree, std::variant_size_v<T>> children;

  private:
  template <std::size_t Idx>
  constexpr void initializeInternal() {
    if constexpr (Idx < std::variant_size_v<T>) {
      using LTSType = LtsT<std::variant_alternative_t<Idx, T>>;
      children[Idx].attachLTS(std::forward<LTSType>(LTSType()));
      initializeInternal<Idx + 1>();
    }
  }

  template <std::size_t Idx, typename F>
  constexpr void visitInternal(F&& visitor) {
    if constexpr (Idx < std::variant_size_v<T>) {
      using LTSType = LtsT<std::variant_alternative_t<Idx, T>>;
      visitor(LTSRef<LTSType>(Idx, children[Idx]));
      visitInternal<Idx + 1>(std::forward<F>(visitor));
    }
  }

  template <std::size_t Idx, typename F>
  constexpr void visitIdxInternal(std::size_t idx, F&& visitor) {
    if constexpr (Idx < std::variant_size_v<T>) {
      if (Idx == idx) {
        using LTSType = LtsT<std::variant_alternative_t<Idx, T>>;
        visitor(LTSRef<LTSType>(Idx, children[Idx]));
      } else {
        visitIdxInternal<Idx + 1>(idx, std::forward<F>(visitor));
      }
    } else {
      logError() << "Invalid index" << idx << "/" << Idx;
    }
  }

  template <std::size_t Idx, typename F>
  constexpr void visitLayersInternal(F&& visitor) {
    if constexpr (Idx < std::variant_size_v<T>) {
      using LTSType = LtsT<std::variant_alternative_t<Idx, T>>;
      for (std::size_t i = 0; i < children[Idx].numChildren(); ++i) {
        auto& child = children[Idx].child(i);
        visitor(std::forward<LTSLayerRef<LTSType>>(
            LTSLayerRef<LTSType>(Idx, i, Ghost, children[Idx], child.template child<Ghost>())));
        visitor(std::forward<LTSLayerRef<LTSType>>(
            LTSLayerRef<LTSType>(Idx, i, Copy, children[Idx], child.template child<Copy>())));
        visitor(std::forward<LTSLayerRef<LTSType>>(LTSLayerRef<LTSType>(
            Idx, i, Interior, children[Idx], child.template child<Interior>())));
      }
      visitLayersInternal<Idx + 1>(std::forward<F>(visitor));
    }
  }

  template <std::size_t Idx, typename F, template <typename> typename Lts2T>
  constexpr void visitTwoInternal(F&& visitor, LTSForest<T, Lts2T>& two) {
    if constexpr (Idx < std::variant_size_v<T>) {
      using LTSType = LtsT<std::variant_alternative_t<Idx, T>>;
      using LTS2Type = Lts2T<std::variant_alternative_t<Idx, T>>;
      visitor(LTSRef<LTSType>(children[Idx]), LTSRef<LTS2Type>(two.children[Idx]));
      visitTwoInternal<Idx + 1, F, Lts2T>(std::forward<F>(visitor), two);
    }
  }

  template <std::size_t Idx, typename F, template <typename> typename Lts2T>
  constexpr void visitTwoLayersInternal(F&& visitor, LTSForest<T, Lts2T>& two) {
    if constexpr (Idx < std::variant_size_v<T>) {
      using LTSType = LtsT<std::variant_alternative_t<Idx, T>>;
      using LTS2Type = Lts2T<std::variant_alternative_t<Idx, T>>;
      for (std::size_t i = 0; i < children[Idx].numChildren(); ++i) {
        auto& child = children[Idx].child(i);
        visitor(LTSLayerRef<LTSType>(Idx, i, Ghost, children[Idx], child.template child<Ghost>()),
                LTSLayerRef<LTS2Type>(Idx, i, Ghost, children[Idx], child.template child<Ghost>()));
        visitor(LTSLayerRef<LTSType>(Idx, i, Copy, children[Idx], child.template child<Copy>()),
                LTSLayerRef<LTS2Type>(Idx, i, Copy, children[Idx], child.template child<Copy>()));
        visitor(
            LTSLayerRef<LTSType>(Idx, i, Interior, children[Idx], child.template child<Interior>()),
            LTSLayerRef<LTS2Type>(
                Idx, i, Interior, children[Idx], child.template child<Interior>()));
      }
      visitTwoLayersInternal<Idx + 1, F, Lts2T>(std::forward<F>(visitor), two);
    }
  }
};

using ClusterLTSForest = LTSForest<SupportedConfigs, LTS>;
using DynRupLTSForest = LTSForest<SupportedConfigs, DynamicRupture>;
using BoundaryLTSForest = LTSForest<SupportedConfigs, Boundary>;
} // namespace initializers
} // namespace seissol

// a generalized variant of LTSForest should take variadic templates
// also introduce untyped layers
