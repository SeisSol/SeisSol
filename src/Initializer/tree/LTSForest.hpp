#pragma once

#include <Common/configs.hpp>
#include <Initializer/Boundary.h>
#include <Initializer/DynamicRupture.h>
#include <variant>
#include <memory>

#include "LTSTree.hpp"
#include "Initializer/LTS.h"

namespace seissol {
    namespace initializers{
        template<typename LtsT>
        struct LTSRef {
            LTSTree& tree;
            const LtsT& lts;

            LTSRef(LTSTree& tree) : tree(tree), lts(dynamic_cast<const LtsT&>(tree.lts())) {}

            template<typename T>
            T* var(const Variable<T>& variable) {
                return tree.var(variable);
            }
        };

        template<typename T, template<typename> typename LtsT>
        class LTSForest {
        public:
            void initialize(unsigned numberTimeClusters) {
                initializeInternal<0>();
                visit([&](LTSTree& tree, const auto& lts){
                    tree.setNumberOfTimeClusters(numberTimeClusters);
                    tree.fixate();
                });
            }

            // F == void(LTSTree, LTS<typed>)
            template<typename F>
            constexpr void visit(F&& visitor) {
                visitInternal<0>(std::forward<F>(visitor));
            }

            template<typename F, template<typename> typename Lts2T>
            constexpr void visitTwo(F&& visitor, LTSForest<T, Lts2T>& two) {
                visitTwoInternal<0, F, Lts2T>(std::forward<F>(visitor), two);
            }
            
            std::array<LTSTree, std::variant_size_v<T>> children;
        private:

            template<std::size_t Idx>
            constexpr void initializeInternal() {
                if constexpr(Idx < std::variant_size_v<T>) {
                    using LTSType = LtsT<std::variant_alternative_t<Idx, T>>;
                    children[Idx].attachLTS(std::forward<LTSType>(LTSType()));
                    initializeInternal<Idx+1>();
                }
            }

            template<std::size_t Idx, typename F>
            constexpr void visitInternal(F&& visitor) {
                if constexpr(Idx < std::variant_size_v<T>) {
                    using LTSType = LtsT<std::variant_alternative_t<Idx, T>>;
                    visitor(LTSRef<LTSType>(children[Idx]));
                    visitInternal<Idx+1>(std::forward<F>(visitor));
                }
            }

            template<std::size_t Idx, typename F, template<typename> typename Lts2T>
            constexpr void visitTwoInternal(F&& visitor, LTSForest<T, Lts2T>& two) {
                if constexpr(Idx < std::variant_size_v<T>) {
                    using LTSType = LtsT<std::variant_alternative_t<Idx, T>>;
                    using LTS2Type = Lts2T<std::variant_alternative_t<Idx, T>>;
                    visitor(LTSRef<LTSType>(children[Idx]), LTSRef<LTS2Type>(two.children[Idx]));
                    visitTwoInternal<Idx+1, F, Lts2T>(std::forward<F>(visitor), two);
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
