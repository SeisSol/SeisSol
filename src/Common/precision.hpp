
namespace seissol {
    template<typename RealT>
    struct RealTypeInfo {};

    template<>
    struct RealTypeInfo<float> {
        using Type = float;
        static constexpr RealType RealType = RealType::F32;
    };

    template<>
    struct RealTypeInfo<double> {
        using Type = double;
        static constexpr RealType RealType = RealType::F64;
    };
}
