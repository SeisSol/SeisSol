#ifndef SRC_SOLVER_MULTIPLESIMULATIONS_H_
#define SRC_SOLVER_MULTIPLESIMULATIONS_H_

namespace seissol::multiplesimulations {
#ifdef MULTIPLE_SIMULATIONS
constexpr unsigned int NumSimulations = MULTIPLE_SIMULATIONS;
constexpr unsigned int BasisFunctionDimension = 1;
template <typename F, typename... Args>
auto multisimWrap(F&& function, size_t sim, Args&&... args) {
  return std::invoke(std::forward<F>(function), sim, std::forward<Args>(args)...);
}
template<typename T, typename F, typename ...Args>
void multisimObjectWrap(F&& func, T& obj, int sim, Args&&... args) {
    std::invoke(std::forward<F>(func), obj, sim, std::forward<Args>(args)...);
}
constexpr size_t MultisimStart = init::QAtPoint::Start[0];
constexpr size_t MultisimEnd = init::QAtPoint::Stop[0];
constexpr bool MultipleSimulations = true;
#else
constexpr unsigned int NumSimulations = 1;
constexpr unsigned int BasisFunctionDimension = 0;
template <typename F, typename... Args>
auto multisimWrap(F&& function, size_t sim, Args&&... args) {
  return std::invoke(std::forward<F>(function), std::forward<Args>(args)...);
}
template<typename T, typename F, typename ...Args>
void multisimObjectWrap(F&& func, T& obj, int sim, Args&&... args) {
    std::invoke(std::forward<F>(func), obj, std::forward<Args>(args)...);
}
constexpr size_t MultisimStart = 0;
constexpr size_t MultisimEnd = 1;
constexpr bool MultipleSimulations = false;
#endif

} // namespace seissol::multiplesimulations

#endif
