#ifndef MODEL_VISITOR_H_
#define MODEL_VISITOR_H_

namespace seissol {
  template<typename... Types>
  struct make_visitor : Types... {
    using Types::operator()...;
  };
  // make_visitor template deduction guide
  template<typename... Types>
  make_visitor(Types...) -> make_visitor<Types...>;
}

#endif
