// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef JSONLEXER_20231019_HPP
#define JSONLEXER_20231019_HPP

#include "JsonToken.hpp"

#include <string>

namespace seissol {

class JsonLexer {
  public:
  JsonLexer(char const* input);
  JsonTokenKind operator()(JsonToken& tok);

  private:
  long lexIntegerNumber(char const* s, char const* e, JsonLexerError& errc);
  double lexFloatingNumber(char const* s, char const* e, JsonLexerError& errc);
  std::string lexString(char const* s, char const* e, JsonLexerError& errc);

  char const* YYCURSOR;
};

} // namespace seissol

#endif // JSONLEXER_20231019_HPP

