// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

// re2c -o src/JsonLexer.cpp src/JsonLexer.re --case-ranges

#include "JsonLexer.hpp"

#include <cerrno>
#include <cstdlib>

namespace seissol {

JsonLexer::JsonLexer(char const* input) : YYCURSOR{input} {}

JsonTokenKind JsonLexer::operator()(JsonToken& tok) {
  char const* YYMARKER;
  JsonLexerError errc = JsonLexerError::None;
lex:
  char const* b = YYCURSOR;
  /*!re2c
      re2c:yyfill:enable = 0;
      re2c:define:YYCTYPE = char;

      // regexp
      whitespace            = [ \t\v\r\n]+;

      digit                 = [0-9];
      integer_number        = "-"? ("0" | [1-9] digit*);
      fractional            = "." digit+;
      exponent              = [eE] [-+]? digit+;
      floating_number       = integer_number fractional? exponent?;

      char                  = [^\"] | "\\\"";
      string                = "\"" char* "\"";

      // structural
      "["                 { return JsonTokenKind::BeginArray; }
      "]"                 { return JsonTokenKind::EndArray; }
      "{"                 { return JsonTokenKind::BeginObject; }
      "}"                 { return JsonTokenKind::EndObject; }
      ":"                 { return JsonTokenKind::NameSeparator; }
      ","                 { return JsonTokenKind::ValueSeparator; }

      // value
      "false"             { return JsonTokenKind::False; }
      "true"              { return JsonTokenKind::True; }
      "null"              { return JsonTokenKind::Null; }
      integer_number      {
          tok.val = lexIntegerNumber(b, YYCURSOR, errc);
          if (errc != JsonLexerError::None) {
              tok.val = errc;
              return JsonTokenKind::Unknown;
          }
          return JsonTokenKind::IntegerNumber;
      }
      floating_number     {
          tok.val = lexFloatingNumber(b, YYCURSOR, errc);
          if (errc != JsonLexerError::None) {
              tok.val = errc;
              return JsonTokenKind::Unknown;
          }
          return JsonTokenKind::FloatingNumber;
      }
      string              {
          tok.val = lexString(b + 1, YYCURSOR, errc);
          if (errc != JsonLexerError::None) {
              tok.val = errc;
              return JsonTokenKind::Unknown;
          }
          return JsonTokenKind::String;
      }

      whitespace          { goto lex; }
      [\x00]              { return JsonTokenKind::EOI; }
      *                   { tok.val = JsonLexerError::None; return JsonTokenKind::Unknown; }
   */
}

long JsonLexer::lexIntegerNumber(char const* s, char const* e, JsonLexerError& errc) {
  auto number = std::string(s, e);
  double d = strtol(number.c_str(), nullptr, 10);
  errc = errno == ERANGE ? JsonLexerError::IntegerOverflow : JsonLexerError::None;
  return d;
}

double JsonLexer::lexFloatingNumber(char const* s, char const* e, JsonLexerError& errc) {
  auto number = std::string(s, e);
  double d = strtod(number.c_str(), nullptr);
  errc = errno == ERANGE ? JsonLexerError::FloatingOutOfRange : JsonLexerError::None;
  return d;
}

std::string JsonLexer::lexString(char const* s, char const* e, JsonLexerError& errc) {
  char const* YYMARKER;
  auto str = std::string{};
  str.reserve(e - s);
  errc = JsonLexerError::None;
lex:
  char const* b = s;
  /*!re2c
      re2c:yyfill:enable = 0;
      re2c:define:YYCURSOR = s;

      unescaped             = [\x20-\x21] | [\x23-\x5B] | [\x5D-\x7F];
      escaped               = "\\" [\x22\x5C\x2F\x62\x66\x6E\x72\x74];

      unescaped  { str.push_back(*b); goto lex; }
      escaped    {
          ++b; // skip escape char
          char c;
          switch(*b) {
          case 0x22:
          case 0x5C:
          case 0x2F:
              c = *b;
              break;
          case 0x62: c = 0x08; break;
          case 0x66: c = 0x0C; break;
          case 0x6E: c = 0x0A; break;
          case 0x72: c = 0x0D; break;
          case 0x74: c = 0x09; break;
          case 0x64: c = 0x22; break;
          };
          str.push_back(c);
          goto lex;
      }
      "\""       { return str; }
      *          { errc = JsonLexerError::InvalidString; return {}; }
  */
  return str;
}

} // namespace seissol
