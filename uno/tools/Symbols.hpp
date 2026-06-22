// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMBOLS_H
#define UNO_SYMBOLS_H

#include <string_view>

namespace uno {
   namespace symbols {
      constexpr bool use_ascii =
#ifdef _WIN32
         true;
#else
         false;
#endif

      // picks the ASCII fallback on Windows, the UTF-8 glyph elsewhere
      constexpr std::string_view glyph(std::string_view utf8, std::string_view ascii) {
         return use_ascii ? ascii : utf8;
      }

      static constexpr std::string_view hyphen = glyph(u8"\u2500", "-"); // ─

      static constexpr std::string_view pipe = glyph(u8"\u2502", "|"); // │
      static constexpr std::string_view top_pipe = glyph(u8"\u250C", "|");  // ┌
      static constexpr std::string_view bottom_pipe = glyph(u8"\u2514", "|");  // └

      static constexpr std::string_view top_left_corner = glyph(u8"\u250C", "+");  // ┌
      static constexpr std::string_view top_tee = glyph(u8"\u252C", "+");  // ┬
      static constexpr std::string_view top_right_corner = glyph(u8"\u2510", "+");  // ┐
      static constexpr std::string_view left_tee = glyph(u8"\u251C", "+");  // ├
      static constexpr std::string_view cross = glyph(u8"\u253C", "+");  // ┼
      static constexpr std::string_view right_tee = glyph(u8"\u2524", "+");  // ┤
      static constexpr std::string_view bottom_left_corner = glyph(u8"\u2514", "+");  // └
      static constexpr std::string_view bottom_tee = glyph(u8"\u2534", "+");  // ┴
      static constexpr std::string_view bottom_right_corner = glyph(u8"\u2518", "+");  // ┘

      static constexpr std::string_view check = glyph(u8"\u2714", "v"); // ✔
      static constexpr std::string_view fail = glyph(u8"\u2718", "x"); // ✘
   }
} // namespace

#endif // UNO_SYMBOLS_H