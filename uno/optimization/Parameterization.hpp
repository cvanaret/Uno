// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PARAMETERIZATION_H
#define UNO_PARAMETERIZATION_H

#include <string_view>
#include <unordered_map>

namespace uno {
   class Parameterization {
   public:
      Parameterization() = default;

      void set(std::string_view key, double value);
      [[nodiscard]] double get(std::string_view key) const;

   protected:
      std::unordered_map<std::string_view, double> parameters;
   };
}

#endif // UNO_PARAMETERIZATION_H