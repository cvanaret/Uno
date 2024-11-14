// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DEFAULTOPTIONS_H
#define UNO_DEFAULTOPTIONS_H

#include <optional>
#include <string>
#include "Options.hpp"

namespace uno {
   class DefaultOptions {
   public:
      [[nodiscard]] static Options load();
      [[nodiscard]] static Options determine_solvers();
   };
} // namespace

#endif // UNO_DEFAULTOPTIONS_H
