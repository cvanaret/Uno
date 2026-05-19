// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRESETS_H
#define UNO_PRESETS_H

#include <string>

namespace uno {
   // forward declaration
   class Options;

   class Presets {
   public:
      static void set_default(Options& options);
      static void set(Options& options, const std::string& preset_name);
   };
} // namespace

#endif // UNO_PRESETS_H
