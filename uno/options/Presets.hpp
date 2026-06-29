// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRESETS_H
#define UNO_PRESETS_H

#include <string>

namespace uno {
   // forward declarations
   class Model;
   class Options;

   enum class Preset {
      FILTERSQP = 0,
      IPOPT,
      FUNNELSQP,
      FILTERSLP,
   };

   class Presets {
   public:
      static void set(Options& options, const std::string& preset);
      static void set(const Model& model, Options& options, const std::string& preset);

   protected:
      [[nodiscard]] static Preset from_string(const std::string& preset);
      [[nodiscard]] static Preset pick_auto_preset(const Model& model);
      static void set(Options& options, Preset preset);
   };
} // namespace

#endif // UNO_PRESETS_H
