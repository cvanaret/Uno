// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNOSOLVERWRAPPER_H
#define UNO_UNOSOLVERWRAPPER_H

#include "Uno.hpp"
#include "options/Options.hpp"
#include "PythonModel.hpp"

namespace uno {
   class UnoSolverWrapper {
   public:
      Uno uno_solver{};
      Options options{};

      UnoSolverWrapper();

      [[nodiscard]] Result optimize(const PythonUserModel& user_model);
   };
} // namespace

#endif // UNO_UNOSOLVERWRAPPER_H