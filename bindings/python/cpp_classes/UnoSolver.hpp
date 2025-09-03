// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_UNOSOLVER_H
#define UNO_UNOSOLVER_H

#include "Uno.hpp"
#include "model/Model.hpp"
#include "optimization/Result.hpp"
#include "options/Options.hpp"

namespace uno {
   class UnoSolver {
   public:
      Uno uno_solver{};
      Options options{};

      UnoSolver();

      void optimize(const Model& model);
   };
} // namespace

#endif // UNO_UNOSOLVER_H