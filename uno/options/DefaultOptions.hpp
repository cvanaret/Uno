// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DEFAULTOPTIONS_H
#define UNO_DEFAULTOPTIONS_H

#include "Options.hpp"

namespace uno {
   class DefaultOptions {
   public:
      static void load(Options& options);

   protected:
      static void determine_subproblem_solvers(Options& options);
   };
} // namespace

#endif // UNO_DEFAULTOPTIONS_H