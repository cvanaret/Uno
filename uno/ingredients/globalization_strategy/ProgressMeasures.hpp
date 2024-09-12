// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PROGRESSMEASURES_H
#define UNO_PROGRESSMEASURES_H

#include <functional>
#include "tools/Infinity.hpp"

namespace uno {
   struct ProgressMeasures {
      double infeasibility{}; // constraint violation
      std::function<double(double objective_multiplier)> objective{}; // objective measure (scaled by penalty parameter): objective, Lagrangian
      double auxiliary{}; // auxiliary terms (independent of penalty parameter): barrier terms, proximal term, ...

      void reset() {
         this->infeasibility = INF<double>;
         this->objective = [](double) { return INF<double>; };
         this->auxiliary = INF<double>;
      }
   };
} // namespace

#endif // UNO_PROGRESSMEASURES_H
