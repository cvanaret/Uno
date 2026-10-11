// Copyright (c) 2018-2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_RESIDUALS_H
#define UNO_RESIDUALS_H

#include "linear_algebra/Vector.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   class Residuals {
   public:
      Residuals(size_t number_variables):
         lagrangian_gradient(number_variables) {
      }

      double primal_infeasibility{Inf};
      double stationarity{Inf};
      double complementarity{Inf};

      double stationarity_scaling{Inf};
      double complementarity_scaling{Inf};

      Vector<double> lagrangian_gradient;
   };
} // namespace

#endif // UNO_RESIDUALS_H