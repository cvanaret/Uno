// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_DUALRESIDUALS_H
#define UNO_DUALRESIDUALS_H

#include "linear_algebra/Vector.hpp"
#include "tools/Infinity.hpp"

namespace uno {
   class DualResiduals {
   public:
      DualResiduals(size_t number_variables, size_t number_constraints):
         lagrangian_gradient(number_variables),
         constraints_buffer(number_constraints) {
      }

      double stationarity{INF<double>};
      double complementarity{INF<double>};

      double stationarity_scaling{INF<double>};
      double complementarity_scaling{INF<double>};

      Vector<double> lagrangian_gradient;
      Vector<double> constraints_buffer;
   };
} // namespace

#endif // UNO_DUALRESIDUALS_H
