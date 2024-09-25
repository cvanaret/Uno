// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALRESIDUALS_H
#define UNO_PRIMALDUALRESIDUALS_H

#include "tools/Infinity.hpp"

namespace uno {
   class PrimalDualResiduals {
   public:
      explicit PrimalDualResiduals(size_t number_variables): lagrangian_gradient(number_variables) { }

      double stationarity{INF<double>};
      double primal_feasibility{INF<double>};
      double complementarity{INF<double>};

      double stationarity_scaling{INF<double>};
      double complementarity_scaling{INF<double>};

      LagrangianGradient<double> lagrangian_gradient;
   };
} // namespace

#endif // UNO_PRIMALDUALRESIDUALS_H
