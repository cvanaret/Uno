// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALRESIDUALS_H
#define UNO_PRIMALDUALRESIDUALS_H

#include "tools/Infinity.hpp"

namespace uno {
   class PrimalDualResiduals {
   public:
      double stationarity{INF<double>};
      double primal_feasibility{INF<double>};
      double dual_feasibility{INF<double>};
      double complementarity{INF<double>};

      double stationarity_scaling{INF<double>};
      double complementarity_scaling{INF<double>};
   };
} // namespace

#endif // UNO_PRIMALDUALRESIDUALS_H
