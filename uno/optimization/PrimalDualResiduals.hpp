// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALRESIDUALS_H
#define UNO_PRIMALDUALRESIDUALS_H

#include "tools/Infinity.hpp"

class PrimalDualResiduals {
public:
   double stationarity{INF<double>};
   double FJ_stationarity{INF<double>};
   double feasibility_stationarity{INF<double>};
   double infeasibility{INF<double>};
   double complementarity{INF<double>};
   double feasibility_complementarity{INF<double>};
   double stationarity_scaling{INF<double>};
   double complementarity_scaling{INF<double>};
};

#endif // UNO_PRIMALDUALRESIDUALS_H