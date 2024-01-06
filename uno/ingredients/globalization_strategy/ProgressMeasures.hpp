// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PROGRESSMEASURES_H
#define UNO_PROGRESSMEASURES_H

#include <functional>

struct ProgressMeasures {
   double infeasibility{}; // constraint violation
   std::function<double(double objective_multiplier)> optimality{}; // optimality measure (scaled by penalty parameter): objective, Lagrangian
   double auxiliary_terms{}; // auxiliary terms (independent of penalty parameter): barrier terms, proximal term, ...
};

#endif // UNO_PROGRESSMEASURES_H
