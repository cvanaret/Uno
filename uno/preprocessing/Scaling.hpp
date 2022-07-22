// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SCALING_H
#define UNO_SCALING_H

#include <vector>
#include "linear_algebra/SparseVector.hpp"

class Scaling {
public:
   Scaling(size_t number_constraints, double gradient_threshold);
   void compute(SparseVector<double>& objective_gradient, std::vector<SparseVector<double>>& constraint_jacobian);
   [[nodiscard]] double get_objective_scaling() const;
   [[nodiscard]] double get_constraint_scaling(size_t j) const;

protected:
   const double gradient_threshold;
   double objective_scaling;
   std::vector<double> constraint_scaling;
};

#endif // UNO_SCALING_H