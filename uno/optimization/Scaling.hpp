#ifndef SCALING_H
#define SCALING_H

#include <vector>
#include "linear_algebra/SparseVector.hpp"

class Scaling {
public:
   Scaling(size_t number_constraints, double gradient_threshold);
   void compute(SparseVector<double>& objective_gradient, std::vector<SparseVector<double>>& constraints_jacobian);
   [[nodiscard]] double get_objective_scaling() const;
   [[nodiscard]] double get_constraint_scaling(size_t j) const;

protected:
   const double gradient_threshold;
   double objective_scaling;
   std::vector<double> constraints_scaling;
};

#endif // SCALING_H