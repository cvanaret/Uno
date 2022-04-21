#include "Scaling.hpp"
#include "linear_algebra/Vector.hpp"

Scaling::Scaling(size_t number_constraints, double gradient_threshold):
      gradient_threshold(gradient_threshold),
      objective_scaling(1.),
      constraint_scaling(number_constraints, 1.) {
}

void Scaling::compute(SparseVector<double>& objective_gradient, std::vector<SparseVector<double>>& constraint_jacobian) {
   // set the objective scaling
   this->objective_scaling = std::min(1., this->gradient_threshold / norm_inf(objective_gradient));

   // set the constraints scaling
   for (size_t j = 0; j < this->constraint_scaling.size(); j++) {
      this->constraint_scaling[j] = std::min(1., this->gradient_threshold / norm_inf(constraint_jacobian[j]));
   }

   // scale the gradients passed as parameters
   scale(objective_gradient, this->get_objective_scaling());
   for (size_t j = 0; j < this->constraint_scaling.size(); j++) {
      scale(constraint_jacobian[j], this->get_constraint_scaling(j));
   }
   DEBUG << "Objective scaling: " << this->objective_scaling << '\n';
   DEBUG << "Constraint scaling: "; print_vector(DEBUG, this->constraint_scaling);
}

double Scaling::get_objective_scaling() const {
   return this->objective_scaling;
}

double Scaling::get_constraint_scaling(size_t j) const {
   return this->constraint_scaling[j];
}