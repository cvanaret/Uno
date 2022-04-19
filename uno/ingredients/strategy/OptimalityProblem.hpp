#ifndef UNO_OPTIMALITYPROBLEM_H
#define UNO_OPTIMALITYPROBLEM_H

#include <vector>
#include <cmath>
#include "NonlinearProblem.hpp"
#include "ingredients/constraint_relaxation/ElasticVariables.hpp"

class OptimalityProblem: public NonlinearProblem {
public:
   explicit OptimalityProblem(const Model& model);

   [[nodiscard]] double get_objective_multiplier() const override;
   [[nodiscard]] double evaluate_objective(Iterate& iterate) const override;
   void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const override;
   void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const override;
   void evaluate_constraint_jacobian(Iterate& iterate, std::vector<SparseVector<double>>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers, SymmetricMatrix& hessian) const override;

   [[nodiscard]] double predicted_reduction_contribution(const Iterate& current_iterate, const Direction& direction, double step_length) const override;

   [[nodiscard]] double get_variable_lower_bound(size_t i) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t i) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t j) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t j) const override;

   [[nodiscard]] size_t get_maximum_number_hessian_nonzeros() const override;
};

inline OptimalityProblem::OptimalityProblem(const Model& model):
      NonlinearProblem(model, model.number_variables, model.number_constraints) {
   // figure out bounded variables
   for (size_t i: this->model.lower_bounded_variables) {
      this->lower_bounded_variables.push_back(i);
   }
   for (size_t i: this->model.upper_bounded_variables) {
      this->upper_bounded_variables.push_back(i);
   }
}

inline double OptimalityProblem::get_objective_multiplier() const {
   return 1.;
}

// return rho*f(x) + coeff*(e^T p + e^T n) + proximal
inline double OptimalityProblem::evaluate_objective(Iterate& iterate) const {
   iterate.evaluate_objective(this->model);
   return iterate.original_evaluations.objective;
}

inline void OptimalityProblem::evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const {
   iterate.evaluate_objective_gradient(this->model);
   objective_gradient = iterate.original_evaluations.objective_gradient;
}

inline void OptimalityProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
   iterate.evaluate_constraints(this->model);
   copy_from(constraints, iterate.original_evaluations.constraints);
}

inline void OptimalityProblem::evaluate_constraint_jacobian(Iterate& iterate, std::vector<SparseVector<double>>& constraint_jacobian) const {
   iterate.evaluate_constraint_jacobian(this->model);
   constraint_jacobian = iterate.original_evaluations.constraint_jacobian;
}

inline void OptimalityProblem::evaluate_lagrangian_hessian(const std::vector<double>& x, const std::vector<double>& multipliers,
      SymmetricMatrix& hessian) const {
   const double objective_multiplier = 1.;
   this->model.evaluate_lagrangian_hessian(x, objective_multiplier, multipliers, hessian);
   hessian.dimension = this->number_variables;
}

inline double OptimalityProblem::predicted_reduction_contribution(const Iterate& /*current_iterate*/, const Direction& /*direction*/,
      double /*step_length*/) const {
   return 0.;
}

inline double OptimalityProblem::get_variable_lower_bound(size_t i) const {
   return this->model.get_variable_lower_bound(i);
}

inline double OptimalityProblem::get_variable_upper_bound(size_t i) const {
   return this->model.get_variable_upper_bound(i);
}

inline double OptimalityProblem::get_constraint_lower_bound(size_t j) const {
   return this->model.get_constraint_lower_bound(j);
}

inline double OptimalityProblem::get_constraint_upper_bound(size_t j) const {
   return this->model.get_constraint_upper_bound(j);
}

inline size_t OptimalityProblem::get_maximum_number_hessian_nonzeros() const {
   return this->model.get_maximum_number_hessian_nonzeros();
}

#endif // UNO_OPTIMALITYPROBLEM_H