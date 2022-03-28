#ifndef UNO_SCALEDREFORMULATION_H
#define UNO_SCALEDREFORMULATION_H

#include "Problem.hpp"
#include "Scaling.hpp"

class ScaledReformulation: public Problem {
public:
   ScaledReformulation(const Problem& original_problem, const Scaling& scaling);

   [[nodiscard]] size_t get_number_original_variables() const override;
   [[nodiscard]] double get_variable_lower_bound(size_t i) const override;
   [[nodiscard]] double get_variable_upper_bound(size_t i) const override;
   [[nodiscard]] double get_constraint_lower_bound(size_t j) const override;
   [[nodiscard]] double get_constraint_upper_bound(size_t j) const override;
   [[nodiscard]] double evaluate_objective(Iterate& iterate) const override;
   void evaluate_objective_gradient(Iterate& iterate) const override;
   void evaluate_constraints(Iterate& iterate) const override;
   void evaluate_constraint_jacobian(Iterate& iterate) const override;
   void evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier, const std::vector<double>& multipliers,
         SymmetricMatrix& hessian) const override;

   [[nodiscard]] ConstraintType get_variable_status(size_t i) const override;
   [[nodiscard]] FunctionType get_constraint_type(size_t j) const override;
   [[nodiscard]] ConstraintType get_constraint_status(size_t j) const override;
   [[nodiscard]] size_t get_hessian_maximum_number_nonzeros() const override;

   void get_initial_primal_point(std::vector<double>& x) const override;
   void get_initial_dual_point(std::vector<double>& multipliers) const override;

private:
   const Problem& original_problem;
   const Scaling& scaling;
};

inline ScaledReformulation::ScaledReformulation(const Problem& original_problem, const Scaling& scaling):
      Problem(original_problem.name + "_scaled", // name
            original_problem.number_variables, // number of variables
            original_problem.number_constraints, // number of constraints
            original_problem.problem_type), // problem type
      original_problem(original_problem),
      scaling(scaling) {
   // check the scaling factors
   assert(0 <= this->scaling.get_objective_scaling() && "Objective scaling failed.");
   for (size_t j = 0; j < this->number_constraints; j++) {
      assert(0 <= this->scaling.get_constraint_scaling(j) && "Constraint scaling failed.");
   }

   // the constraint repartition (inequality/equality, linear) is the same as in the original problem
   this->original_problem.equality_constraints.for_each([&](size_t j, size_t i) {
      this->equality_constraints.insert(j, i);
   });
   this->original_problem.inequality_constraints.for_each([&](size_t j, size_t i) {
      this->inequality_constraints.insert(j, i);
   });
   this->original_problem.linear_constraints.for_each([&](size_t j, size_t i) {
      this->linear_constraints.insert(j, i);
   });

   // figure out bounded variables
   for (size_t i: this->original_problem.lower_bounded_variables) {
      this->lower_bounded_variables.push_back(i);
   }
   for (size_t i: this->original_problem.upper_bounded_variables) {
      this->upper_bounded_variables.push_back(i);
   }
}

inline size_t ScaledReformulation::get_number_original_variables() const {
   return this->original_problem.get_number_original_variables();
}

inline double ScaledReformulation::get_variable_lower_bound(size_t i) const {
   return this->original_problem.get_variable_lower_bound(i);
}

inline double ScaledReformulation::get_variable_upper_bound(size_t i) const {
   return this->original_problem.get_variable_upper_bound(i);
}

inline double ScaledReformulation::get_constraint_lower_bound(size_t j) const {
   const double lb = this->original_problem.get_constraint_lower_bound(j);
   // scale
   return this->scaling.get_constraint_scaling(j)*lb;
}

inline double ScaledReformulation::get_constraint_upper_bound(size_t j) const {
   const double ub = this->original_problem.get_constraint_upper_bound(j);
   // scale
   return this->scaling.get_constraint_scaling(j)*ub;
}

inline double ScaledReformulation::evaluate_objective(Iterate& iterate) const {
   const double objective = this->original_problem.evaluate_objective(iterate);
   // scale
   return this->scaling.get_objective_scaling()*objective;
}

inline void ScaledReformulation::evaluate_objective_gradient(Iterate& iterate) const {
   this->original_problem.evaluate_objective_gradient(iterate);
   // scale
   scale(iterate.problem_evaluations.objective_gradient, this->scaling.get_objective_scaling());
}

inline void ScaledReformulation::evaluate_constraints(Iterate& iterate) const {
   this->original_problem.evaluate_constraints(iterate);
   // scale
   for (size_t j = 0; j < this->number_constraints; j++) {
      iterate.problem_evaluations.constraints[j] *= this->scaling.get_constraint_scaling(j);
   }
}

inline void ScaledReformulation::evaluate_constraint_jacobian(Iterate& iterate) const {
   // evaluate
   this->original_problem.evaluate_constraint_jacobian(iterate);
   // scale
   for (size_t j = 0; j < this->number_constraints; j++) {
      scale(iterate.problem_evaluations.constraint_jacobian[j], this->scaling.get_constraint_scaling(j));
   }
}

inline void ScaledReformulation::evaluate_lagrangian_hessian(const std::vector<double>& x, double objective_multiplier,
      const std::vector<double>& multipliers, SymmetricMatrix& hessian) const {
   // scale the objective and constraint multipliers
   const double scaled_objective_multiplier = objective_multiplier*this->scaling.get_objective_scaling();
   // TODO preallocate this vector
   static std::vector<double> scaled_multipliers(this->number_constraints);
   for (size_t j = 0; j < this->number_constraints; j++) {
      scaled_multipliers[j] = scaling.get_constraint_scaling(j)*multipliers[j];
   }
   this->original_problem.evaluate_lagrangian_hessian(x, scaled_objective_multiplier, scaled_multipliers, hessian);
}

inline ConstraintType ScaledReformulation::get_variable_status(size_t i) const {
   return this->original_problem.get_variable_status(i);
}

inline FunctionType ScaledReformulation::get_constraint_type(size_t j) const {
   return this->original_problem.get_constraint_type(j);
}

inline ConstraintType ScaledReformulation::get_constraint_status(size_t j) const {
   return this->original_problem.get_constraint_status(j);
}

inline size_t ScaledReformulation::get_hessian_maximum_number_nonzeros() const {
   return this->original_problem.get_hessian_maximum_number_nonzeros();
}

inline void ScaledReformulation::get_initial_primal_point(std::vector<double>& x) const {
   this->original_problem.get_initial_primal_point(x);
}

inline void ScaledReformulation::get_initial_dual_point(std::vector<double>& multipliers) const {
   this->original_problem.get_initial_dual_point(multipliers);
}

#endif // UNO_SCALEDREFORMULATION_H