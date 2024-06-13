// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_L1RELAXEDPROBLEM_H
#define UNO_L1RELAXEDPROBLEM_H

#include "OptimizationProblem.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "model/Model.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/LagrangianGradient.hpp"
#include "symbolic/Expression.hpp"
#include "symbolic/VectorExpression.hpp"
#include "symbolic/Concatenation.hpp"
#include "symbolic/Range.hpp"
#include "tools/Infinity.hpp"

struct ElasticVariables {
   SparseVector<size_t> positive;
   SparseVector<size_t> negative;
   explicit ElasticVariables(size_t capacity): positive(capacity), negative(capacity) {}
   [[nodiscard]] size_t size() const { return this->positive.size() + this->negative.size(); }
};

class l1RelaxedProblem: public OptimizationProblem {
public:
   l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient);

   [[nodiscard]] double get_objective_multiplier() const override;
   void evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const override;
   void evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const override;
   void evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const override;
   void evaluate_lagrangian_hessian(const Vector<double>& x, const Vector<double>& multipliers, SymmetricMatrix<double>& hessian) const override;

   [[nodiscard]] double variable_lower_bound(size_t variable_index) const override;
   [[nodiscard]] double variable_upper_bound(size_t variable_index) const override;
   [[nodiscard]] const Collection<size_t>& get_lower_bounded_variables() const override;
   [[nodiscard]] const Collection<size_t>& get_upper_bounded_variables() const override;
   [[nodiscard]] const Collection<size_t>& get_single_lower_bounded_variables() const override;
   [[nodiscard]] const Collection<size_t>& get_single_upper_bounded_variables() const override;

   [[nodiscard]] double constraint_lower_bound(size_t constraint_index) const override;
   [[nodiscard]] double constraint_upper_bound(size_t constraint_index) const override;

   [[nodiscard]] size_t number_objective_gradient_nonzeros() const override;
   [[nodiscard]] size_t number_jacobian_nonzeros() const override;
   [[nodiscard]] size_t number_hessian_nonzeros() const override;

   [[nodiscard]] double complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, Norm residual_norm) const override;

   // parameterization
   void set_objective_multiplier(double new_objective_multiplier);

   void set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t, double)>& elastic_setting_function) const;

protected:
   double objective_multiplier;
   const double constraint_violation_coefficient;
   ElasticVariables elastic_variables;
   const Concatenation<const Collection<size_t>&, ForwardRange> lower_bounded_variables; // model variables + elastic variables
   const Concatenation<const Collection<size_t>&, ForwardRange> single_lower_bounded_variables; // model variables + elastic variables

   [[nodiscard]] static size_t count_elastic_variables(const Model& model);
   void generate_elastic_variables();
};

inline l1RelaxedProblem::l1RelaxedProblem(const Model& model, double objective_multiplier, double constraint_violation_coefficient):
      OptimizationProblem(model, model.number_variables + l1RelaxedProblem::count_elastic_variables(model), model.number_constraints),
      objective_multiplier(objective_multiplier),
      constraint_violation_coefficient(constraint_violation_coefficient),
      elastic_variables(this->number_constraints),
      // lower bounded variables are the model variables + the elastic variables
      lower_bounded_variables(concatenate(this->model.get_lower_bounded_variables(), Range(model.number_variables,
            model.number_variables + this->elastic_variables.size()))),
      single_lower_bounded_variables(concatenate(this->model.get_single_lower_bounded_variables(),
            Range(model.number_variables, model.number_variables + this->elastic_variables.size()))) {
   this->generate_elastic_variables();
}

inline double l1RelaxedProblem::get_objective_multiplier() const {
   return this->objective_multiplier;
}

inline void l1RelaxedProblem::evaluate_objective_gradient(Iterate& iterate, SparseVector<double>& objective_gradient) const {
   // scale nabla f(x) by rho
   if (this->objective_multiplier != 0.) {
      iterate.evaluate_objective_gradient(this->model);
      // TODO change this
      objective_gradient = iterate.evaluations.objective_gradient;
      scale(objective_gradient, this->objective_multiplier);
   }
   else {
      objective_gradient.clear();
   }

   // elastic contribution
   for (const auto [_, elastic_index]: this->elastic_variables.positive) {
      objective_gradient.insert(elastic_index, this->constraint_violation_coefficient);
   }
   for (const auto [_, elastic_index]: this->elastic_variables.negative) {
      objective_gradient.insert(elastic_index, this->constraint_violation_coefficient);
   }
}

inline void l1RelaxedProblem::evaluate_constraints(Iterate& iterate, std::vector<double>& constraints) const {
   iterate.evaluate_constraints(this->model);
   constraints = iterate.evaluations.constraints;
   // add the contribution of the elastics
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
      constraints[constraint_index] -= iterate.primals[elastic_index];
   }
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
      constraints[constraint_index] += iterate.primals[elastic_index];
   }
}

inline void l1RelaxedProblem::evaluate_constraint_jacobian(Iterate& iterate, RectangularMatrix<double>& constraint_jacobian) const {
   iterate.evaluate_constraint_jacobian(this->model);
   // TODO change this
   constraint_jacobian = iterate.evaluations.constraint_jacobian;
   // add the contribution of the elastics
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
      constraint_jacobian[constraint_index].insert(elastic_index, -1.);
   }
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
      constraint_jacobian[constraint_index].insert(elastic_index, 1.);
   }
}

inline void l1RelaxedProblem::evaluate_lagrangian_hessian(const Vector<double>& x, const Vector<double>& multipliers,
      SymmetricMatrix<double>& hessian) const {
   this->model.evaluate_lagrangian_hessian(x, this->objective_multiplier, multipliers, hessian);

   // extend the dimension of the Hessian by finalizing the remaining columns (note: the elastics do not enter the Hessian)
   for (size_t constraint_index: Range(this->model.number_variables, this->number_variables)) {
      hessian.finalize_column(constraint_index);
   }
}

// complementary slackness error: expression for violated constraints depends on the definition of the relaxed problem
inline double l1RelaxedProblem::complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
      const Multipliers& multipliers, Norm residual_norm) const {
   // complementarity for variable bounds
   const VectorExpression variable_complementarity(Range(this->model.number_variables), [&](size_t variable_index) {
      if (0. < multipliers.lower_bounds[variable_index]) {
         return multipliers.lower_bounds[variable_index] * (primals[variable_index] - this->variable_lower_bound(variable_index));
      }
      if (multipliers.upper_bounds[variable_index] < 0.) {
         return multipliers.upper_bounds[variable_index] * (primals[variable_index] - this->variable_upper_bound(variable_index));
      }
      return 0.;
   });

   // complementarity for constraint bounds
   const VectorExpression constraint_complementarity(Range(constraints.size()), [&](size_t constraint_index) {
      // violated constraints
      if (constraints[constraint_index] < this->constraint_lower_bound(constraint_index)) { // lower violated
         return (this->constraint_violation_coefficient - multipliers.constraints[constraint_index]) * (constraints[constraint_index] -
               this->constraint_lower_bound(constraint_index));
      }
      else if (this->constraint_upper_bound(constraint_index) < constraints[constraint_index]) { // upper violated
         return (this->constraint_violation_coefficient + multipliers.constraints[constraint_index]) * (constraints[constraint_index] -
               this->constraint_upper_bound(constraint_index));
      }
      // satisfied constraints
      else if (0. < multipliers.constraints[constraint_index]) { // lower bound
         return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->constraint_lower_bound(constraint_index));
      }
      else if (multipliers.constraints[constraint_index] < 0.) { // upper bound
         return multipliers.constraints[constraint_index] * (constraints[constraint_index] - this->constraint_upper_bound(constraint_index));
      }
      return 0.;
   });
   return norm(residual_norm, variable_complementarity, constraint_complementarity);
}

inline double l1RelaxedProblem::variable_lower_bound(size_t variable_index) const {
   if (variable_index < this->model.number_variables) { // model variable
      return this->model.variable_lower_bound(variable_index);
   }
   else { // elastic variable in [0, +inf[
      return 0.;
   }
}

inline double l1RelaxedProblem::variable_upper_bound(size_t variable_index) const {
   if (variable_index < this->model.number_variables) { // model variable
      return this->model.variable_upper_bound(variable_index);
   }
   else { // elastic variable in [0, +inf[
      return INF<double>;
   }
}

inline double l1RelaxedProblem::constraint_lower_bound(size_t constraint_index) const {
   return this->model.constraint_lower_bound(constraint_index);
}

inline double l1RelaxedProblem::constraint_upper_bound(size_t constraint_index) const {
   return this->model.constraint_upper_bound(constraint_index);
}

inline const Collection<size_t>& l1RelaxedProblem::get_lower_bounded_variables() const {
   return this->lower_bounded_variables;
}

inline const Collection<size_t>& l1RelaxedProblem::get_upper_bounded_variables() const {
   // same set as the model
   return this->model.get_upper_bounded_variables();
}

inline const Collection<size_t>& l1RelaxedProblem::get_single_lower_bounded_variables() const {
   return this->single_lower_bounded_variables;
}

inline const Collection<size_t>& l1RelaxedProblem::get_single_upper_bounded_variables() const {
   // same set as the model
   return this->model.get_single_upper_bounded_variables();
}

inline size_t l1RelaxedProblem::number_objective_gradient_nonzeros() const {
   // elastic contribution
   size_t number_nonzeros = this->elastic_variables.size();

   // objective contribution
   if (this->objective_multiplier != 0.) {
      number_nonzeros += this->model.number_objective_gradient_nonzeros();
   }
   return number_nonzeros;
}

inline size_t l1RelaxedProblem::number_jacobian_nonzeros() const {
   return this->model.number_jacobian_nonzeros() + this->elastic_variables.size();
}

inline size_t l1RelaxedProblem::number_hessian_nonzeros() const {
   return this->model.number_hessian_nonzeros();
}

inline void l1RelaxedProblem::set_objective_multiplier(double new_objective_multiplier) {
   assert(0. <= new_objective_multiplier && "The objective multiplier should be non-negative");
   this->objective_multiplier = new_objective_multiplier;
}

inline size_t l1RelaxedProblem::count_elastic_variables(const Model& model) {
   size_t number_elastic_variables = 0;
   // if the subproblem uses slack variables, the bounds of the constraints are [0, 0]
   for (size_t constraint_index: Range(model.number_constraints)) {
      if (is_finite(model.constraint_lower_bound(constraint_index))) {
         number_elastic_variables++;
      }
      if (is_finite(model.constraint_upper_bound(constraint_index))) {
         number_elastic_variables++;
      }
   }
   return number_elastic_variables;
}

inline void l1RelaxedProblem::generate_elastic_variables() {
   // generate elastic variables to relax the constraints
   size_t elastic_index = this->model.number_variables;
   for (size_t constraint_index: Range(this->model.number_constraints)) {
      if (is_finite(this->model.constraint_upper_bound(constraint_index))) {
         // nonnegative variable p that captures the positive part of the constraint violation
         this->elastic_variables.positive.insert(constraint_index, elastic_index);
         elastic_index++;
      }
      if (is_finite(this->model.constraint_lower_bound(constraint_index))) {
         // nonpositive variable n that captures the negative part of the constraint violation
         this->elastic_variables.negative.insert(constraint_index, elastic_index);
         elastic_index++;
      }
   }
}

inline void l1RelaxedProblem::set_elastic_variable_values(Iterate& iterate, const std::function<void(Iterate&, size_t, size_t, double)>&
      elastic_setting_function) const {
   iterate.set_number_variables(this->number_variables);
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.positive) {
      elastic_setting_function(iterate, constraint_index, elastic_index, -1.);
   }
   for (const auto [constraint_index, elastic_index]: this->elastic_variables.negative) {
      elastic_setting_function(iterate, constraint_index, elastic_index, 1.);
   }
}

#endif // UNO_L1RELAXEDPROBLEM_H
