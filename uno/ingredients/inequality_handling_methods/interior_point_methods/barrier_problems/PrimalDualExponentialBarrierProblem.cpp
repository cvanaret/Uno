// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalDualExponentialBarrierProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   PrimalDualExponentialBarrierProblem::PrimalDualExponentialBarrierProblem(const OptimizationProblem& problem, double barrier_parameter,
      const InteriorPointParameters& parameters):
         // no slacks: as many constraints as the number of equality constraints of the problem
         BarrierProblem(problem.model, problem.number_variables + 2*PrimalDualExponentialBarrierProblem::count_number_extra_variables(problem),
            0),
         problem(problem), number_extra_variables(PrimalDualExponentialBarrierProblem::count_number_extra_variables(problem)),
         barrier_parameter(barrier_parameter), parameters(parameters) {
      DEBUG << "The exponential barrier problem has " << this->number_variables << " variables\n";
   }

   void PrimalDualExponentialBarrierProblem::generate_initial_iterate(Iterate& initial_iterate) const {
      // set the bound multipliers
      for (size_t constraint_index: Range(this->problem.number_constraints)) {
         initial_iterate.multipliers.constraints[constraint_index] = this->parameters.default_multiplier;
      }
      for (const size_t variable_index: Range(this->problem.number_variables)) {
         initial_iterate.multipliers.lower_bounds[variable_index] = this->parameters.default_multiplier;
         initial_iterate.multipliers.upper_bounds[variable_index] = this->parameters.default_multiplier;
      }
   }

   double PrimalDualExponentialBarrierProblem::get_objective_multiplier() const {
      return this->problem.get_objective_multiplier();
   }

   void PrimalDualExponentialBarrierProblem::evaluate_constraints(Iterate& /*iterate*/, std::vector<double>& /*constraints*/) const {
      // no constraints
   }

   void PrimalDualExponentialBarrierProblem::evaluate_objective_gradient(Iterate& iterate, const EvaluationSpace& evaluation_space,
         double* objective_gradient) const {
      // original objective gradient
      this->problem.evaluate_objective_gradient(iterate, evaluation_space, objective_gradient);

      // contributions of inequality constraints
      std::vector<double> original_constraints(this->problem.number_constraints);
      this->problem.evaluate_constraints(iterate, original_constraints);

      //const auto& constraints = evaluation_space.get_constraints();
      size_t gradient_index = this->problem.number_variables;
      // TODO we need to maintain two sets of constraint multipliers
      for (size_t constraint_index: Range(this->problem.number_constraints)) {
         const double lower_bound = this->problem.constraint_lower_bound(constraint_index);
         const double upper_bound = this->problem.constraint_upper_bound(constraint_index);
         if (is_finite(lower_bound)) {
            objective_gradient[gradient_index] = lower_bound - original_constraints[constraint_index] - this->barrier_parameter*
               std::log(iterate.multipliers.constraints[constraint_index]);
            ++gradient_index;
         }
         if (is_finite(upper_bound)) {
            objective_gradient[gradient_index] = original_constraints[constraint_index] - upper_bound - this->barrier_parameter*
               std::log(iterate.multipliers.constraints[constraint_index]);
            ++gradient_index;
         }
      }
      for (size_t variable_index: Range(this->problem.number_variables)) {
         const double lower_bound = this->problem.variable_lower_bound(variable_index);
         const double upper_bound = this->problem.variable_upper_bound(variable_index);
         if (is_finite(lower_bound)) {
            objective_gradient[gradient_index] = lower_bound - iterate.primals[variable_index] - this->barrier_parameter*
               std::log(iterate.multipliers.lower_bounds[variable_index]);
            ++gradient_index;
         }
         if (is_finite(upper_bound)) {
            objective_gradient[gradient_index] = iterate.primals[variable_index] - upper_bound - this->barrier_parameter*
               std::log(iterate.multipliers.upper_bounds[variable_index]);
            ++gradient_index;
         }
      }
      std::cout << "PrimalDualExponentialBarrierProblem::evaluate_objective_gradient\n";
      for (size_t index: Range(gradient_index)) {
         std::cout << objective_gradient[index] << ' ';
      }
      std::cout << '\n';
   }

   size_t PrimalDualExponentialBarrierProblem::number_jacobian_nonzeros() const {
      return 0;
   }

   bool PrimalDualExponentialBarrierProblem::has_curvature(const HessianModel& hessian_model) const {
      return hessian_model.has_curvature(this->model);
   }

   size_t PrimalDualExponentialBarrierProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      return this->problem.number_hessian_nonzeros(hessian_model) + 2*this->number_extra_variables;
   }

   void PrimalDualExponentialBarrierProblem::compute_constraint_jacobian_sparsity(int* /*row_indices*/, int* /*column_indices*/,
         int /*solver_indexing*/, MatrixOrder /*matrix_order*/) const {
   }

   void PrimalDualExponentialBarrierProblem::compute_hessian_sparsity(const HessianModel& hessian_model, int* row_indices,
         int* column_indices, int solver_indexing) const {
      // original Lagrangian Hessian
      this->problem.compute_hessian_sparsity(hessian_model, row_indices, column_indices, solver_indexing);

      // two copies of the Jacobian
      const size_t number_hessian_nonzeros = this->problem.number_hessian_nonzeros(hessian_model);
      // sparsity of the Jacobian
      this->problem.compute_constraint_jacobian_sparsity(row_indices + number_hessian_nonzeros, column_indices +
         number_hessian_nonzeros, solver_indexing, MatrixOrder::COLUMN_MAJOR); // TODO
      const size_t number_jacobian_nonzeros = this->problem.number_jacobian_nonzeros();
      // copy it a second time
      for (size_t index: Range(number_jacobian_nonzeros)) {
         row_indices[number_hessian_nonzeros + number_jacobian_nonzeros + index] = row_indices[number_hessian_nonzeros + index];
         column_indices[number_hessian_nonzeros + number_jacobian_nonzeros + index] = column_indices[number_hessian_nonzeros + index];
      }
      // shift the row indices
      for (size_t index: Range(number_jacobian_nonzeros)) {
         row_indices[number_hessian_nonzeros + index] += static_cast<int>(this->problem.number_variables);
         row_indices[number_hessian_nonzeros + number_jacobian_nonzeros + index] += static_cast<int>(this->problem.number_variables +
            this->problem.number_constraints);
      }
      // diagonal contributions
      row_indices[11] = 7;
      column_indices[11] = 1;
      row_indices[12] = 8;
      column_indices[12] = 1;
   }

   void PrimalDualExponentialBarrierProblem::evaluate_constraint_jacobian(Iterate& /*iterate*/, double* /*jacobian_values*/) const {
   }

   void PrimalDualExponentialBarrierProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient,
         const InequalityHandlingMethod& inequality_handling_method, const EvaluationSpace& evaluation_space, Iterate& iterate) const {
      this->problem.evaluate_lagrangian_gradient(lagrangian_gradient, inequality_handling_method,
         evaluation_space, iterate);
   }
   
   void PrimalDualExponentialBarrierProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model,
         const Vector<double>& primal_variables, const Multipliers& multipliers, double* hessian_values) const {
      // original Lagrangian Hessian
      this->problem.evaluate_lagrangian_hessian(statistics, hessian_model, primal_variables, multipliers, hessian_values);
   }

   void PrimalDualExponentialBarrierProblem::compute_hessian_vector_product(HessianModel& hessian_model, const double* x,
         const double* vector, const Multipliers& multipliers, double* result) const {
      // original Lagrangian Hessian
      this->problem.compute_hessian_vector_product(hessian_model, x, vector, multipliers, result);
   }

   double PrimalDualExponentialBarrierProblem::variable_lower_bound(size_t /*variable_index*/) const {
      return -INF<double>;
   }

   double PrimalDualExponentialBarrierProblem::variable_upper_bound(size_t /*variable_index*/) const {
      return INF<double>;
   }

   const Vector<size_t>& PrimalDualExponentialBarrierProblem::get_fixed_variables() const {
      return this->fixed_variables;
   }

   const Collection<size_t>& PrimalDualExponentialBarrierProblem::get_primal_regularization_variables() const {
      return this->problem.get_primal_regularization_variables();
   }

   double PrimalDualExponentialBarrierProblem::constraint_lower_bound(size_t constraint_index) const {
      return this->problem.constraint_lower_bound(constraint_index);
   }

   double PrimalDualExponentialBarrierProblem::constraint_upper_bound(size_t constraint_index) const {
      return this->problem.constraint_upper_bound(constraint_index);
   }

   const Collection<size_t>& PrimalDualExponentialBarrierProblem::get_equality_constraints() const {
      return this->problem.get_equality_constraints();
   }

   const Collection<size_t>& PrimalDualExponentialBarrierProblem::get_inequality_constraints() const {
      return this->empty_set;
   }

   const Collection<size_t>& PrimalDualExponentialBarrierProblem::get_dual_regularization_constraints() const {
      if (this->problem.get_dual_regularization_constraints().empty()) {
         // this is an indication that the constraints (if there is any) were already regularized in a previous
         // reformulation (e.g. l1 relaxation). In that case, we stick to an empty set
         return this->problem.get_dual_regularization_constraints();
      }
      // otherwise, we pick the set of equality constraints
      return this->problem.get_equality_constraints();
   }

   void PrimalDualExponentialBarrierProblem::assemble_primal_dual_direction(const Iterate& /*current_iterate*/, const Vector<double>& solution,
         Direction& direction) const {
      // form the primal-dual direction
      direction.primals = view(solution, 0, this->problem.number_variables);
      // retrieve the duals with correct signs (note the minus sign)
      direction.multipliers.constraints = view(-solution, this->problem.number_variables,
         this->problem.number_variables + this->problem.number_constraints);
   }

   double PrimalDualExponentialBarrierProblem::push_variable_to_interior(double variable_value, double /*lower_bound*/, double /*upper_bound*/) const {
      return variable_value;
   }

   void PrimalDualExponentialBarrierProblem::set_auxiliary_measure(Iterate& iterate) const {
      // auxiliary measure: barrier terms
      double barrier_terms = 0.;
      iterate.progress.auxiliary = barrier_terms;
   }

   double PrimalDualExponentialBarrierProblem::complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const {
      return this->problem.complementarity_error(primals, constraints, multipliers, shift_value, residual_norm);
   }

   double PrimalDualExponentialBarrierProblem::dual_regularization_factor() const {
      return std::pow(this->barrier_parameter, this->parameters.dual_regularization_exponent);
   }

   // protected member functions

   double PrimalDualExponentialBarrierProblem::compute_barrier_term_directional_derivative(const Iterate& /*current_iterate*/,
         const Vector<double>& /*primal_direction*/) const {
      double directional_derivative = 0.;
      return directional_derivative;
   }

   void PrimalDualExponentialBarrierProblem::postprocess_iterate(Iterate& /*iterate*/) const {
      // do nothing
   }

   double PrimalDualExponentialBarrierProblem::compute_centrality_error(const Vector<double>& primals,
         const Multipliers& multipliers, double shift) const {
      const Range variables_range = Range(this->problem.number_variables);
      const VectorExpression shifted_bound_complementarity{variables_range, [&](size_t variable_index) {
         double result = 0.;
         if (0. < multipliers.lower_bounds[variable_index]) { // lower bound
            result = std::max(result, std::abs(multipliers.lower_bounds[variable_index] *
               (primals[variable_index] - this->problem.variable_lower_bound(variable_index)) - shift));
         }
         if (multipliers.upper_bounds[variable_index] < 0.) { // upper bound
            result = std::max(result, std::abs(multipliers.upper_bounds[variable_index] *
               (primals[variable_index] - this->problem.variable_upper_bound(variable_index)) - shift));
         }
         return result;
      }};
      return norm_inf(shifted_bound_complementarity); // TODO use a generic norm
   }

   size_t PrimalDualExponentialBarrierProblem::count_number_extra_variables(const OptimizationProblem& problem) {
      size_t number_variables = 0;

      // count the number of variables associated to constraint bounds
      for (size_t constraint_index: Range(problem.number_constraints)) {
         if (is_finite(problem.constraint_lower_bound(constraint_index))) {
            number_variables++;
         }
         if (is_finite(problem.constraint_upper_bound(constraint_index))) {
            number_variables++;
         }
      }
      for (size_t variable_index: Range(problem.number_variables)) {
         if (is_finite(problem.variable_lower_bound(variable_index))) {
            number_variables++;
         }
         if (is_finite(problem.variable_upper_bound(variable_index))) {
            number_variables++;
         }
      }
      return number_variables;
   }
} // namespace