// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "PrimalExponentialBarrierProblem.hpp"
#include "ingredients/hessian_models/HessianModel.hpp"
#include "linear_algebra/Indexing.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Iterate.hpp"
#include "symbolic/UnaryNegation.hpp"
#include "symbolic/VectorView.hpp"
#include "tools/Infinity.hpp"
#include "tools/Logger.hpp"

namespace uno {
   PrimalExponentialBarrierProblem::PrimalExponentialBarrierProblem(const OptimizationProblem& problem, double barrier_parameter,
      const InteriorPointParameters& parameters):
         BarrierProblem(problem.model, problem.number_variables, 0),
         problem(problem), barrier_parameter(barrier_parameter), parameters(parameters) {
   }

   void PrimalExponentialBarrierProblem::generate_initial_iterate(Iterate& initial_iterate) const {
      // set the bound multipliers
      for (size_t constraint_index: Range(this->problem.number_constraints)) {
         initial_iterate.multipliers.constraints[constraint_index] = this->parameters.default_multiplier;
      }
      for (const size_t variable_index: Range(this->problem.number_variables)) {
         initial_iterate.multipliers.lower_bounds[variable_index] = this->parameters.default_multiplier;
         initial_iterate.multipliers.upper_bounds[variable_index] = this->parameters.default_multiplier;
      }
   }

   double PrimalExponentialBarrierProblem::get_objective_multiplier() const {
      return this->problem.get_objective_multiplier();
   }

   void PrimalExponentialBarrierProblem::evaluate_constraints(Iterate& /*iterate*/, std::vector<double>& /*constraints*/) const {
      // no constraints
   }

   void PrimalExponentialBarrierProblem::evaluate_objective_gradient(Iterate& iterate, const EvaluationSpace& evaluation_space,
         double* objective_gradient) const {
      // original objective gradient
      this->problem.evaluate_objective_gradient(iterate, evaluation_space, objective_gradient);

      // constraints
      std::vector<double> constraints(this->problem.number_constraints);
      this->problem.evaluate_constraints(iterate, constraints);

      // multipliers
      std::vector<double> constraint_multipliers(this->problem.number_constraints);
      for (size_t constraint_index: Range(this->problem.number_constraints)) {
         const double lower_bound = this->problem.constraint_lower_bound(constraint_index);
         const double upper_bound = this->problem.constraint_upper_bound(constraint_index);
         constraint_multipliers[constraint_index] = 0.;
         if (is_finite(lower_bound)) {
            constraint_multipliers[constraint_index] += std::exp(-(constraints[constraint_index] - lower_bound)/this->barrier_parameter);
         }
         if (is_finite(upper_bound)) {
            constraint_multipliers[constraint_index] -= std::exp(-(upper_bound - constraints[constraint_index])/this->barrier_parameter);
         }
      }
      std::vector<double> lower_bound_multipliers(this->problem.number_variables);
      std::vector<double> upper_bound_multipliers(this->problem.number_variables);
      for (size_t variable_index: Range(this->problem.number_variables)) {
         const double lower_bound = this->problem.variable_lower_bound(variable_index);
         const double upper_bound = this->problem.variable_upper_bound(variable_index);
         lower_bound_multipliers[variable_index] = 0.;
         upper_bound_multipliers[variable_index] = 0.;
         if (is_finite(lower_bound)) {
            lower_bound_multipliers[variable_index] = std::exp(-(iterate.primals[variable_index] - lower_bound)/this->barrier_parameter);
         }
         if (is_finite(upper_bound)) {
            upper_bound_multipliers[variable_index] = -std::exp(-(upper_bound - iterate.primals[variable_index])/this->barrier_parameter);
         }
      }
      std::cout << "Constraint multipliers: "; print_vector(std::cout, constraint_multipliers);
      std::cout << "LB multipliers: "; print_vector(std::cout, lower_bound_multipliers);
      std::cout << "UB multipliers: "; print_vector(std::cout, upper_bound_multipliers);

      // Jacobian
      const size_t number_jacobian_nonzeros = this->problem.number_jacobian_nonzeros();
      std::vector<int> jacobian_row_indices(number_jacobian_nonzeros);
      std::vector<int> jacobian_column_indices(number_jacobian_nonzeros);
      this->problem.compute_constraint_jacobian_sparsity(jacobian_row_indices.data(), jacobian_column_indices.data(),
         Indexing::C_indexing, MatrixOrder::COLUMN_MAJOR);
      std::vector<double> jacobian_values(number_jacobian_nonzeros);
      this->problem.evaluate_constraint_jacobian(iterate, jacobian_values.data());

      for (size_t nonzero_index: Range(number_jacobian_nonzeros)) {
         const size_t constraint_index = static_cast<size_t>(jacobian_row_indices[nonzero_index]);
         const size_t variable_index = static_cast<size_t>(jacobian_column_indices[nonzero_index]);
         const double derivative = jacobian_values[nonzero_index];

         objective_gradient[variable_index] -= derivative * constraint_multipliers[constraint_index];
      }

      // bound constraints

      std::cout << "PrimalExponentialBarrierProblem::evaluate_objective_gradient\n";
      for (size_t index: Range(this->problem.number_variables)) {
         std::cout << objective_gradient[index] << ' ';
      }
      std::cout << '\n';
   }

   size_t PrimalExponentialBarrierProblem::number_jacobian_nonzeros() const {
      return 0;
   }

   bool PrimalExponentialBarrierProblem::has_curvature(const HessianModel& hessian_model) const {
      return hessian_model.has_curvature(this->model);
   }

   size_t PrimalExponentialBarrierProblem::number_hessian_nonzeros(const HessianModel& hessian_model) const {
      return this->problem.number_hessian_nonzeros(hessian_model);
   }

   void PrimalExponentialBarrierProblem::compute_constraint_jacobian_sparsity(int* /*row_indices*/, int* /*column_indices*/,
         int /*solver_indexing*/, MatrixOrder /*matrix_order*/) const {
   }

   void PrimalExponentialBarrierProblem::compute_hessian_sparsity(const HessianModel& hessian_model, int* row_indices,
         int* column_indices, int solver_indexing) const {
      this->problem.compute_hessian_sparsity(hessian_model, row_indices, column_indices, solver_indexing);
   }

   void PrimalExponentialBarrierProblem::evaluate_constraint_jacobian(Iterate& /*iterate*/, double* /*jacobian_values*/) const {
   }

   void PrimalExponentialBarrierProblem::evaluate_lagrangian_gradient(LagrangianGradient<double>& lagrangian_gradient,
         const InequalityHandlingMethod& inequality_handling_method, const EvaluationSpace& evaluation_space, Iterate& iterate) const {
      this->problem.evaluate_lagrangian_gradient(lagrangian_gradient, inequality_handling_method,
         evaluation_space, iterate);
   }
   
   void PrimalExponentialBarrierProblem::evaluate_lagrangian_hessian(Statistics& statistics, HessianModel& hessian_model,
         const Vector<double>& primal_variables, const Multipliers& multipliers, double* hessian_values) const {
      // original Lagrangian Hessian
      this->problem.evaluate_lagrangian_hessian(statistics, hessian_model, primal_variables, multipliers, hessian_values);
   }

   void PrimalExponentialBarrierProblem::compute_hessian_vector_product(HessianModel& hessian_model, const double* vector,
         const Multipliers& multipliers, double* result) const {
      // original Lagrangian Hessian
      this->problem.compute_hessian_vector_product(hessian_model, vector, multipliers, result);
   }

   double PrimalExponentialBarrierProblem::variable_lower_bound(size_t /*variable_index*/) const {
      return -INF<double>;
   }

   double PrimalExponentialBarrierProblem::variable_upper_bound(size_t /*variable_index*/) const {
      return INF<double>;
   }

   const Vector<size_t>& PrimalExponentialBarrierProblem::get_fixed_variables() const {
      return this->fixed_variables;
   }

   const Collection<size_t>& PrimalExponentialBarrierProblem::get_primal_regularization_variables() const {
      return this->problem.get_primal_regularization_variables();
   }

   double PrimalExponentialBarrierProblem::constraint_lower_bound(size_t constraint_index) const {
      return this->problem.constraint_lower_bound(constraint_index);
   }

   double PrimalExponentialBarrierProblem::constraint_upper_bound(size_t constraint_index) const {
      return this->problem.constraint_upper_bound(constraint_index);
   }

   const Collection<size_t>& PrimalExponentialBarrierProblem::get_equality_constraints() const {
      return this->problem.get_equality_constraints();
   }

   const Collection<size_t>& PrimalExponentialBarrierProblem::get_inequality_constraints() const {
      return this->empty_set;
   }

   const Collection<size_t>& PrimalExponentialBarrierProblem::get_dual_regularization_constraints() const {
      if (this->problem.get_dual_regularization_constraints().empty()) {
         // this is an indication that the constraints (if there is any) were already regularized in a previous
         // reformulation (e.g. l1 relaxation). In that case, we stick to an empty set
         return this->problem.get_dual_regularization_constraints();
      }
      // otherwise, we pick the set of equality constraints
      return this->problem.get_equality_constraints();
   }

   void PrimalExponentialBarrierProblem::assemble_primal_dual_direction(const Iterate& /*current_iterate*/, const Vector<double>& solution,
         Direction& direction) const {
      // form the primal-dual direction
      direction.primals = view(solution, 0, this->problem.number_variables);
      // retrieve the duals with correct signs (note the minus sign)
      direction.multipliers.constraints = view(-solution, this->problem.number_variables,
         this->problem.number_variables + this->problem.number_constraints);
   }

   double PrimalExponentialBarrierProblem::push_variable_to_interior(double variable_value, double /*lower_bound*/, double /*upper_bound*/) const {
      return variable_value;
   }

   void PrimalExponentialBarrierProblem::set_auxiliary_measure(Iterate& iterate) const {
      // auxiliary measure: barrier terms
      double barrier_terms = 0.;
      iterate.progress.auxiliary = barrier_terms;
   }

   double PrimalExponentialBarrierProblem::complementarity_error(const Vector<double>& primals, const std::vector<double>& constraints,
         const Multipliers& multipliers, double shift_value, Norm residual_norm) const {
      return this->problem.complementarity_error(primals, constraints, multipliers, shift_value, residual_norm);
   }

   double PrimalExponentialBarrierProblem::dual_regularization_factor() const {
      return std::pow(this->barrier_parameter, this->parameters.dual_regularization_exponent);
   }

   // protected member functions

   double PrimalExponentialBarrierProblem::compute_barrier_term_directional_derivative(const Iterate& /*current_iterate*/,
         const Vector<double>& /*primal_direction*/) const {
      double directional_derivative = 0.;
      return directional_derivative;
   }

   void PrimalExponentialBarrierProblem::postprocess_iterate(Iterate& /*iterate*/) const {
      // do nothing
   }

   double PrimalExponentialBarrierProblem::compute_centrality_error(const Vector<double>& primals,
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
} // namespace