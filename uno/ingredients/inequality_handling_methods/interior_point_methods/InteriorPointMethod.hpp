// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_PRIMALDUALINTERIORPOINTMETHOD_H
#define UNO_PRIMALDUALINTERIORPOINTMETHOD_H

#include <cmath>
#include <memory>
#include "../InequalityHandlingMethod.hpp"
#include "BarrierParameterUpdateStrategy.hpp"
#include "InteriorPointParameters.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/constraint_relaxation_strategies/l1RelaxedProblem.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "linear_algebra/SparseVector.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationSpace.hpp"
#include "optimization/Iterate.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   // forward references
   class DualResiduals;
   class Subproblem;

   template <typename BarrierProblem>
   class InteriorPointMethod : public InequalityHandlingMethod {
   public:
      explicit InteriorPointMethod(const Options& options);

      void initialize(const OptimizationProblem& problem, Iterate& current_iterate,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius) override;
      void initialize_statistics(Statistics& statistics, const Options& options) override;
      void generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) override;
      void set_initial_point(const Vector<double>& point) override;

      void initialize_feasibility_problem(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      void exit_feasibility_problem(const OptimizationProblem& problem, Iterate& trial_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) override;
      [[nodiscard]] double proximal_coefficient() const override;

      void solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         Direction& direction, HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy,
         double trust_region_radius, WarmstartInformation& warmstart_information) override;

      // matrix computations
      [[nodiscard]] EvaluationSpace& get_evaluation_space() const override;
      void evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate) override;
      void compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      void compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const override;
      [[nodiscard]] double compute_hessian_quadratic_product(const Vector<double>& vector) const override;

      void set_auxiliary_measure(const OptimizationProblem& problem, Iterate& iterate) override;
      [[nodiscard]] double compute_predicted_auxiliary_reduction_model(const OptimizationProblem& problem, const Iterate& current_iterate,
         const Vector<double>& primal_direction, double step_length) const override;

      void postprocess_iterate(const OptimizationProblem& problem, Iterate& iterate) override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      const std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
      BarrierParameterUpdateStrategy<BarrierProblem> barrier_parameter_update_strategy;
      double previous_barrier_parameter;
      const double default_multiplier;
      const InteriorPointParameters parameters;
      const double least_square_multiplier_max_norm;
      const double l1_constraint_violation_coefficient; // (rho in Section 3.3.1 in IPOPT paper)

      bool solving_feasibility_problem{false};
      bool first_feasibility_iteration{false};

      [[nodiscard]] double barrier_parameter() const;
      void update_barrier_parameter(const BarrierProblem& barrier_problem, const Iterate& current_iterate,
         const DualResiduals& residuals);
      [[nodiscard]] bool is_small_step(const OptimizationProblem& problem, const Vector<double>& current_primals, const Vector<double>& direction_primals) const;
      [[nodiscard]] double evaluate_subproblem_objective(const Direction& direction) const;
   };

   template <typename BarrierProblem>
   InteriorPointMethod<BarrierProblem>::InteriorPointMethod(const Options& options):
         InequalityHandlingMethod(),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"))),
         barrier_parameter_update_strategy(options),
         previous_barrier_parameter(options.get_double("barrier_initial_parameter")),
         default_multiplier(options.get_double("barrier_default_multiplier")),
         parameters({
               options.get_double("barrier_tau_min"),
               options.get_double("barrier_k_sigma"),
               options.get_double("barrier_regularization_exponent"),
               options.get_double("barrier_small_direction_factor"),
               options.get_double("barrier_push_variable_to_interior_k1"),
               options.get_double("barrier_push_variable_to_interior_k2"),
               options.get_double("barrier_damping_factor")
         }),
         least_square_multiplier_max_norm(options.get_double("least_square_multiplier_max_norm")),
         l1_constraint_violation_coefficient(options.get_double("l1_constraint_violation_coefficient")) {
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::initialize(const OptimizationProblem& problem, Iterate& current_iterate,
         HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy, double trust_region_radius) {
      if (!problem.get_fixed_variables().empty()) {
         throw std::runtime_error("The problem has fixed variables. Move them to the set of general constraints.");
      }
      const BarrierProblem barrier_problem(problem, this->barrier_parameter(), this->parameters);
      const Subproblem subproblem{barrier_problem, current_iterate, hessian_model, regularization_strategy, trust_region_radius};
      this->linear_solver->initialize_augmented_system(subproblem);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::initialize_statistics(Statistics& statistics, const Options& options) {
      statistics.add_column("barrier", Statistics::double_width - 5, options.get_int("statistics_barrier_parameter_column_order"));
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::generate_initial_iterate(const OptimizationProblem& problem, Iterate& initial_iterate) {
      // TODO: enforce linear constraints at initial point
      //if (options.get_bool("enforce_linear_constraints")) {
      //   Preprocessing::enforce_linear_constraints(problem.model, initial_iterate.primals, initial_iterate.multipliers, this->solver);
      //}

      const BarrierProblem barrier_problem(problem, this->barrier_parameter(), this->parameters);

      // add the slacks to the initial iterate
      initial_iterate.set_number_variables(problem.number_variables);
      // make the initial point strictly feasible wrt the bounds
      for (size_t variable_index: Range(problem.number_variables)) {
         initial_iterate.primals[variable_index] = barrier_problem.push_variable_to_interior(initial_iterate.primals[variable_index],
            problem.variable_lower_bound(variable_index), problem.variable_upper_bound(variable_index));
      }

      // set the slack variables (if any)
      if (!problem.model.get_slacks().is_empty()) {
         // set the slacks to the constraint values
         initial_iterate.evaluate_constraints(problem.model);
         for (const auto [constraint_index, slack_index]: problem.model.get_slacks()) {
            initial_iterate.primals[slack_index] =
               barrier_problem.push_variable_to_interior(initial_iterate.evaluations.constraints[constraint_index],
               problem.variable_lower_bound(slack_index), problem.variable_upper_bound(slack_index));
         }
         // since the slacks have been set, the function evaluations should also be updated
         initial_iterate.is_objective_gradient_computed = false;
         initial_iterate.are_constraints_computed = false;
         initial_iterate.is_constraint_jacobian_computed = false;
      }

      // set the bound multipliers
      for (size_t variable_index: Range(problem.number_variables)) {
         const double lower_bound = problem.variable_lower_bound(variable_index);
         const double upper_bound = problem.variable_upper_bound(variable_index);
         if (is_finite(lower_bound)) {
            initial_iterate.multipliers.lower_bounds[variable_index] = this->default_multiplier;
         }
         if (is_finite(upper_bound)) {
            initial_iterate.multipliers.upper_bounds[variable_index] = -this->default_multiplier;
         }
      }

      // compute least-square multipliers
      if (0 < problem.number_constraints) {
         // TODO
      }
   }

   template <typename BarrierProblem>
   double InteriorPointMethod<BarrierProblem>::barrier_parameter() const {
      return this->barrier_parameter_update_strategy.get_barrier_parameter();
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         Direction& direction, HessianModel& hessian_model, RegularizationStrategy<double>& regularization_strategy,
         double trust_region_radius, WarmstartInformation& warmstart_information) {
      if (is_finite(trust_region_radius)) {
         throw std::runtime_error("The interior-point subproblem has a trust region. This is not implemented yet");
      }

      // possibly update the barrier parameter
      if (!this->first_feasibility_iteration) {
         const BarrierProblem barrier_problem(problem, this->barrier_parameter(), this->parameters);
         this->update_barrier_parameter(barrier_problem, current_iterate, current_iterate.residuals);
      }
      else {
         this->first_feasibility_iteration = false;
      }
      statistics.set("barrier", this->barrier_parameter());

      // crate the subproblem
      const BarrierProblem barrier_problem(problem, this->barrier_parameter(), this->parameters);
      const Subproblem subproblem{barrier_problem, current_iterate, hessian_model, regularization_strategy,
         trust_region_radius};

      // compute the primal-dual solution
      this->linear_solver->solve_indefinite_system(statistics, subproblem, direction, warmstart_information);
      ++this->number_subproblems_solved;

      // check whether the augmented matrix was singular, in which case the subproblem is infeasible
      if (this->linear_solver->matrix_is_singular()) {
         direction.status = SubproblemStatus::INFEASIBLE;
         return;
      }
      direction.subproblem_objective = this->evaluate_subproblem_objective(direction);

      // determine if the direction is a "small direction" (Section 3.9 of the Ipopt paper) TODO
      const bool is_small_step = InteriorPointMethod::is_small_step(problem, current_iterate.primals,
         direction.primals);
      if (is_small_step) {
         DEBUG << "This is a small step\n";
      }
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::initialize_feasibility_problem(const l1RelaxedProblem& /*problem*/, Iterate& current_iterate) {
      this->solving_feasibility_problem = true;
      this->first_feasibility_iteration = true;
      this->subproblem_definition_changed = true;

      // temporarily update the objective multiplier
      this->previous_barrier_parameter = this->barrier_parameter();
      const double new_barrier_parameter = std::max(this->barrier_parameter(), current_iterate.primal_feasibility);
      this->barrier_parameter_update_strategy.set_barrier_parameter(new_barrier_parameter);
      DEBUG << "Barrier parameter mu temporarily updated to " << this->barrier_parameter() << '\n';

      // set the bound multipliers
      /*
      for (const size_t variable_index: problem.get_lower_bounded_variables()) {
         current_iterate.multipliers.lower_bounds[variable_index] = std::min(this->default_multiplier, problem.constraint_violation_coefficient);
      }
      for (const size_t variable_index: problem.get_upper_bounded_variables()) {
         current_iterate.multipliers.upper_bounds[variable_index] = -this->default_multiplier;
      }
       */
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::exit_feasibility_problem(const OptimizationProblem& /*problem*/, Iterate& /*trial_iterate*/) {
      //assert(this->solving_feasibility_problem && "The barrier subproblem did not know it was solving the feasibility problem.");
      this->barrier_parameter_update_strategy.set_barrier_parameter(this->previous_barrier_parameter);
      this->solving_feasibility_problem = false;
      /*
      Preprocessing::compute_least_square_multipliers(problem.model, *this->linear_solver, trial_iterate,
         trial_iterate.multipliers.constraints, this->least_square_multiplier_max_norm);
      */
   }

   // set the elastic variables of the current iterate
   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) {
      DEBUG << "IPM: setting the elastic variables and their duals\n";

      for (size_t variable_index: Range(problem.number_variables)) {
         const double lower_bound = problem.variable_lower_bound(variable_index);
         const double upper_bound = problem.variable_upper_bound(variable_index);
         if (is_finite(lower_bound)) {
            current_iterate.multipliers.lower_bounds[variable_index] = this->default_multiplier;
         }
         if (is_finite(upper_bound)) {
            current_iterate.multipliers.upper_bounds[variable_index] = -this->default_multiplier;
         }
      }

      // c(x) - p + n = 0
      // analytical expression for p and n:
      // (mu_over_rho - jacobian_coefficient*this->barrier_constraints[j] + std::sqrt(radical))/2.
      // where jacobian_coefficient = -1 for p, +1 for n
      // Note: IPOPT uses a '+' sign because they define the Lagrangian as f(x) + \lambda^T c(x)
      const double mu = this->barrier_parameter();
      const auto elastic_setting_function = [&](Iterate& iterate, size_t /*constraint_index*/, size_t elastic_index, double jacobian_coefficient) {
         // precomputations
         const double constraint_j = 0.; // TODO the constraints are not stored here any more... this->constraints[constraint_index];
         const double rho = this->l1_constraint_violation_coefficient;
         const double mu_over_rho = mu / rho;
         const double radical = std::pow(constraint_j, 2) + std::pow(mu_over_rho, 2);
         const double sqrt_radical = std::sqrt(radical);

         iterate.primals[elastic_index] = (mu_over_rho - jacobian_coefficient * constraint_j + sqrt_radical) / 2.;
         iterate.multipliers.lower_bounds[elastic_index] = mu / iterate.primals[elastic_index];
         iterate.multipliers.upper_bounds[elastic_index] = 0.;
         assert(0. < iterate.primals[elastic_index] && "The elastic variable is not strictly positive.");
         assert(0. < iterate.multipliers.lower_bounds[elastic_index] && "The elastic dual is not strictly positive.");
      };
      problem.set_elastic_variable_values(current_iterate, elastic_setting_function);
   }

   template <typename BarrierProblem>
   double InteriorPointMethod<BarrierProblem>::proximal_coefficient() const {
      return std::sqrt(this->barrier_parameter());
   }

   template <typename BarrierProblem>
   EvaluationSpace& InteriorPointMethod<BarrierProblem>::get_evaluation_space() const {
      return this->linear_solver->get_evaluation_space();
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::evaluate_constraint_jacobian(const OptimizationProblem& problem, Iterate& iterate) {
      // create the subproblem
      const BarrierProblem barrier_problem(problem, this->barrier_parameter(), this->parameters);
      auto& evaluation_space = this->linear_solver->get_evaluation_space();
      evaluation_space.evaluate_constraint_jacobian(barrier_problem, iterate);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::compute_constraint_jacobian_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      const auto& evaluation_space = this->linear_solver->get_evaluation_space();
      evaluation_space.compute_constraint_jacobian_vector_product(vector, result);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::compute_constraint_jacobian_transposed_vector_product(const Vector<double>& vector, Vector<double>& result) const {
      const auto& evaluation_space = this->linear_solver->get_evaluation_space();
      evaluation_space.compute_constraint_jacobian_transposed_vector_product(vector, result);
   }

   template <typename BarrierProblem>
   double InteriorPointMethod<BarrierProblem>::compute_hessian_quadratic_product(const Vector<double>& vector) const {
      const auto& evaluation_space = this->linear_solver->get_evaluation_space();
      return evaluation_space.compute_hessian_quadratic_product(vector);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::set_auxiliary_measure(const OptimizationProblem& problem, Iterate& iterate) {
      // auxiliary measure: barrier terms
      const BarrierProblem barrier_problem(problem, this->barrier_parameter(), this->parameters);
      barrier_problem.set_auxiliary_measure(iterate);
   }

   template <typename BarrierProblem>
   double InteriorPointMethod<BarrierProblem>::compute_predicted_auxiliary_reduction_model(const OptimizationProblem& problem,
         const Iterate& current_iterate, const Vector<double>& primal_direction, double step_length) const {
      const BarrierProblem barrier_problem(problem, this->barrier_parameter(), this->parameters);
      const double directional_derivative = barrier_problem.compute_barrier_term_directional_derivative(current_iterate, primal_direction);
      return step_length * (-directional_derivative);
      // }, "α*(μ*X^{-1} e^T d)"};
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::update_barrier_parameter(const BarrierProblem& barrier_problem,
         const Iterate& current_iterate, const DualResiduals& residuals) {
      const bool barrier_parameter_updated = this->barrier_parameter_update_strategy.update_barrier_parameter(barrier_problem,
         current_iterate, residuals);
      // the barrier parameter may have been changed earlier when entering restoration
      this->subproblem_definition_changed = this->subproblem_definition_changed || barrier_parameter_updated;
   }

   // Section 3.9 in IPOPT paper
   template <typename BarrierProblem>
   bool InteriorPointMethod<BarrierProblem>::is_small_step(const OptimizationProblem& problem, const Vector<double>& current_primals,
         const Vector<double>& direction_primals) const {
      const Range variables_range = Range(problem.number_variables);
      const VectorExpression relative_direction_size{variables_range, [&](size_t variable_index) {
         return direction_primals[variable_index] / (1 + std::abs(current_primals[variable_index]));
      }};
      static double machine_epsilon = std::numeric_limits<double>::epsilon();
      return (norm_inf(relative_direction_size) <= this->parameters.small_direction_factor * machine_epsilon);
   }

   template <typename BarrierProblem>
   double InteriorPointMethod<BarrierProblem>::evaluate_subproblem_objective(const Direction& /*direction*/) const {
      return 0.; // TODO (only used in l1Relaxation at the moment)
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::postprocess_iterate(const OptimizationProblem& problem, Iterate& iterate) {
      const BarrierProblem barrier_problem(problem, this->barrier_parameter(), this->parameters);
      barrier_problem.postprocess_iterate(iterate);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::set_initial_point(const Vector<double>& /*point*/) {
      // do nothing
   }

   template <typename BarrierProblem>
   std::string InteriorPointMethod<BarrierProblem>::get_name() const {
      return "primal-dual interior-point method";
   }
} // namespace

#endif // UNO_PRIMALDUALINTERIORPOINTMETHOD_H