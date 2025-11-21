// Copyright (c) 2025 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_INTERIORPOINTMETHOD_H
#define UNO_INTERIORPOINTMETHOD_H

#include <memory>
#include "../InequalityHandlingMethod.hpp"
#include "BarrierParameterUpdateStrategy.hpp"
#include "InteriorPointParameters.hpp"
#include "ingredients/constraint_relaxation_strategies/relaxed_problems/l1RelaxedProblem.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationSpace.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   template <typename BarrierProblem>
   class InteriorPointMethod : public InequalityHandlingMethod {
   public:
      explicit InteriorPointMethod(const Options& options);

      void initialize(const OptimizationProblem& problem, Iterate& current_iterate, HessianModel& hessian_model,
         InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius) override;
      void initialize_statistics(Statistics& statistics) override;
      void generate_initial_iterate(Iterate& initial_iterate) override;
      void solve(Statistics& statistics, Iterate& current_iterate, Direction& direction, double trust_region_radius,
         WarmstartInformation& warmstart_information) override;

      void initialize_feasibility_problem(Iterate& current_iterate) override;
      void set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& iterate) override;
      [[nodiscard]] double proximal_coefficient() const override;

      // matrix computations
      [[nodiscard]] EvaluationSpace& get_evaluation_space() const override;
      void evaluate_constraint_jacobian(Iterate& iterate) override;

      // acceptance
      [[nodiscard]] bool is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         UserCallbacks& user_callbacks) override;

      void postprocess_iterate(Iterate& iterate) override;

      void set_initial_point(const Vector<double>& point) override;

      [[nodiscard]] std::string get_name() const override;

   protected:
      const OptimizationProblem* problem{};
      std::unique_ptr<BarrierProblem> barrier_problem{}; // generic barrier problem
      std::unique_ptr<Subproblem> subproblem{};
      const std::unique_ptr<DirectSymmetricIndefiniteLinearSolver<double>> linear_solver;
      BarrierParameterUpdateStrategy<BarrierProblem> barrier_parameter_update_strategy;
      double previous_barrier_parameter;
      const InteriorPointParameters parameters;
      const double least_square_multiplier_max_norm;
      const double l1_constraint_violation_coefficient; // (rho in Section 3.3.1 in IPOPT paper)

      bool first_feasibility_iteration{false};

      [[nodiscard]] double barrier_parameter() const;
      void update_barrier_parameter(const Iterate& current_iterate, const DualResiduals& residuals);
      [[nodiscard]] bool is_small_step(const Vector<double>& current_primals, const Vector<double>& direction_primals) const;
      [[nodiscard]] double evaluate_subproblem_objective(const Direction& direction) const;
   };

   // class template implementation

   template <typename BarrierProblem>
   InteriorPointMethod<BarrierProblem>::InteriorPointMethod(const Options& options):
         InequalityHandlingMethod(options),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"))),
         barrier_parameter_update_strategy(options),
         previous_barrier_parameter(options.get_double("barrier_initial_parameter")),
         parameters({
               options.get_double("barrier_tau_min"),
               options.get_double("barrier_k_sigma"),
               options.get_double("barrier_regularization_exponent"),
               options.get_double("barrier_small_direction_factor"),
               options.get_double("barrier_push_variable_to_interior_k1"),
               options.get_double("barrier_push_variable_to_interior_k2"),
               options.get_double("barrier_damping_factor"),
               options.get_double("barrier_default_multiplier")
         }),
         least_square_multiplier_max_norm(options.get_double("least_square_multiplier_max_norm")),
         l1_constraint_violation_coefficient(options.get_double("l1_constraint_violation_coefficient")) {
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::initialize(const OptimizationProblem& problem, Iterate& current_iterate,
         HessianModel& hessian_model, InertiaCorrectionStrategy<double>& inertia_correction_strategy, double trust_region_radius) {
      if (trust_region_radius < INF<double>) {
         throw std::runtime_error("A trust-region radius is not supported.");
      }
      if (!problem.get_inequality_constraints().empty()) {
         throw std::runtime_error("The problem has inequality constraints. Create an instance of HomogeneousEqualityConstrainedModel");
      }
      if (!problem.get_fixed_variables().empty()) {
         throw std::runtime_error("The problem has fixed variables. Move them to the set of general constraints.");
      }
      this->problem = &problem;
      // reformulate the problem into a barrier problem
      this->barrier_problem = std::make_unique<BarrierProblem>(problem, this->parameters);
      this->barrier_problem->set_barrier_parameter(this->barrier_parameter());
      this->subproblem = std::make_unique<Subproblem>(*this->barrier_problem, current_iterate, hessian_model,
         inertia_correction_strategy);
      this->linear_solver->initialize_augmented_system(*this->subproblem);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::initialize_statistics(Statistics& statistics) {
      statistics.add_column("Barrier", Statistics::double_width, 2, Statistics::column_order.at("Barrier"));
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::generate_initial_iterate(Iterate& initial_iterate) {
      // TODO: enforce linear constraints at initial point
      // resize the initial iterate
      initial_iterate.set_number_variables(this->barrier_problem->number_variables);
      this->barrier_problem->generate_initial_iterate(initial_iterate);
      this->evaluate_progress_measures(*this->barrier_problem, initial_iterate);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::solve(Statistics& statistics, Iterate& current_iterate, Direction& direction,
         double trust_region_radius, WarmstartInformation& warmstart_information) {
      if (is_finite(trust_region_radius)) {
         throw std::runtime_error("The interior-point subproblem has a trust region. This is not implemented yet");
      }

      // possibly update the barrier parameter
      if (!this->first_feasibility_iteration) {
         this->update_barrier_parameter(current_iterate, current_iterate.residuals);
      }
      else {
         this->first_feasibility_iteration = false;
      }
      statistics.set("Barrier", this->barrier_parameter());

      // compute the primal-dual solution
      this->linear_solver->solve_indefinite_system(statistics, *this->subproblem, direction, warmstart_information);
      ++this->number_subproblems_solved;

      // check whether the augmented matrix was singular, in which case the subproblem is infeasible
      if (this->linear_solver->matrix_is_singular()) {
         direction.status = SubproblemStatus::INFEASIBLE;
         return;
      }
      direction.subproblem_objective = this->evaluate_subproblem_objective(direction);

      // determine if the direction is a "small direction" (Section 3.9 of the Ipopt paper) TODO
      if (InteriorPointMethod<BarrierProblem>::is_small_step(current_iterate.primals, direction.primals)) {
         DEBUG << "This is a small step\n";
      }
   }

   template <typename BarrierProblem>
   double InteriorPointMethod<BarrierProblem>::barrier_parameter() const {
      return this->barrier_parameter_update_strategy.get_barrier_parameter();
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::initialize_feasibility_problem(Iterate& current_iterate) {
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
   // set the elastic variables of the current iterate
   void InteriorPointMethod<BarrierProblem>::set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& iterate) {
      DEBUG << "IPM: setting the elastic variables and their duals\n";

      for (size_t variable_index: Range(problem.number_variables)) {
         const double lower_bound = problem.variable_lower_bound(variable_index);
         const double upper_bound = problem.variable_upper_bound(variable_index);
         if (is_finite(lower_bound)) {
            iterate.multipliers.lower_bounds[variable_index] = this->parameters.default_multiplier;
         }
         if (is_finite(upper_bound)) {
            iterate.multipliers.upper_bounds[variable_index] = -this->parameters.default_multiplier;
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
      problem.set_elastic_variable_values(iterate, elastic_setting_function);
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
   void InteriorPointMethod<BarrierProblem>::evaluate_constraint_jacobian(Iterate& iterate) {
      // create the subproblem
      auto& evaluation_space = this->linear_solver->get_evaluation_space();
      evaluation_space.evaluate_constraint_jacobian(*this->barrier_problem, iterate);
   }

   template <typename BarrierProblem>
   bool InteriorPointMethod<BarrierProblem>::is_iterate_acceptable(Statistics& statistics, GlobalizationStrategy& globalization_strategy,
         Iterate& current_iterate, Iterate& trial_iterate, const Direction& direction, double step_length,
         UserCallbacks& user_callbacks) {
      return InequalityHandlingMethod::is_iterate_acceptable(statistics, globalization_strategy, *this->subproblem,
         this->get_evaluation_space(), current_iterate, trial_iterate, direction, step_length, user_callbacks);
   }

   template <typename BarrierProblem>
   void InteriorPointMethod<BarrierProblem>::update_barrier_parameter(const Iterate& current_iterate, const DualResiduals& residuals) {
      const bool barrier_parameter_updated = this->barrier_parameter_update_strategy.update_barrier_parameter(*this->barrier_problem,
         current_iterate, residuals);
      // the barrier parameter may have been changed earlier when entering restoration
      this->subproblem_definition_changed = this->subproblem_definition_changed || barrier_parameter_updated;
      this->barrier_problem->set_barrier_parameter(this->barrier_parameter());
   }

   // Section 3.9 in IPOPT paper
   template <typename BarrierProblem>
   bool InteriorPointMethod<BarrierProblem>::is_small_step(const Vector<double>& current_primals, const Vector<double>& direction_primals) const {
      const Range variables_range = Range(this->problem->number_variables);
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
   void InteriorPointMethod<BarrierProblem>::postprocess_iterate(Iterate& iterate) {
      this->barrier_problem->postprocess_iterate(iterate);
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

#endif // UNO_INTERIORPOINTMETHOD_H