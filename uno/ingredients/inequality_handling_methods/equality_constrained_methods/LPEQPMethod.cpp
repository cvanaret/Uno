// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "LPEQPMethod.hpp"
#include "FixedActiveSetProblem.hpp"
#include "optimization/Iterate.hpp"
#include "ingredients/constraint_relaxation_strategies/l1RelaxedProblem.hpp"
#include "ingredients/hessian_models/HessianModelFactory.hpp"
#include "ingredients/regularization_strategies/NoRegularization.hpp"
#include "ingredients/regularization_strategies/RegularizationStrategy.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "ingredients/subproblem_solvers/LPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/QPSolverFactory.hpp"
#include "ingredients/subproblem_solvers/SubproblemStatus.hpp"
#include "optimization/Direction.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   LPEQPMethod::LPEQPMethod(const Options& options):
         InequalityHandlingMethod(),
         enforce_linear_constraints_at_initial_iterate(options.get_bool("enforce_linear_constraints")),
         LP_solver(LPSolverFactory::create(options)),
         QP_solver(QPSolverFactory::create(options)),
         activity_tolerance(options.get_double("TR_activity_tolerance")) {
   }

   void LPEQPMethod::initialize(const OptimizationProblem& problem, const HessianModel& hessian_model,
         RegularizationStrategy<double>& regularization_strategy) {
      this->initial_point.resize(problem.number_variables);
      regularization_strategy.initialize_memory(problem, hessian_model);
      this->LP_direction = Direction(problem.number_variables, problem.number_constraints);
      this->LP_solver->initialize_memory(problem, this->LP_hessian_model, this->LP_regularization_strategy);
      this->QP_solver->initialize_memory(problem, hessian_model, regularization_strategy);
   }

   void LPEQPMethod::initialize_statistics(Statistics& /*statistics*/, const Options& /*options*/) {
      // do nothing
   }

   void LPEQPMethod::generate_initial_iterate(const OptimizationProblem& /*problem*/, Iterate& /*initial_iterate*/) {
      /*
      if (this->enforce_linear_constraints_at_initial_iterate) {
         Preprocessing::enforce_linear_constraints(problem.model, initial_iterate.primals, initial_iterate.multipliers, *this->QP_solver);
      }
      */
   }

   void LPEQPMethod::solve(Statistics& statistics, const OptimizationProblem& problem, Iterate& current_iterate,
         const Multipliers& current_multipliers, Direction& direction, HessianModel& hessian_model,
         RegularizationStrategy<double>& regularization_strategy, double trust_region_radius,
         WarmstartInformation& warmstart_information) {
      //std::cout << "Before solving LP:\n"; warmstart_information.display();
      // solve LP subproblem (within trust-region to avoid unboundedness)
      if (trust_region_radius == INF<double>) {
         trust_region_radius = 100.; // TODO option
      }

      Subproblem LP_subproblem{problem, current_iterate, current_multipliers, this->LP_hessian_model,
         this->LP_regularization_strategy, trust_region_radius};
      this->solve_LP(statistics, LP_subproblem, current_multipliers, warmstart_information);
      DEBUG << "d^*(LP) = " << this->LP_direction << '\n';

      if (this->LP_direction.status == SubproblemStatus::INFEASIBLE) {
         DEBUG << "Infeasible LP, EQP direction will not be computed.\n";
         direction = this->LP_direction;
         // reset the initial point
         this->initial_point.fill(0.);
         return;
      }
      else {
         this->initial_point = this->LP_direction.primals;
      }

      // TODO compute Cauchy point

      // reformulate the problem by setting active constraints as equations and inactive bounds as +/-INF
      const FixedActiveSetProblem fixed_active_set_problem(problem, this->LP_direction.active_set, trust_region_radius);

      // compute EQP direction
      Subproblem EQP_subproblem{fixed_active_set_problem, current_iterate, current_multipliers, hessian_model,
         regularization_strategy, trust_region_radius};
      this->solve_EQP(statistics, EQP_subproblem, current_multipliers, direction, warmstart_information);
      DEBUG << "d^*(EQP) = " << direction << '\n';
      // reset the initial point
      this->initial_point.fill(0.);

      // TODO compute convex combination of the Cauchy and EQP directions

   }

   void LPEQPMethod::initialize_feasibility_problem(const l1RelaxedProblem& /*problem*/, Iterate& /*current_iterate*/) {
      // do nothing
   }

   void LPEQPMethod::exit_feasibility_problem(const OptimizationProblem& /*problem*/, Iterate& /*trial_iterate*/) {
      // do nothing
   }

   void LPEQPMethod::set_elastic_variable_values(const l1RelaxedProblem& problem, Iterate& current_iterate) {
      problem.set_elastic_variable_values(current_iterate, [&](Iterate& iterate, size_t /*j*/, size_t elastic_index, double /*jacobian_coefficient*/) {
         iterate.primals[elastic_index] = 0.;
         iterate.feasibility_multipliers.lower_bounds[elastic_index] = 1.;
         iterate.feasibility_multipliers.upper_bounds[elastic_index] = 0.;
      });
   }

   double LPEQPMethod::proximal_coefficient() const {
      return 0.;
   }

   // progress measures
   double LPEQPMethod::hessian_quadratic_product(const Vector<double>& vector) const {
      return this->QP_solver->hessian_quadratic_product(vector);
   }

   void LPEQPMethod::set_auxiliary_measure(const OptimizationProblem& /*problem*/, Iterate& iterate) {
      iterate.progress.auxiliary = 0.;
   }

   double LPEQPMethod::compute_predicted_auxiliary_reduction_model(const OptimizationProblem& /*problem*/,
         const Iterate& /*iterate*/, const Vector<double>& /*primal_direction*/, double /*step_length*/) const {
      return 0.;
   }

   void LPEQPMethod::postprocess_iterate(const OptimizationProblem& /*problem*/, Vector<double>& /*primals*/, Multipliers& /*multipliers*/) {
      // do nothing
   }

   void LPEQPMethod::set_initial_point(const Vector<double>& point) {
      // copy the point into the member
      this->initial_point = point;
   }

   std::string LPEQPMethod::get_name() const {
      return "LP-EQP";
   }

   // protected member functions
   void LPEQPMethod::solve_LP(Statistics& statistics, Subproblem& subproblem, const Multipliers& current_multipliers,
         const WarmstartInformation& warmstart_information) {
      this->LP_solver->solve(statistics, subproblem, this->initial_point, this->LP_direction, warmstart_information);
      InequalityHandlingMethod::compute_dual_displacements(current_multipliers, this->LP_direction.multipliers);
      this->number_subproblems_solved++;

      // TODO compare radius and original bound
      // reset multipliers for bound constraints active at trust region (except if one of the original bounds is active)
      for (size_t variable_index: Range(subproblem.number_variables)) {
         if (std::abs(this->LP_direction.primals[variable_index] + subproblem.trust_region_radius) <= this->activity_tolerance) {
            this->LP_direction.multipliers.lower_bounds[variable_index] = 0.;
         }
         if (std::abs(this->LP_direction.primals[variable_index] - subproblem.trust_region_radius) <= this->activity_tolerance) {
            this->LP_direction.multipliers.upper_bounds[variable_index] = 0.;
         }
      }
      this->number_subproblems_solved++;
   }
   
   void LPEQPMethod::solve_EQP(Statistics& statistics, Subproblem& subproblem, const Multipliers& current_multipliers,
         Direction& direction, const WarmstartInformation& warmstart_information) {
      this->QP_solver->solve(statistics, subproblem, this->initial_point, direction, warmstart_information);

      // TODO compare radius and original bound
      // correct EQP multipliers (the QP solver has no knowledge of the original bounds of fixed variables)
      for (size_t variable_index: this->LP_direction.active_set.bounds.at_lower_bound) {
         direction.multipliers.lower_bounds[variable_index] = direction.multipliers.upper_bounds[variable_index];
         direction.multipliers.upper_bounds[variable_index] = 0.;
      }
      for (size_t variable_index: this->LP_direction.active_set.bounds.at_upper_bound) {
         direction.multipliers.upper_bounds[variable_index] = direction.multipliers.lower_bounds[variable_index];
         direction.multipliers.lower_bounds[variable_index] = 0.;
      }
      InequalityHandlingMethod::compute_dual_displacements(current_multipliers, direction.multipliers);
      this->number_subproblems_solved++;
   }
}