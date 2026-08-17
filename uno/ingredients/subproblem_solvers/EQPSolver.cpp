// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "EQPSolver.hpp"
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "SymmetricIndefiniteLinearSolverFactory.hpp"
#include "LinearSystem.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"

namespace uno {
   EQPSolver::EQPSolver(const Options& options):
         SubproblemSolver(),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"),
            options.get_string_optional("libhsl_path").value_or(""))) {
   }

   void EQPSolver::initialize_memory(const Subproblem& subproblem) {
      if (!subproblem.has_hessian_matrix()) {
         throw std::runtime_error("The subproblem does not have an explicit Hessian matrix and cannot be solved with a direct linear solver");
      }
      this->direction = Direction(subproblem.number_variables, subproblem.number_constraints);
      // access the linear system of the linear solver
      auto& linear_system = this->linear_solver->get_linear_system();
      linear_system.initialize_augmented_system(subproblem);
      this->linear_solver->initialize_memory();
   }

   void EQPSolver::generate_initial_iterate(const Subproblem& subproblem, Iterate& initial_iterate, Evaluations& evaluations) {
      compute_least_squares_multipliers(subproblem, initial_iterate, evaluations, 1000. /* TODO add option */);
   }

   void EQPSolver::compute_least_squares_multipliers(const Subproblem& subproblem, Iterate& iterate, Evaluations& evaluations,
         double multipliers_threshold) {
      INFO << "Computing least-squares multipliers at initial point\n";

      // compute least-square multipliers
      auto& linear_system = this->linear_solver->get_linear_system();

      // set up the linear system
      const size_t number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
      const size_t number_primal_inertia_correction_nonzeros = subproblem.number_variables; // full block
      const size_t number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
      const size_t number_dual_inertia_correction_nonzeros = subproblem.number_dual_inertia_correction_nonzeros();
      View hessian(linear_system.matrix_values.data(), number_hessian_nonzeros);
      View primal_inertia_correction(hessian.end(), number_primal_inertia_correction_nonzeros);
      View jacobian(primal_inertia_correction.end(), number_jacobian_nonzeros);
      View dual_inertia_correction(jacobian.end(), number_dual_inertia_correction_nonzeros);

      hessian.fill(0.); // no Hessian contribution
      primal_inertia_correction.fill(1.); // Identity block
      subproblem.evaluate_jacobian(iterate, jacobian, evaluations);
      dual_inertia_correction.fill(0.); // no dual regularization

      // perform the symbolic analysis once and for all
      if (!this->analysis_performed) {
         DEBUG << "Performing symbolic analysis of the indefinite system\n";
         this->linear_solver->do_symbolic_analysis();
         this->analysis_performed = true;
      }

      // factorize the matrix
      this->linear_solver->do_numerical_factorization(false);

      // assemble the RHS
      linear_system.rhs.fill(0.);
      evaluations.evaluate_objective_gradient(subproblem.problem.model, iterate.primals);
      view(linear_system.rhs.data(), subproblem.number_variables) = evaluations.objective_gradient;
      for (size_t variable_index: Range(subproblem.number_variables)) {
         linear_system.rhs[variable_index] -= (iterate.multipliers.lower_bounds[variable_index] +
            iterate.multipliers.upper_bounds[variable_index]);
      }

      // solve the linear system
      DEBUG << "KKT matrix values: " << linear_system.matrix_values << '\n';
      this->linear_solver->solve_indefinite_system(linear_system.solution.data());

      // set the constraint multipliers if their norm is reasonable
      const auto least_squares_multipliers = view(linear_system.solution.data(), subproblem.number_variables,
         subproblem.number_variables + subproblem.number_constraints);
      if (norm_inf(least_squares_multipliers) <= multipliers_threshold) {
         iterate.multipliers.constraints = least_squares_multipliers;
         DEBUG << "Least-squares multipliers set to " << least_squares_multipliers << '\n';
      }
      else {
         DEBUG << "Least-squares multipliers too large\n";
      }
   }

   Direction& EQPSolver::solve(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
         double trust_region_radius, const Vector<double>& /*initial_point*/, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      if (is_finite(trust_region_radius)) {
         throw std::runtime_error("The direct linear solver does not support a trust region");
      }
      this->direction.reset();

      // access the linear system
      auto& linear_system = this->linear_solver->get_linear_system();

      // set up the linear system by evaluating the functions at the current iterate
      if (warmstart_information.new_iterate) {
         // assemble the augmented matrix
         const size_t number_hessian_nonzeros = subproblem.number_hessian_nonzeros();
         const size_t number_primal_inertia_correction_nonzeros = subproblem.number_variables; // full block
         const size_t number_jacobian_nonzeros = subproblem.number_jacobian_nonzeros();
         const size_t number_dual_inertia_correction_nonzeros = subproblem.number_dual_inertia_correction_nonzeros();
         View hessian(linear_system.matrix_values.data(), number_hessian_nonzeros);
         View primal_inertia_correction(hessian.end(), number_primal_inertia_correction_nonzeros);
         View jacobian(primal_inertia_correction.end(), number_jacobian_nonzeros);
         View dual_inertia_correction(jacobian.end(), number_dual_inertia_correction_nonzeros);

         subproblem.evaluate_lagrangian_hessian(statistics, current_iterate, hessian);
         subproblem.evaluate_jacobian(current_iterate, jacobian, current_evaluations);

         // perform the symbolic analysis once and for all
         if (!this->analysis_performed) {
            DEBUG << "Performing symbolic analysis of the indefinite system\n";
            this->linear_solver->do_symbolic_analysis();
            this->analysis_performed = true;
         }

         // regularize the augmented matrix (this calls the analysis and the factorization)
         subproblem.regularize_augmented_matrix(statistics, primal_inertia_correction, dual_inertia_correction,
            subproblem.dual_regularization_factor(), *this->linear_solver);

         // assemble the RHS
         subproblem.assemble_augmented_rhs(current_iterate, current_evaluations, linear_system.rhs);
      }

      // solve the linear system
      DEBUG << "KKT matrix values: " << linear_system.matrix_values << '\n';
      this->linear_solver->solve_indefinite_system(linear_system.solution.data());
      if (this->linear_solver->matrix_is_singular()) {
         this->direction.status = SubproblemStatus::INFEASIBLE;
      }
      else {
         // assemble the full primal-dual direction
         subproblem.assemble_primal_dual_direction(current_iterate, linear_system.solution, this->direction);
      }
      return this->direction;
   }

   bool EQPSolver::has_second_order_corrections() const {
      return true;
   }

   // precondition: the constraints have been evaluated at the trial iterate in trial_evaluations
   const Direction& EQPSolver::compute_second_order_correction(const Subproblem& subproblem, const Iterate& current_iterate,
         const Vector<double>& constraints_SOC) {
      // initialize upon the first time
      if (!this->SOC_initialized) {
         this->direction_SOC = Direction(subproblem.number_variables, subproblem.number_constraints);
         this->SOC_initialized = true;
      }

      // access the linear system
      auto& linear_system = this->linear_solver->get_linear_system();

      // copy the constraints at the end of the RHS
      auto rhs_constraints = view(linear_system.rhs, subproblem.number_variables, subproblem.number_variables +
         subproblem.number_constraints);
      rhs_constraints = constraints_SOC;
      // shift the bound (lb == ub)
      for (size_t constraint_index: Range(subproblem.number_constraints)) {
         rhs_constraints[constraint_index] -= subproblem.problem.get_constraints_lower_bounds()[constraint_index];
      }
      // flip sign
      rhs_constraints.scale(-1.);

      // solve the linear system
      this->linear_solver->solve_indefinite_system(linear_system.solution.data());
      if (this->linear_solver->matrix_is_singular()) {
         this->direction_SOC.status = SubproblemStatus::INFEASIBLE;
      }
      // assemble the full primal-dual direction
      subproblem.assemble_primal_dual_direction(current_iterate, linear_system.solution, this->direction_SOC);
      return this->direction_SOC;
   }

   const SolverWorkspace& EQPSolver::get_workspace() const {
      return this->linear_solver->get_linear_system();
   }
} // namespace