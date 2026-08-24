// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <cassert>
#include "WoodburyEQPSolver.hpp"
#include "LinearSystem.hpp"
#include "SymmetricIndefiniteLinearSolverFactory.hpp"
#include "ingredients/hessian_models/quasi_newton/direct/DirectQuasiNewtonHessian.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "model/Model.hpp"
#include "optimization/Direction.hpp"
#include "optimization/Evaluations.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "options/Options.hpp"
#include "symbolic/Multiplication.hpp"
#include "symbolic/Transpose.hpp"

namespace uno {
   WoodburyEQPSolver::WoodburyEQPSolver(const DirectQuasiNewtonHessian& hessian_model, const Options& options):
         SubproblemSolver(),
         hessian_model(hessian_model),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"),
            options.get_string_optional("libhsl_path").value_or(""))) {
      assert(!this->hessian_model.has_hessian_matrix());
   }

   void WoodburyEQPSolver::initialize_memory(const Subproblem& subproblem) {
      this->direction = Direction(subproblem.number_variables, subproblem.number_constraints);
      // access the linear system of the linear solver
      auto& linear_system = this->linear_solver->get_linear_system();
      linear_system.initialize_augmented_system(subproblem);
      this->linear_solver->initialize_memory();
   }

   void WoodburyEQPSolver::generate_initial_iterate(const Subproblem& subproblem, Iterate& initial_iterate,
         Evaluations& evaluations) {
      compute_least_squares_multipliers(subproblem, initial_iterate, evaluations, 1000. /* TODO add option */);
   }

   void WoodburyEQPSolver::compute_least_squares_multipliers(const Subproblem& subproblem, Iterate& iterate,
         Evaluations& evaluations, double multipliers_threshold) {
      DEBUG << "Computing least-squares multipliers\n";

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
      view(linear_system.rhs.data(), subproblem.problem.model.number_variables) = evaluations.objective_gradient;
      for (size_t variable_index: Range(subproblem.number_variables)) {
         linear_system.rhs[variable_index] -= (iterate.multipliers.lower_bounds[variable_index] +
            iterate.multipliers.upper_bounds[variable_index]);
      }

      // solve the linear system
      DEBUG3 << "KKT matrix values: " << linear_system.matrix_values << "\n\n";
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

   const Direction& WoodburyEQPSolver::solve(Statistics& statistics, const Subproblem& subproblem, const Iterate& current_iterate,
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
         // assemble the augmented matrix with only the diagonal part of the quasi-Newton Hessian
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

      // solve the linear system with only the diagonal part and store the result in solution_diagonal_part
      Vector<double> solution_diagonal_part(subproblem.number_variables + subproblem.number_constraints);
      DEBUG3 << "KKT matrix values: " << linear_system.matrix_values << "\n\n";
      this->linear_solver->solve_indefinite_system(solution_diagonal_part.data());
      if (this->linear_solver->matrix_is_singular()) {
         this->direction.status = SubproblemStatus::INFEASIBLE;
      }
      else {
         // compute the low-rank correction
         this->compute_low_rank_correction(subproblem, solution_diagonal_part);

         // assemble the full primal-dual direction
         subproblem.assemble_primal_dual_direction(current_iterate, solution_diagonal_part, this->direction);
      }
      direction.norm = norm_inf(view(direction.primals, 0, subproblem.problem.get_number_original_variables()));
      return this->direction;
   }

   bool WoodburyEQPSolver::has_second_order_corrections() const {
      return false;
   }

   void WoodburyEQPSolver::initialize_second_order_corrections(const Subproblem& /*subproblem*/, const Iterate& /*current_iterate*/,
         const Iterate& /*trial_iterate*/, Evaluations& /*current_evaluations*/, Evaluations& /*trial_evaluations*/) {
      throw std::runtime_error("No SOC implemented in WoodburyEQPSolver");
   }

   const Direction& WoodburyEQPSolver::compute_second_order_correction(const Subproblem& /*subproblem*/,
         const Iterate& /*current_iterate*/) {
      throw std::runtime_error("No SOC implemented in WoodburyEQPSolver");
   }

   void WoodburyEQPSolver::update_second_order_corrections(const Subproblem& /*subproblem*/, const Iterate& /*trial_iterate*/,
         Evaluations& /*trial_evaluations*/) {
      throw std::runtime_error("No SOC implemented in WoodburyEQPSolver");
   }

   const SolverWorkspace& WoodburyEQPSolver::get_workspace() const {
      return this->linear_solver->get_linear_system();
   }

   // protected members

   void WoodburyEQPSolver::compute_low_rank_correction(const Subproblem& subproblem, Vector<double>& b) const {
      DEBUG << "b = " << b << '\n';
      const size_t correction_rank = this->hessian_model.get_correction_rank();
      DEBUG << "Correction rank: " << correction_rank << '\n';
      if (0 < correction_rank) {
         // compute correction_rank backsolves with the correction columns as RHS
         DenseMatrix<double> E(subproblem.number_variables + subproblem.number_constraints, correction_rank);
         DenseMatrix<double> H(subproblem.number_variables + subproblem.number_constraints, correction_rank);
         // assemble the correction columns into E (E is taller than the correction matrix; the extra rows stay 0)
         for (size_t column_index: Range(correction_rank)) {
            const auto correction_column = this->hessian_model.get_correction_column(column_index);
            for (size_t row_index: Range(subproblem.problem.model.number_variables)) {
               E.entry(row_index, column_index) = correction_column[row_index];
            }
         }
         // solve A H = E with all correction columns as right-hand sides at once (column-major blocks)
         this->linear_solver->solve_indefinite_system(E.data(), H.data(), correction_rank);
         DEBUG2 << "E = " << E;
         DEBUG2 << "H = " << H;
         // compute c = Eᵀ b
         Vector<double> c(correction_rank);
         c = transpose(E)*b; // TODO move to constructor
         DEBUG2 << "c = " << c << '\n';
         // construct T = P⁻¹ + Eᵀ H: first set P⁻¹ into T, then add Eᵀ H
         DenseMatrix<double> T(correction_rank, correction_rank);
         for (size_t column_index: Range(correction_rank)) {
            T.entry(column_index, column_index) = 1./this->hessian_model.get_correction_column_scaling(column_index);
         }
         T += transpose(E)*H;
         DEBUG2 << "T = " << T;
         // solve T d = c by computing a Bunch-Kaufman factorization of T
         Vector<double> d(correction_rank);
         const bool success = WoodburyEQPSolver::solve_dense_indefinite_system(T, c, d);
         DEBUG2 << "Bunch-Kaufman success: " << success << '\n';
         DEBUG2 << "d = " << d << '\n';
         // add the correction to b: b := b - H d
         b -= H * d;
         DEBUG2 << "b = " << b << '\n';
      }
   }

   bool WoodburyEQPSolver::solve_dense_indefinite_system(DenseMatrix<double>& T, const Vector<double>& c, Vector<double>& d) {
      auto [success, ipiv] = T.compute_bunch_kaufman_factorization();
      if (success) {
         success = solve_bunch_kaufman(T, c, d, ipiv);
      }
      return success;
   }
} // namespace