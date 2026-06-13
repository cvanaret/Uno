// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "SCPSolver.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationCache.hpp"
#include "DirectSymmetricIndefiniteLinearSolver.hpp"
#include "SymmetricIndefiniteLinearSolverFactory.hpp"
#include "LinearSystem.hpp"
#include "optimization/WarmstartInformation.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "options/Options.hpp"
#include "tools/Statistics.hpp"
#include "tools/Logger.hpp"

namespace uno {
   SCPSolver::SCPSolver(size_t number_variables, size_t number_constraints, const Options& options) :
         SubproblemSolver(),
         number_variables(number_variables),
         number_constraints(number_constraints),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"))) {
   }

   void SCPSolver::initialize_memory(const Subproblem& subproblem) {
      auto& linear_system = this->linear_solver->get_linear_system();
      linear_system.initialize_augmented_system(subproblem);
      this->linear_solver->initialize_memory();
   }

   void SCPSolver::solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& initial_point, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      if (is_finite(trust_region_radius)) {
         throw std::runtime_error("SCPSolvers do not support a trust region");
      }

      auto& linear_system = this->linear_solver->get_linear_system();
      
      direction.set_dimensions(subproblem.number_variables, subproblem.number_constraints);
      direction.status = SubproblemStatus::OPTIMAL;

      if (warmstart_information.new_iterate) {
         // Populate sparsity pattern
         subproblem.evaluate_lagrangian_hessian(statistics, linear_system.matrix_values.data());
         size_t num_hessian_nonzeros = subproblem.number_hessian_nonzeros();
         subproblem.evaluate_jacobian(linear_system.matrix_values.data() + num_hessian_nonzeros, current_evaluations);
         
         // Generate Dx from the derived class
         std::vector<double> Dx(this->number_variables, 0.0);
         this->compute_diagonal_hessian(initial_point, current_evaluations, Dx);

         // We assume subproblem uses diagonal Hessian since model=zero or model=diagonal.
         // If num_hessian_nonzeros is not 0, and corresponds to the variables, we'd inject Dx here.
         // For a strictly diagonal model, the first `number_variables` elements of matrix_values map to the diagonals!
         if (num_hessian_nonzeros >= this->number_variables) {
             for (size_t j = 0; j < this->number_variables; ++j) {
                 linear_system.matrix_values[j] = Dx[j];
             }
         }

         if (!this->analysis_performed) {
            this->linear_solver->do_symbolic_analysis();
            this->analysis_performed = true;
         }
         subproblem.regularize_augmented_matrix(statistics, linear_system.matrix_values.data(), subproblem.dual_regularization_factor(), *this->linear_solver);
         subproblem.assemble_augmented_rhs(current_evaluations, linear_system.rhs);
      }

      this->linear_solver->solve_indefinite_system(linear_system.solution.data());
      subproblem.assemble_primal_dual_direction(linear_system.solution, direction);
      this->iterations++;
   }

   SolverWorkspace& SCPSolver::get_workspace() {
      return this->linear_solver->get_linear_system();
   }
} // namespace
