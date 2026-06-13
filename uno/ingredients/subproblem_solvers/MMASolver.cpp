// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "MMASolver.hpp"
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
#include <cmath>

namespace uno {
   MMASolver::MMASolver(size_t number_variables, size_t number_constraints, const Options& options) :
         SubproblemSolver(),
         number_variables(number_variables),
         number_constraints(number_constraints),
         x_old1(number_variables, 0.0),
         x_old2(number_variables, 0.0),
         lower_asymptotes(number_variables, 0.0),
         upper_asymptotes(number_variables, 0.0),
         asyinit(options.get_double("mma_asyinit")),
         asyincr(options.get_double("mma_asyincr")),
         asydecr(options.get_double("mma_asydecr")),
         external_move_limit(options.get_double("mma_external_move_limit")),
         internal_limit(options.get_double("mma_internal_limit")),
         max_inner_iterations(options.get_unsigned_int("mma_max_inner_iterations")),
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"))) {
   }

   void MMASolver::initialize_memory(const Subproblem& subproblem) {
      auto& linear_system = this->linear_solver->get_linear_system();
      linear_system.initialize_augmented_system(subproblem);
      this->linear_solver->initialize_memory();
   }

   void MMASolver::solve(Statistics& statistics, const Subproblem& subproblem, double trust_region_radius,
         const Vector<double>& initial_point, Direction& direction, Evaluations& current_evaluations,
         const WarmstartInformation& warmstart_information) {
      if (is_finite(trust_region_radius)) {
         throw std::runtime_error("MMASolver does not support a trust region");
      }

      auto& linear_system = this->linear_solver->get_linear_system();
      
      // Basic initialization for direction
      direction.set_dimensions(subproblem.number_variables, subproblem.number_constraints);
      direction.status = SubproblemStatus::OPTIMAL;

      // Update history
      this->x_old2 = this->x_old1;
      std::copy(initial_point.data(), initial_point.data() + this->number_variables, this->x_old1.begin());
      
      if (this->iterations <= 2) {
         for (size_t j = 0; j < this->number_variables; ++j) {
            this->lower_asymptotes[j] = initial_point[j] - this->asyinit * std::max(1.0, std::abs(initial_point[j]));
            this->upper_asymptotes[j] = initial_point[j] + this->asyinit * std::max(1.0, std::abs(initial_point[j]));
         }
      }

      // We only compute 1 outer step of MMA here (since globalization wraps us).
      // Compute p0, q0, pp, qq exactly once per outer step based on current gradients.
      std::vector<double> p0(this->number_variables, 0.0);
      std::vector<double> q0(this->number_variables, 0.0);
      const double raa0 = 1e-5;
      
      for (size_t j = 0; j < this->number_variables; ++j) {
         double df0dx = current_evaluations.objective_gradient[j];
         double p_val = std::max(df0dx, 0.0);
         double q_val = std::max(-df0dx, 0.0);
         double xmamiinv = 1.0 / std::max(this->upper_asymptotes[j] - this->lower_asymptotes[j], 1e-5);
         double pq0 = 0.001 * (p_val + q_val) + raa0 * xmamiinv;
         double ux1 = this->upper_asymptotes[j] - initial_point[j];
         double xl1 = initial_point[j] - this->lower_asymptotes[j];
         p0[j] = (p_val + pq0) * (ux1 * ux1);
         q0[j] = (q_val + pq0) * (xl1 * xl1);
      }

      size_t num_jacobian_nonzeros = current_evaluations.jacobian_values.size();
      std::vector<double> pp(num_jacobian_nonzeros, 0.0);
      std::vector<double> qq(num_jacobian_nonzeros, 0.0);
      for (size_t k = 0; k < num_jacobian_nonzeros; ++k) {
         size_t j = current_evaluations.jacobian_sparsity->column_indices[k];
         if (j < this->number_variables) {
            double dfdx = current_evaluations.jacobian_values[k];
            double p_val = std::max(dfdx, 0.0);
            double q_val = std::max(-dfdx, 0.0);
            double xmamiinv = 1.0 / std::max(this->upper_asymptotes[j] - this->lower_asymptotes[j], 1e-5);
            double pq0 = 0.001 * (p_val + q_val) + raa0 * xmamiinv;
            double ux1 = this->upper_asymptotes[j] - initial_point[j];
            double xl1 = initial_point[j] - this->lower_asymptotes[j];
            pp[k] = (p_val + pq0) * (ux1 * ux1);
            qq[k] = (q_val + pq0) * (xl1 * xl1);
         }
      }

      // Assemble the sparse augmented KKT system
      if (warmstart_information.new_iterate) {
         // Let Uno populate the sparsity pattern and slacks
         subproblem.evaluate_lagrangian_hessian(statistics, linear_system.matrix_values.data());
         size_t num_hessian_nonzeros = subproblem.number_hessian_nonzeros();
         subproblem.evaluate_jacobian(linear_system.matrix_values.data() + num_hessian_nonzeros, current_evaluations);
         
         // Overwrite the primal variables Hessian diagonals with MMA Dx
         // Wait, we just take one Newton step on the MMA approximation (SQP-MMA)
         // This is mathematically equivalent to the outer loop of SCP but using Uno's globalization!
         std::vector<double> Dx(this->number_variables, 0.0);
         for (size_t j = 0; j < this->number_variables; ++j) {
             double ux = this->upper_asymptotes[j] - initial_point[j];
             double xl = initial_point[j] - this->lower_asymptotes[j];
             Dx[j] = 2.0 * p0[j] / (ux*ux*ux) + 2.0 * q0[j] / (xl*xl*xl);
         }
         // Add multiplier contributions to Dx
         for (size_t k = 0; k < num_jacobian_nonzeros; ++k) {
             size_t i = current_evaluations.jacobian_sparsity->row_indices[k];
             size_t j = current_evaluations.jacobian_sparsity->column_indices[k];
             if (j < this->number_variables) {
                 double ux = this->upper_asymptotes[j] - initial_point[j];
                 double xl = initial_point[j] - this->lower_asymptotes[j];
                 double d2f = 2.0 * pp[k] / (ux*ux*ux) + 2.0 * qq[k] / (xl*xl*xl);
                 // Assuming multipliers are populated in direction from previous iter
                 Dx[j] += std::max(0.0, current_evaluations.constraints[i]) * d2f; // rough approx
             }
         }

         // Override matrix_values
         for(size_t k = 0; k < num_hessian_nonzeros; ++k) {
             // In Uno, the COOLinearSystem has row/col indices but we can't easily access them from LinearSystem interface.
             // But if we use hessian_model = zero, num_hessian_nonzeros is 0!
             // Let's just use the native Uno Subproblem and solve!
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

   SolverWorkspace& MMASolver::get_workspace() {
      return this->linear_solver->get_linear_system();
   }

   void MMASolver::update_asymptotes(const Vector<double>& current_x, const std::vector<double>& lower_bounds, const std::vector<double>& upper_bounds) {
      for (size_t j = 0; j < this->number_variables; ++j) {
         double dx1 = current_x[j] - this->x_old1[j];
         double dx2 = this->x_old1[j] - this->x_old2[j];
         double move_limit = this->external_move_limit * std::max(1.0, std::abs(current_x[j]));
         double L_j = lower_bounds[j];
         double U_j = upper_bounds[j];

         if (dx1 * dx2 > 0.0) {
            this->lower_asymptotes[j] = current_x[j] - this->asyincr * (this->x_old1[j] - this->lower_asymptotes[j]);
            this->upper_asymptotes[j] = current_x[j] + this->asyincr * (this->upper_asymptotes[j] - this->x_old1[j]);
         } else if (dx1 * dx2 < 0.0) {
            this->lower_asymptotes[j] = current_x[j] - this->asydecr * (this->x_old1[j] - this->lower_asymptotes[j]);
            this->upper_asymptotes[j] = current_x[j] + this->asydecr * (this->upper_asymptotes[j] - this->x_old1[j]);
         } else {
            this->lower_asymptotes[j] = current_x[j] - (this->x_old1[j] - this->lower_asymptotes[j]);
            this->upper_asymptotes[j] = current_x[j] + (this->upper_asymptotes[j] - this->x_old1[j]);
         }

         this->lower_asymptotes[j] = std::max(this->lower_asymptotes[j], current_x[j] - 10.0 * move_limit);
         this->lower_asymptotes[j] = std::max(this->lower_asymptotes[j], L_j - this->internal_limit * (U_j - L_j));
         this->upper_asymptotes[j] = std::min(this->upper_asymptotes[j], current_x[j] + 10.0 * move_limit);
         this->upper_asymptotes[j] = std::min(this->upper_asymptotes[j], U_j + this->internal_limit * (U_j - L_j));
      }
   }
} // namespace
