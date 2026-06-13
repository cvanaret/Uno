// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "MMASolver.hpp"
#include "optimization/Direction.hpp"
#include "optimization/EvaluationCache.hpp"
#include "LinearSystem.hpp"
#include "ingredients/subproblem/Subproblem.hpp"
#include "options/Options.hpp"
#include <cmath>

namespace uno {
   MMASolver::MMASolver(size_t number_variables, size_t number_constraints, const Options& options) :
         SCPSolver(number_variables, number_constraints, options),
         x_old1(number_variables, 0.0),
         x_old2(number_variables, 0.0),
         lower_asymptotes(number_variables, 0.0),
         upper_asymptotes(number_variables, 0.0),
         asyinit(options.get_double("mma_asyinit")),
         asyincr(options.get_double("mma_asyincr")),
         asydecr(options.get_double("mma_asydecr")),
         external_move_limit(options.get_double("mma_external_move_limit")),
         internal_limit(options.get_double("mma_internal_limit")),
         max_inner_iterations(options.get_unsigned_int("mma_max_inner_iterations")) {
   }

   void MMASolver::compute_diagonal_hessian(const Vector<double>& initial_point, Evaluations& current_evaluations, std::vector<double>& Dx) {
      // Update history
      this->x_old2 = this->x_old1;
      std::copy(initial_point.data(), initial_point.data() + this->number_variables, this->x_old1.begin());
      
      if (this->iterations <= 2) {
         for (size_t j = 0; j < this->number_variables; ++j) {
            this->lower_asymptotes[j] = initial_point[j] - this->asyinit * std::max(1.0, std::abs(initial_point[j]));
            this->upper_asymptotes[j] = initial_point[j] + this->asyinit * std::max(1.0, std::abs(initial_point[j]));
         }
      }

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

      for (size_t j = 0; j < this->number_variables; ++j) {
          double ux = this->upper_asymptotes[j] - initial_point[j];
          double xl = initial_point[j] - this->lower_asymptotes[j];
          Dx[j] = 2.0 * p0[j] / (ux*ux*ux) + 2.0 * q0[j] / (xl*xl*xl);
      }
      for (size_t k = 0; k < num_jacobian_nonzeros; ++k) {
          size_t i = current_evaluations.jacobian_sparsity->row_indices[k];
          size_t j = current_evaluations.jacobian_sparsity->column_indices[k];
          if (j < this->number_variables) {
              double ux = this->upper_asymptotes[j] - initial_point[j];
              double xl = initial_point[j] - this->lower_asymptotes[j];
              double d2f = 2.0 * pp[k] / (ux*ux*ux) + 2.0 * qq[k] / (xl*xl*xl);
              Dx[j] += std::max(0.0, current_evaluations.constraints[i]) * d2f;
          }
      }
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
