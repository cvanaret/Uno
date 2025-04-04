// Copyright (c) 2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "ConvexifiedHessian.hpp"
#include "ingredients/constraint_relaxation_strategies/OptimizationProblem.hpp"
#include "ingredients/hessian_models/UnstableRegularization.hpp"
#include "ingredients/subproblem_solvers/DirectSymmetricIndefiniteLinearSolver.hpp"
#include "ingredients/subproblem_solvers/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "options/Options.hpp"
#include "tools/Logger.hpp"
#include "tools/Statistics.hpp"

namespace uno {
   ConvexifiedHessian::ConvexifiedHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options):
         HessianModel(),
         // inertia-based convexification needs a linear solver
         linear_solver(SymmetricIndefiniteLinearSolverFactory::create(dimension, maximum_number_nonzeros, options)),
         regularization_initial_value(options.get_double("regularization_initial_value")),
         regularization_increase_factor(options.get_double("regularization_increase_factor")),
         regularization_failure_threshold(options.get_double("regularization_failure_threshold")) {
   }

   void ConvexifiedHessian::initialize_statistics(Statistics& statistics, const Options& options) const {
      statistics.add_column("regulariz", Statistics::double_width - 4, options.get_int("statistics_regularization_column_order"));
   }

   void ConvexifiedHessian::evaluate_hessian(Statistics& statistics, const OptimizationProblem& problem, const Vector<double>& primal_variables,
         const Vector<double>& constraint_multipliers, SymmetricMatrix<size_t, double>& hessian) {
      // evaluate Lagrangian Hessian
      hessian.set_dimension(problem.number_variables);
      problem.evaluate_lagrangian_hessian(primal_variables, constraint_multipliers, hessian);
      this->evaluation_count++;
      // regularize (only on the original variables) to convexify the problem
      this->regularize(statistics, hessian, problem.get_number_original_variables());
   }

   // evaluate_hessian() should have been called prior to calling compute_hessian_vector_product()
   // this->regularization_factor was set
   void ConvexifiedHessian::compute_hessian_vector_product(const OptimizationProblem& problem,
         const Vector<double>& vector, const Vector<double>& constraint_multipliers, Vector<double>& result) {
      // (H + lambda I)*x = H*x + lambda*x
      problem.compute_hessian_vector_product(vector, constraint_multipliers, result);
      // add regularization contribution
      for (size_t variable_index: Range(problem.get_number_original_variables())) {
         result[variable_index] += this->regularization_factor * vector[variable_index];
      }
   }

   // Nocedal and Wright, p51
   void ConvexifiedHessian::regularize(Statistics& statistics, SymmetricMatrix<size_t, double>& hessian, size_t number_original_variables) {
      DEBUG << "Current Hessian:\n" << hessian << '\n';
      const double smallest_diagonal_entry = hessian.smallest_diagonal_entry(number_original_variables);
      DEBUG << "The minimal diagonal entry of the matrix is " << smallest_diagonal_entry << '\n';

      this->regularization_factor = (smallest_diagonal_entry > 0.) ? 0. : this->regularization_initial_value - smallest_diagonal_entry;
      bool good_inertia = false;
      bool symbolic_analysis_performed = false;
      while (!good_inertia) {
         DEBUG << "Testing factorization with regularization factor " << this->regularization_factor << '\n';
         if (0. < this->regularization_factor) {
            hessian.set_regularization([=](size_t variable_index) {
               return (variable_index < number_original_variables) ? this->regularization_factor : 0.;
            });
         }
         DEBUG << "Current Hessian:\n" << hessian << '\n';

         // perform the symbolic analysis only once
         if (!symbolic_analysis_performed) {
            this->linear_solver->do_symbolic_analysis(hessian);
            symbolic_analysis_performed = true;
         }
         this->linear_solver->do_numerical_factorization(hessian);
         if (this->linear_solver->rank() == number_original_variables && this->linear_solver->number_negative_eigenvalues() == 0) {
            good_inertia = true;
            DEBUG << "Factorization was a success\n";
         }
         else {
            DEBUG << "rank: " << this->linear_solver->rank() << ", negative eigenvalues: " << this->linear_solver->number_negative_eigenvalues() << '\n';
            this->regularization_factor = (this->regularization_factor == 0.) ? this->regularization_initial_value :
               this->regularization_increase_factor * this->regularization_factor;
            if (this->regularization_factor > this->regularization_failure_threshold) {
               throw UnstableRegularization();
            }
         }
      }
      statistics.set("regulariz", this->regularization_factor);
   }
} // namespace
