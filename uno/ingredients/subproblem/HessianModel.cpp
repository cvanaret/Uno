// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include "HessianModel.hpp"
#include "linear_algebra/SymmetricMatrixFactory.hpp"
#include "reformulation/OptimizationProblem.hpp"
#include "solvers/linear/SymmetricIndefiniteLinearSolverFactory.hpp"
#include "tools/Options.hpp"
#include "tools/Statistics.hpp"

HessianModel::HessianModel(size_t dimension, size_t maximum_number_nonzeros, const std::string& sparse_format, bool use_regularization) :
      hessian(SymmetricMatrixFactory<double>::create(sparse_format, dimension, maximum_number_nonzeros, use_regularization)) {
}

// Exact Hessian
ExactHessian::ExactHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options) :
   HessianModel(dimension, maximum_number_nonzeros, options.get_string("sparse_format"), /* use_regularization = */false) {
}

void ExactHessian::evaluate(Statistics& /*statistics*/, const OptimizationProblem& problem, const Vector<double>& primal_variables,
      const Vector<double>& constraint_multipliers) {
   // evaluate Lagrangian Hessian
   this->hessian->dimension = problem.number_variables;
   problem.evaluate_lagrangian_hessian(primal_variables, constraint_multipliers, *this->hessian);
   this->evaluation_count++;
}

// Convexified Hessian
ConvexifiedHessian::ConvexifiedHessian(size_t dimension, size_t maximum_number_nonzeros, const Options& options):
      HessianModel(dimension, maximum_number_nonzeros, options.get_string("sparse_format"), /* use_regularization = */true),
      // inertia-based convexification needs a linear solver
      linear_solver(SymmetricIndefiniteLinearSolverFactory::create(options.get_string("linear_solver"), dimension, maximum_number_nonzeros)),
      regularization_initial_value(options.get_double("regularization_initial_value")),
      regularization_increase_factor(options.get_double("regularization_increase_factor")) {
}

void ConvexifiedHessian::evaluate(Statistics& statistics, const OptimizationProblem& problem, const Vector<double>& primal_variables,
      const Vector<double>& constraint_multipliers) {
   // evaluate Lagrangian Hessian
   this->hessian->dimension = problem.number_variables;
   problem.evaluate_lagrangian_hessian(primal_variables, constraint_multipliers, *this->hessian);
   this->evaluation_count++;
   // regularize (only on the original variables) to convexify the problem
   DEBUG2 << "hessian before convexification: " << *this->hessian;
   this->regularize(statistics, *this->hessian, problem.get_number_original_variables());
}

// Nocedal and Wright, p51
void ConvexifiedHessian::regularize(Statistics& statistics, SymmetricMatrix<double>& hessian, size_t number_original_variables) {
   const double smallest_diagonal_entry = hessian.smallest_diagonal_entry(number_original_variables);
   DEBUG << "The minimal diagonal entry of the matrix is " << smallest_diagonal_entry << '\n';

   double regularization_factor = (smallest_diagonal_entry <= 0.) ? this->regularization_initial_value - smallest_diagonal_entry : 0.;
   bool good_inertia = false;
   bool symbolic_factorization_performed = false;
   while (not good_inertia) {
      DEBUG << "Testing factorization with regularization factor " << regularization_factor << '\n';
      if (0. < regularization_factor) {
         hessian.set_regularization([=](size_t variable_index) {
            return (variable_index < number_original_variables) ? regularization_factor : 0.;
         });
      }
      // perform the symbolic factorization only once
      if (not symbolic_factorization_performed) {
         this->linear_solver->do_symbolic_factorization(hessian);
         symbolic_factorization_performed = true;
      }
      this->linear_solver->do_numerical_factorization(hessian);

      if (this->linear_solver->rank() == number_original_variables && this->linear_solver->number_negative_eigenvalues() == 0) {
         good_inertia = true;
         DEBUG << "Factorization was a success\n";
      }
      else {
         DEBUG << "rank: " << this->linear_solver->rank() << ", negative eigenvalues: " << this->linear_solver->number_negative_eigenvalues() << '\n';
         regularization_factor = (regularization_factor == 0.) ? this->regularization_initial_value : this->regularization_increase_factor * regularization_factor;
         assert(is_finite(regularization_factor) && "The regularization coefficient diverged");
      }
   }
   statistics.set("regularization", regularization_factor);
}