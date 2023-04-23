// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSYSTEM_H
#define UNO_SYMMETRICINDEFINITELINEARSYSTEM_H

#include <memory>
#include "SymmetricMatrix.hpp"
#include "SymmetricMatrixFactory.hpp"
#include "RectangularMatrix.hpp"
#include "optimization/Model.hpp"
#include "solvers/linear/SymmetricIndefiniteLinearSolver.hpp"
#include "tools/Options.hpp"
#include "tools/Statistics.hpp"

struct UnstableRegularization : public std::exception {

   [[nodiscard]] const char* what() const noexcept override {
      return "The inertia correction got unstable (delta_w > threshold)";
   }
};

template <typename T>
class SymmetricIndefiniteLinearSystem {
public:
   std::unique_ptr<SymmetricMatrix<T>> matrix;
   std::vector<T> rhs{};
   std::vector<T> solution{};

   SymmetricIndefiniteLinearSystem(const std::string& sparse_format, size_t max_dimension, size_t max_number_non_zeros, bool use_regularization,
         const Options& options);
   void assemble_matrix(const SymmetricMatrix<double>& hessian, const RectangularMatrix<double>& constraint_jacobian,
         size_t number_variables, size_t number_constraints);
   void factorize_matrix(const Model& model, SymmetricIndefiniteLinearSolver<T>& linear_solver);
   void regularize_matrix(Statistics& statistics, const Model& model, SymmetricIndefiniteLinearSolver<T>& linear_solver, size_t size_primal_block,
         size_t size_dual_block, T dual_regularization_parameter);
   void solve(SymmetricIndefiniteLinearSolver<T>& linear_solver);
   // [[nodiscard]] T get_primal_regularization() const;

protected:
   size_t number_factorizations{0};
   T primal_regularization{0.};
   T dual_regularization{0.};
   T previous_primal_regularization{0.};
   const T regularization_failure_threshold;
   const T primal_regularization_initial_factor;
   const T dual_regularization_fraction;
   const T primal_regularization_lb;
   const T primal_regularization_decrease_factor;
   const T primal_regularization_fast_increase_factor;
   const T primal_regularization_slow_increase_factor;
   const size_t threshold_unsuccessful_attempts;
};

template <typename T>
SymmetricIndefiniteLinearSystem<T>::SymmetricIndefiniteLinearSystem(const std::string& sparse_format, size_t max_dimension,
      size_t max_number_non_zeros, bool use_regularization, const Options& options):
      matrix(SymmetricMatrixFactory<T>::create(sparse_format, max_dimension, max_number_non_zeros, use_regularization)),
      rhs(max_dimension),
      solution(max_dimension),
      regularization_failure_threshold(T(options.get_double("regularization_failure_threshold"))),
      primal_regularization_initial_factor(T(options.get_double("primal_regularization_initial_factor"))),
      dual_regularization_fraction(T(options.get_double("dual_regularization_fraction"))),
      primal_regularization_lb(T(options.get_double("primal_regularization_lb"))),
      primal_regularization_decrease_factor(T(options.get_double("primal_regularization_decrease_factor"))),
      primal_regularization_fast_increase_factor(T(options.get_double("primal_regularization_fast_increase_factor"))),
      primal_regularization_slow_increase_factor(T(options.get_double("primal_regularization_slow_increase_factor"))),
      threshold_unsuccessful_attempts(options.get_unsigned_int("threshold_unsuccessful_attempts")) {
}

template <typename T>
void SymmetricIndefiniteLinearSystem<T>::assemble_matrix(const SymmetricMatrix<double>& hessian, const RectangularMatrix<double>& constraint_jacobian,
      size_t number_variables, size_t number_constraints) {
   this->matrix->dimension = number_variables + number_constraints;
   this->matrix->reset();
   // copy the Lagrangian Hessian in the top left block
   size_t current_column = 0;
   hessian.for_each([&](size_t i, size_t j, double entry) {
      // finalize all empty columns
      for (size_t column: Range(current_column, j)) {
         this->matrix->finalize_column(column);
         current_column++;
      }
      this->matrix->insert(entry, i, j);
   });

   // Jacobian of general constraints
   for (size_t j: Range(number_constraints)) {
      constraint_jacobian[j].for_each([&](size_t i, double derivative) {
         this->matrix->insert(derivative, i, number_variables + j);
      });
      this->matrix->finalize_column(j);
   }
}

template <typename T>
void SymmetricIndefiniteLinearSystem<T>::factorize_matrix(const Model& model, SymmetricIndefiniteLinearSolver<T>& linear_solver) {
   // compute the symbolic factorization only when:
   // the problem has a non-constant augmented system (ie is not an LP or a QP) or it is the first factorization
   if (true || this->number_factorizations == 0 || model.problem_type == NONLINEAR || not model.fixed_hessian_sparsity) {
      linear_solver.do_symbolic_factorization(*this->matrix);
   }
   linear_solver.do_numerical_factorization(*this->matrix);
   this->number_factorizations++;
}

template <typename T>
void SymmetricIndefiniteLinearSystem<T>::regularize_matrix(Statistics& statistics, const Model& model, SymmetricIndefiniteLinearSolver<T>& linear_solver,
      size_t size_primal_block, size_t size_dual_block, T dual_regularization_parameter) {
   DEBUG << "Original matrix\n" << *this->matrix << '\n';
   this->primal_regularization = T(0.);
   this->dual_regularization = T(0.);
   size_t number_attempts = 1;
   DEBUG << "Testing factorization with regularization factors (" << this->primal_regularization << ", " << this->dual_regularization << ")\n";

   if (not linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_dual_block) {
      DEBUG << "Inertia is good\n";
      statistics.add_statistic("regularization", this->primal_regularization);
      return;
   }
   auto[number_pos_eigenvalues, number_neg_eigenvalues, number_zero_eigenvalues] = linear_solver.get_inertia();
   DEBUG << "Expected inertia (" << size_primal_block << ", " << size_dual_block << ", 0), ";
   DEBUG << "got (" << number_pos_eigenvalues << ", " << number_neg_eigenvalues << ", " << number_zero_eigenvalues << ")\n";
   DEBUG << "Number of attempts: " << number_attempts << "\n\n";

   // set the constraint regularization coefficient
   if (linear_solver.matrix_is_singular()) {
      DEBUG << "Matrix is singular\n";
      this->dual_regularization = this->dual_regularization_fraction * dual_regularization_parameter;
   }
   // set the Hessian regularization coefficient
   if (this->previous_primal_regularization == 0.) {
      this->primal_regularization = this->primal_regularization_initial_factor;
   }
   else {
      this->primal_regularization = std::max(this->primal_regularization_lb, this->previous_primal_regularization / this->primal_regularization_decrease_factor);
   }

   // regularize the augmented matrix
   this->matrix->set_regularization([=](size_t i) {
      return (i < size_primal_block) ? this->primal_regularization : -this->dual_regularization;
   });

   bool good_inertia = false;
   while (not good_inertia) {
      DEBUG << "Testing factorization with regularization factors (" << this->primal_regularization << ", " << this->dual_regularization << ")\n";
      DEBUG << *this->matrix << '\n';
      this->factorize_matrix(model, linear_solver);
      number_attempts++;

      if (not linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_dual_block) {
         good_inertia = true;
         DEBUG << "Factorization was a success\n\n";
         this->previous_primal_regularization = this->primal_regularization;
      }
      else {
         auto[number_pos_eigenvalues, number_neg_eigenvalues, number_zero_eigenvalues] = linear_solver.get_inertia();
         DEBUG << "Expected inertia (" << size_primal_block << ", " << size_dual_block << ", 0), ";
         DEBUG << "got (" << number_pos_eigenvalues << ", " << number_neg_eigenvalues << ", " << number_zero_eigenvalues << ")\n";
         DEBUG << "Number of attempts: " << number_attempts << "\n";
         if (this->previous_primal_regularization == 0. || this->threshold_unsuccessful_attempts < number_attempts) {
            this->primal_regularization *= this->primal_regularization_fast_increase_factor;
         }
         else {
            this->primal_regularization *= this->primal_regularization_slow_increase_factor;
         }

         if (this->primal_regularization <= this->regularization_failure_threshold) {
            // regularize the augmented matrix
            this->matrix->set_regularization([=](size_t i) {
               return (i < size_primal_block) ? this->primal_regularization : -this->dual_regularization;
            });
         }
         else {
            throw UnstableRegularization();
         }
      }
   }
   statistics.add_statistic("regularization", this->primal_regularization);
}

template <typename T>
void SymmetricIndefiniteLinearSystem<T>::solve(SymmetricIndefiniteLinearSolver<T>& linear_solver) {
   linear_solver.solve_indefinite_system(*this->matrix, this->rhs, this->solution);
}

/*
template <typename T>
T SymmetricIndefiniteLinearSystem<T>::get_primal_regularization() const {
   return this->primal_regularization;
}
*/

#endif // UNO_SYMMETRICINDEFINITELINEARSYSTEM_H
