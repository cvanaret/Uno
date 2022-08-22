// Copyright (c) 2022 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMMETRICINDEFINITELINEARSYSTEM_H
#define UNO_SYMMETRICINDEFINITELINEARSYSTEM_H

#include <memory>
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/SymmetricMatrixFactory.hpp"
#include "optimization/Model.hpp"
#include "solvers/linear/SymmetricIndefiniteLinearSolver.hpp"
#include "tools/Options.hpp"

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
   void factorize_matrix(const Model& model, SymmetricIndefiniteLinearSolver<T>& linear_solver);
   void regularize_matrix(const Model& model, SymmetricIndefiniteLinearSolver<T>& linear_solver, size_t size_primal_block, size_t size_dual_block,
         T dual_regularization_parameter);
   void solve(SymmetricIndefiniteLinearSolver<T>& linear_solver);

protected:
   size_t number_factorizations{0};
   T previous_primal_regularization{0.};
   const T regularization_failure_threshold;
   const T primal_regularization_initial_factor;
   const T dual_regularization_fraction;
   const T primal_regularization_lb;
   const T primal_regularization_decrease_factor;
   const T primal_regularization_fast_increase_factor;
   const T primal_regularization_slow_increase_factor;
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
      primal_regularization_slow_increase_factor(T(options.get_double("primal_regularization_slow_increase_factor"))) {
}

template <typename T>
void SymmetricIndefiniteLinearSystem<T>::factorize_matrix(const Model& model, SymmetricIndefiniteLinearSolver<T>& linear_solver) {
   // compute the symbolic factorization only when:
   // the problem has a non-constant augmented system (ie is not an LP or a QP) or it is the first factorization
   if (true || this->number_factorizations == 0 || model.problem_type == NONLINEAR || !model.fixed_hessian_sparsity) {
      linear_solver.do_symbolic_factorization(*this->matrix);
   }
   linear_solver.do_numerical_factorization(*this->matrix);
   this->number_factorizations++;
}

template <typename T>
void SymmetricIndefiniteLinearSystem<T>::regularize_matrix(const Model& model, SymmetricIndefiniteLinearSolver<T>& linear_solver,
      size_t size_primal_block, size_t size_dual_block, T dual_regularization_parameter) {
   DEBUG << "Original matrix\n" << *this->matrix << '\n';
   T primal_regularization = 0.;
   T dual_regularization = 0.;
   DEBUG << "Testing factorization with regularization factors (" << primal_regularization << ", " << dual_regularization << ")\n";

   if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_dual_block) {
      DEBUG << "Inertia is good\n";
      return;
   }
   DEBUG << "Inertia is not good\n";
   // set the constraint regularization coefficient
   if (linear_solver.matrix_is_singular()) {
      DEBUG << "Matrix is singular\n";
      dual_regularization = this->dual_regularization_fraction * dual_regularization_parameter;
   }
   else {
      dual_regularization = 0.;
   }
   // set the Hessian regularization coefficient
   if (this->previous_primal_regularization == 0.) {
      primal_regularization = this->primal_regularization_initial_factor;
   }
   else {
      primal_regularization = std::max(this->primal_regularization_lb, this->previous_primal_regularization / this->primal_regularization_decrease_factor);
   }

   // regularize the augmented matrix
   this->matrix->set_regularization([=](size_t i) {
      return (i < size_primal_block) ? primal_regularization : -dual_regularization;
   });

   bool good_inertia = false;
   while (!good_inertia) {
      DEBUG << "Testing factorization with regularization factors (" << primal_regularization << ", " << dual_regularization << ")\n";
      DEBUG << *this->matrix << '\n';
      this->factorize_matrix(model, linear_solver);

      if (!linear_solver.matrix_is_singular() && linear_solver.number_negative_eigenvalues() == size_dual_block) {
         good_inertia = true;
         DEBUG << "Factorization was a success\n\n";
         this->previous_primal_regularization = primal_regularization;
      }
      else {
         auto[number_pos_eigenvalues, number_neg_eigenvalues, number_zero_eigenvalues] = linear_solver.get_inertia();
         DEBUG << "Expected inertia (" << size_primal_block << ", " << size_dual_block << ", 0), ";
         DEBUG << "got (" << number_pos_eigenvalues << ", " << number_neg_eigenvalues << ", " << number_zero_eigenvalues << ")\n";
         if (this->previous_primal_regularization == 0.) {
            primal_regularization *= this->primal_regularization_fast_increase_factor;
         }
         else {
            primal_regularization *= this->primal_regularization_slow_increase_factor;
         }

         if (primal_regularization <= this->regularization_failure_threshold) {
            // regularize the augmented matrix
            this->matrix->set_regularization([=](size_t i) {
               return (i < size_primal_block) ? primal_regularization : -dual_regularization;
            });
         }
         else {
            throw UnstableRegularization();
         }
      }
   }
}

template <typename T>
void SymmetricIndefiniteLinearSystem<T>::solve(SymmetricIndefiniteLinearSolver<T>& linear_solver) {
   linear_solver.solve(*this->matrix, this->rhs, this->solution);
}

#endif // UNO_SYMMETRICINDEFINITELINEARSYSTEM_H
