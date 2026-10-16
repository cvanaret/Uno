// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NEWTONSYSTEM_H
#define UNO_NEWTONSYSTEM_H

#include <cstddef>
#include <memory>
#include <utility>
#include "NewtonSystemSymbols.hpp"
#include "Subproblem.hpp"
#include "optimization/Iterate.hpp"

namespace uno {
   // Newton system of the equality-constrained subproblem, seen by the subproblem solvers as an abstract (condensed)
   // symmetric matrix. The concrete systems are symbolic expressions built by the optimization problems
   // (see OptimizationProblem::create_newton_system), of the form
   //    [H + δ_w I    *    ]
   //    [J          -δ_c I ]
   // where H is the (possibly condensed) Hessian block, J the constraint Jacobian and δ_w, δ_c the regularization
   // blocks set by the inertia correction strategy.
   class NewtonSystem {
   public:
      virtual ~NewtonSystem() = default;

      [[nodiscard]] virtual size_t dimension() const = 0;
      [[nodiscard]] virtual size_t number_nonzeros() const = 0;
      // sparsity pattern (lower triangle) of the symmetric matrix
      virtual void compute_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const = 0;
      // values of the matrix at the current iterate (the regularization blocks are left untouched)
      virtual void evaluate(const NewtonSystemContext& context, View<double> values) const = 0;
      // views on the blocks of the matrix values (Hessian block, primal regularization, Jacobian, dual regularization)
      [[nodiscard]] virtual BlockAugmentedMatrix block_augmented_matrix(View<double> values) const = 0;
   };

   // Newton system whose matrix is the expression returned by `factory(iterate)`. The expression only holds references
   // and component functions: it is rebuilt (at negligible cost) whenever it is evaluated at an iterate
   template <typename Factory>
   class SymbolicNewtonSystem: public NewtonSystem {
   public:
      explicit SymbolicNewtonSystem(Factory factory): factory(std::move(factory)) {
         // the structure does not depend on the iterate
         const auto expression = this->factory(this->empty_iterate);
         this->system_dimension = number_rows(expression);
         this->system_number_nonzeros = uno::number_nonzeros(expression);
         this->hessian_number_nonzeros = uno::number_nonzeros(expression.get_a11().get_left());
         this->primal_regularization_number_nonzeros = uno::number_nonzeros(expression.get_a11().get_right());
         this->jacobian_number_nonzeros = uno::number_nonzeros(expression.get_a21());
         this->dual_regularization_number_nonzeros = uno::number_nonzeros(expression.get_a22());
      }

      [[nodiscard]] size_t dimension() const override { return this->system_dimension; }
      [[nodiscard]] size_t number_nonzeros() const override { return this->system_number_nonzeros; }

      void compute_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override {
         uno::compute_sparsity(this->factory(this->empty_iterate), row_indices, column_indices, 0, 0, solver_indexing);
      }

      void evaluate(const NewtonSystemContext& context, View<double> values) const override {
         uno::evaluate(this->factory(context.iterate), context, values);
      }

      [[nodiscard]] BlockAugmentedMatrix block_augmented_matrix(View<double> values) const override {
         double* hessian = values.data();
         double* primal_regularization = hessian + this->hessian_number_nonzeros;
         double* jacobian = primal_regularization + this->primal_regularization_number_nonzeros;
         double* dual_regularization = jacobian + this->jacobian_number_nonzeros;
         return {View<double>(hessian, this->hessian_number_nonzeros),
            View<double>(primal_regularization, this->primal_regularization_number_nonzeros),
            View<double>(jacobian, this->jacobian_number_nonzeros),
            View<double>(dual_regularization, this->dual_regularization_number_nonzeros)};
      }

   protected:
      const Factory factory;
      const Iterate empty_iterate{};
      size_t system_dimension{};
      size_t system_number_nonzeros{};
      size_t hessian_number_nonzeros{};
      size_t primal_regularization_number_nonzeros{};
      size_t jacobian_number_nonzeros{};
      size_t dual_regularization_number_nonzeros{};
   };

   // free function
   template <typename Factory>
   std::unique_ptr<NewtonSystem> make_newton_system(Factory&& factory) {
      return std::make_unique<SymbolicNewtonSystem<std::decay_t<Factory>>>(std::forward<Factory>(factory));
   }
} // namespace

#endif // UNO_NEWTONSYSTEM_H
