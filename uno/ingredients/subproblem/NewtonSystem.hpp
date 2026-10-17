// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NEWTONSYSTEM_H
#define UNO_NEWTONSYSTEM_H

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include "BlockSystem.hpp"
#include "NewtonSystemSymbols.hpp"
#include "Subproblem.hpp"
#include "optimization/Iterate.hpp"

namespace uno {
   // Newton system of the equality-constrained subproblem. The optimization problems declare the uncondensed (possibly
   // unsymmetric) system as a symbolic block system (see OptimizationProblem::create_newton_system), with blocks of
   // unknowns (Δx, -Δy, ...) and a list of blocks to condense. The subproblem solvers see the condensed symmetric system
   //    [H + δ_w I    *    ]
   //    [J          -δ_c I ]
   // where H is the curvature block (+ the Schur complements of the condensed blocks), J the constraint Jacobian and δ_w,
   // δ_c the regularization blocks set by the inertia correction strategy
   class NewtonSystem {
   public:
      virtual ~NewtonSystem() = default;

      // condensed system
      [[nodiscard]] virtual size_t dimension() const = 0;
      [[nodiscard]] virtual size_t number_nonzeros() const = 0;
      // sparsity pattern (lower triangle) of the symmetric condensed matrix
      virtual void compute_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const = 0;
      // values of the condensed matrix at the current iterate (the regularization blocks are left untouched)
      virtual void evaluate(const NewtonSystemContext& context, View<double> values) const = 0;
      // views on the blocks of the matrix values (Hessian block, primal regularization, Jacobian, dual regularization)
      [[nodiscard]] virtual BlockAugmentedMatrix block_augmented_matrix(View<double> values) const = 0;
      // right-hand side of the condensed system at the current iterate (precondition: evaluate() was called at this iterate)
      virtual void assemble_rhs(const NewtonSystemContext& context, View<double> rhs) const = 0;
      // eliminated unknowns (number_eliminated() of them, ordered by condensed block) from the solution of the condensed
      // system (precondition: assemble_rhs() was called at the current iterate)
      [[nodiscard]] virtual size_t number_eliminated() const = 0;
      virtual void back_substitute(View<const double> solution, View<double> eliminated_solution) const = 0;
   };

   // Newton system whose uncondensed system is the expression returned by `factory(iterate)`. The expression only holds
   // references and component functions: it is rebuilt (at negligible cost) whenever it is used at an iterate
   template <typename Factory>
   class SymbolicNewtonSystem: public NewtonSystem {
   public:
      SymbolicNewtonSystem(Factory factory, const std::vector<size_t>& condensed_blocks): factory(std::move(factory)) {
         // the structure does not depend on the iterate
         const auto system = this->factory(this->empty_iterate);
         this->condensation = Condensation(system, condensed_blocks);
         const auto& block_matrix = system.get_matrix();
         if (block_matrix.number_blocks() < 2 || this->condensation.is_condensed(0) || this->condensation.is_condensed(1) ||
               this->condensation.number_kept() != block_matrix.block_dimension(0) + block_matrix.block_dimension(1)) {
            throw std::invalid_argument("SymbolicNewtonSystem: the kept blocks must be the primal (0) and dual (1) unknowns");
         }
         for (size_t position: Range(this->condensation.number_diagonal_entries())) {
            if (block_matrix.block_dimension(0) <= this->condensation.diagonal_index(position)) {
               throw std::invalid_argument("SymbolicNewtonSystem: the condensed blocks must be coupled with the primal unknowns");
            }
         }
         const auto matrix = condensed_matrix(system, this->condensation);
         this->system_dimension = number_rows(matrix);
         this->system_number_nonzeros = uno::number_nonzeros(matrix);
         this->hessian_number_nonzeros = uno::number_nonzeros(matrix.get_a11().get_left());
         this->primal_regularization_number_nonzeros = uno::number_nonzeros(matrix.get_a11().get_right());
         this->jacobian_number_nonzeros = uno::number_nonzeros(matrix.get_a21());
         this->dual_regularization_number_nonzeros = uno::number_nonzeros(matrix.get_a22());
      }

      [[nodiscard]] size_t dimension() const override { return this->system_dimension; }
      [[nodiscard]] size_t number_nonzeros() const override { return this->system_number_nonzeros; }

      void compute_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, uno_int solver_indexing) const override {
         const auto system = this->factory(this->empty_iterate);
         uno::compute_sparsity(condensed_matrix(system, this->condensation), row_indices, column_indices, 0, 0, solver_indexing);
      }

      void evaluate(const NewtonSystemContext& context, View<double> values) const override {
         const auto system = this->factory(context.iterate);
         this->condensation.evaluate_matrix(system);
         uno::evaluate(condensed_matrix(system, this->condensation), context, values);
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

      void assemble_rhs(const NewtonSystemContext& context, View<double> rhs) const override {
         const auto system = this->factory(context.iterate);
         this->condensation.evaluate_rhs(system, context, rhs);
      }

      [[nodiscard]] size_t number_eliminated() const override { return this->condensation.number_eliminated(); }

      void back_substitute(View<const double> solution, View<double> eliminated_solution) const override {
         this->condensation.back_substitute(solution, eliminated_solution);
      }

   protected:
      const Factory factory;
      const Iterate empty_iterate{};
      mutable Condensation condensation{};
      size_t system_dimension{};
      size_t system_number_nonzeros{};
      size_t hessian_number_nonzeros{};
      size_t primal_regularization_number_nonzeros{};
      size_t jacobian_number_nonzeros{};
      size_t dual_regularization_number_nonzeros{};
   };

   // free function: the factory returns the uncondensed system at an iterate; the blocks in condensed_blocks are eliminated
   template <typename Factory>
   std::unique_ptr<NewtonSystem> make_newton_system(Factory&& factory, const std::vector<size_t>& condensed_blocks = {}) {
      return std::make_unique<SymbolicNewtonSystem<std::decay_t<Factory>>>(std::forward<Factory>(factory), condensed_blocks);
   }
} // namespace

#endif // UNO_NEWTONSYSTEM_H
