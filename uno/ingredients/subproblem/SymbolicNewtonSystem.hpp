// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_SYMBOLICNEWTONSYSTEM_H
#define UNO_SYMBOLICNEWTONSYSTEM_H

#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <utility>
#include "BlockSystem.hpp"
#include "NewtonSystem.hpp"
#include "NewtonSystemSymbols.hpp"
#include "optimization/Iterate.hpp"

// the block-system templates are only needed where Newton systems are declared (the optimization problems): the solvers
// use the abstract interface in NewtonSystem.hpp

namespace uno {
   namespace detail {
      // hash of the primal-dual values of an iterate (debug check of the call protocol of a Newton system)
      inline size_t iterate_fingerprint(const Iterate& iterate) {
         size_t hash = 0;
         const auto combine = [&hash](const auto& vector) {
            for (size_t index: Range(vector.size())) {
               hash ^= std::hash<double>{}(vector[index]) + static_cast<size_t>(0x9e3779b9u) + (hash << 6) + (hash >> 2);
            }
         };
         combine(iterate.primals);
         combine(iterate.multipliers.constraints);
         combine(iterate.multipliers.lower_bounds);
         combine(iterate.multipliers.upper_bounds);
         return hash;
      }
   } // namespace detail

   // Newton system whose uncondensed system is the expression returned by `factory(iterate)`; the blocks in Condensed are
   // eliminated. The expression only holds references and component functions: it is rebuilt (at negligible cost)
   // whenever it is used at an iterate
   template <typename Factory, typename Condensed>
   class SymbolicNewtonSystem: public NewtonSystem {
   public:
      using System = decltype(std::declval<const Factory&>()(std::declval<const Iterate&>()));
      using Structure = CondensationStructure<typename System::matrix_type, Condensed>;
      static_assert(!Condensed::contains(0) && !Condensed::contains(1),
         "SymbolicNewtonSystem: the primal (0) and dual (1) unknowns are kept");
      static_assert(Structure::template coupled_only_with<0>(),
         "SymbolicNewtonSystem: the condensed blocks can only be coupled with the primal unknowns (0)");

      SymbolicNewtonSystem(Factory factory, Condensed condensed_blocks): factory(std::move(factory)) {
         // the structure does not depend on the iterate
         const auto system = this->factory(this->empty_iterate);
         this->condensation = Condensation(system, condensed_blocks);
         const auto& block_matrix = system.get_matrix();
         if (block_matrix.number_blocks() < 2 ||
               this->condensation.number_kept() != block_matrix.block_dimension(0) + block_matrix.block_dimension(1)) {
            throw std::invalid_argument("SymbolicNewtonSystem: the kept blocks must be the primal (0) and dual (1) unknowns");
         }
         // layout of the values (the roles of the blocks are checked at compile time by condensed_matrix)
         const auto matrix = condensed_matrix(system, this->condensation);
         this->system_dimension = number_rows(matrix);
         this->system_number_nonzeros = uno::number_nonzeros(matrix);
         this->hessian_number_nonzeros = uno::number_nonzeros(matrix.get_a11().get_left());
         this->curvature_number_nonzeros = uno::number_nonzeros(matrix.get_a11().get_left().get_left());
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
#ifndef NDEBUG
         this->evaluated_iterate = detail::iterate_fingerprint(context.iterate);
#endif
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
#ifndef NDEBUG
         if (detail::iterate_fingerprint(context.iterate) != this->evaluated_iterate) {
            throw std::logic_error("NewtonSystem::assemble_rhs: the matrix was evaluated at another iterate");
         }
#endif
         const auto system = this->factory(context.iterate);
         this->condensation.evaluate_rhs(system, context, rhs);
      }

      [[nodiscard]] size_t number_eliminated() const override { return this->condensation.number_eliminated(); }

      void back_substitute(View<const double> solution, View<double> eliminated_solution) const override {
         this->condensation.back_substitute(solution, eliminated_solution);
      }

      [[nodiscard]] size_t uncondensed_dimension() const override { return this->condensation.uncondensed_dimension(); }

      // value layout: curvature, Schur complement, primal regularization, Jacobian, dual regularization
      [[nodiscard]] std::pair<size_t, size_t> schur_complement_nonzeros() const override {
         return {this->curvature_number_nonzeros, this->hessian_number_nonzeros};
      }

      void uncondense_rhs(View<const double> condensed_rhs, View<double> uncondensed_rhs) const override {
         this->condensation.uncondense_rhs(condensed_rhs, uncondensed_rhs);
      }

      void condense_rhs(View<const double> uncondensed_rhs, View<double> condensed_rhs) const override {
         this->condensation.condense_rhs(uncondensed_rhs, condensed_rhs);
      }

      void expand_solution(View<const double> uncondensed_rhs, View<const double> condensed_solution,
            View<double> uncondensed_solution) const override {
         this->condensation.expand_solution(uncondensed_rhs, condensed_solution, uncondensed_solution);
      }

      void complete_uncondensed_product(View<const double> vector, View<double> product) const override {
         this->condensation.complete_uncondensed_product(vector, product);
      }

   protected:
      const Factory factory;
      const Iterate empty_iterate{};
      // the condensation caches the evaluated coefficients between the calls of the protocol
      mutable Condensation condensation{};
      size_t system_dimension{};
      size_t system_number_nonzeros{};
      size_t hessian_number_nonzeros{}; // curvature + Schur complement
      size_t curvature_number_nonzeros{};
      size_t primal_regularization_number_nonzeros{};
      size_t jacobian_number_nonzeros{};
      size_t dual_regularization_number_nonzeros{};
      mutable size_t evaluated_iterate{0}; // fingerprint of the iterate of the last evaluation (debug check)
   };

   // free function: the factory returns the uncondensed system at an iterate; the blocks in condensed_blocks (e.g.
   // condense<2, 3>) are eliminated
   template <typename Factory, typename Condensed = CondensedBlocks<>>
   std::unique_ptr<NewtonSystem> make_newton_system(Factory&& factory, Condensed condensed_blocks = Condensed{}) {
      return std::make_unique<SymbolicNewtonSystem<std::decay_t<Factory>, Condensed>>(std::forward<Factory>(factory), condensed_blocks);
   }
} // namespace

#endif // UNO_SYMBOLICNEWTONSYSTEM_H
