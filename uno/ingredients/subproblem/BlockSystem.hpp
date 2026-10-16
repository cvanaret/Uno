// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BLOCKSYSTEM_H
#define UNO_BLOCKSYSTEM_H

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include "NewtonSystemSymbols.hpp"

namespace uno {
   // block (Row, Column) of a block matrix
   template <size_t Row, size_t Column, typename Expression>
   class Block {
   public:
      static constexpr size_t row = Row;
      static constexpr size_t column = Column;

      explicit Block(Expression&& expression): expression(std::forward<Expression>(expression)) { }

      UNO_FORWARD_ACCESSOR(get_expression, this->expression)

   protected:
      storage_t<Expression> expression;
   };

   // free function
   template <size_t Row, size_t Column, typename Expression>
   Block<Row, Column, Expression> block(Expression&& expression) {
      return Block<Row, Column, Expression>(std::forward<Expression>(expression));
   }

   // (unsymmetric) block matrix: the blocks of unknowns have given dimensions, the absent blocks are zero
   template <typename... Blocks>
   class BlockMatrix: public SymbolicExpression {
   public:
      using value_type = double;

      BlockMatrix(std::vector<size_t> block_dimensions, Blocks&&... blocks): block_dimensions(std::move(block_dimensions)),
            blocks(std::forward<Blocks>(blocks)...) {
         this->block_offsets.resize(this->block_dimensions.size() + 1);
         this->block_offsets[0] = 0;
         for (size_t block_index: Range(this->block_dimensions.size())) {
            this->block_offsets[block_index + 1] = this->block_offsets[block_index] + this->block_dimensions[block_index];
         }
         this->for_each_block([&](const auto& block) {
            using B = std::decay_t<decltype(block)>;
            if (this->number_blocks() <= B::row || this->number_blocks() <= B::column) {
               throw std::invalid_argument("BlockMatrix: block position out of range");
            }
         });
      }

      [[nodiscard]] size_t number_blocks() const { return this->block_dimensions.size(); }
      [[nodiscard]] size_t block_dimension(size_t block_index) const { return this->block_dimensions[block_index]; }
      [[nodiscard]] size_t block_offset(size_t block_index) const { return this->block_offsets[block_index]; }
      [[nodiscard]] size_t number_rows() const { return this->block_offsets.back(); }
      [[nodiscard]] size_t number_columns() const { return this->block_offsets.back(); }

      template <typename Function>
      void for_each_block(Function&& function) const {
         std::apply([&](const auto&... block) { (function(block), ...); }, this->blocks);
      }

      // block (Row, Column), found at compile time
      template <size_t Row, size_t Column>
      [[nodiscard]] const auto& get_block() const {
         constexpr size_t index = find_block<Row, Column>();
         static_assert(index < sizeof...(Blocks), "BlockMatrix::get_block: no such block");
         return std::get<index>(this->blocks).get_expression();
      }

      // y += A x
      void add_matrix_vector_product(const NewtonSystemContext& context, View<const double> x, View<double> y) const {
         this->for_each_block([&](const auto& block) {
            using B = std::decay_t<decltype(block)>;
            uno::add_matrix_vector_product(block.get_expression(), context,
               View<const double>(x.data() + this->block_offset(B::column), this->block_dimension(B::column)),
               View<double>(y.data() + this->block_offset(B::row), this->block_dimension(B::row)));
         });
      }

   protected:
      std::vector<size_t> block_dimensions;
      std::vector<size_t> block_offsets{};
      std::tuple<storage_t<Blocks>...> blocks;

      template <size_t Row, size_t Column, size_t Index = 0>
      static constexpr size_t find_block() {
         if constexpr (Index == sizeof...(Blocks)) {
            return Index;
         }
         else {
            using B = std::decay_t<std::tuple_element_t<Index, std::tuple<Blocks...>>>;
            if constexpr (B::row == Row && B::column == Column) {
               return Index;
            }
            else {
               return find_block<Row, Column, Index + 1>();
            }
         }
      }
   };

   // free function
   template <typename... Blocks>
   BlockMatrix<Blocks...> block_matrix(std::vector<size_t> block_dimensions, Blocks&&... blocks) {
      return BlockMatrix<Blocks...>(std::move(block_dimensions), std::forward<Blocks>(blocks)...);
   }

   // block vector: one component per block of unknowns. A component is a vector symbol or a NewtonRHS (residual form)
   template <typename... Components>
   class BlockVector {
   public:
      explicit BlockVector(Components&&... components): components(std::forward<Components>(components)...) { }

      // apply function(std::integral_constant<size_t, index>, component) to each component
      template <typename Function>
      void for_each_component(Function&& function) const {
         this->for_each_component(std::forward<Function>(function), std::index_sequence_for<Components...>{});
      }

   protected:
      std::tuple<storage_t<Components>...> components;

      template <typename Function, size_t... Indices>
      void for_each_component(Function&& function, std::index_sequence<Indices...>) const {
         (function(std::integral_constant<size_t, Indices>{}, std::get<Indices>(this->components)), ...);
      }
   };

   // free function
   template <typename... Components>
   BlockVector<Components...> block_vector(Components&&... components) {
      return BlockVector<Components...>(std::forward<Components>(components)...);
   }

   // result -= K(Row, Column) v for a linear term of the right-hand side of block row Row
   template <size_t Row, typename Matrix, typename Term>
   void subtract_linear_term(const Matrix& matrix, const Term& term, const NewtonSystemContext& context, View<double> result) {
      constexpr size_t column = std::decay_t<Term>::column;
      const size_t number_columns = matrix.block_dimension(column);
      const size_t number_rows = matrix.block_dimension(Row);
      std::vector<double> values(number_columns), product(number_rows, 0.);
      for (size_t index: Range(number_columns)) {
         values[index] = term.value(index);
      }
      add_matrix_vector_product(matrix.template get_block<Row, column>(), context, View<const double>(values.data(), number_columns),
         View<double>(product.data(), number_rows));
      for (size_t index: Range(number_rows)) {
         result[index] -= product[index];
      }
   }

   // uncondensed Newton system: block matrix and right-hand side
   template <typename Matrix, typename RHS>
   class KKTSystem {
   public:
      KKTSystem(Matrix&& matrix, RHS&& rhs): matrix(std::forward<Matrix>(matrix)), rhs(std::forward<RHS>(rhs)) { }

      [[nodiscard]] const auto& get_matrix() const { return this->matrix; }
      [[nodiscard]] const auto& get_rhs() const { return this->rhs; }

      // uncondensed right-hand side: for each block row, constant - sum K(row, j) v_j
      void evaluate_rhs(const NewtonSystemContext& context, View<double> result) const {
         this->rhs.for_each_component([&](auto index, const auto& component) {
            this->template evaluate_rhs_row<decltype(index)::value>(component, context, result);
         });
      }

   protected:
      storage_t<Matrix> matrix;
      storage_t<RHS> rhs;

      // note: the block row is a template parameter (not a constexpr variable of an enclosing lambda) for portability (MSVC)
      template <size_t Row, typename Component>
      void evaluate_rhs_row(const Component& component, const NewtonSystemContext& context, View<double> result) const {
         View<double> slice(result.data() + this->matrix.block_offset(Row), this->matrix.block_dimension(Row));
         rhs_constant(component).evaluate(context, slice);
         for_each_linear_term(component, [&](const auto& term) {
            subtract_linear_term<Row>(this->matrix, term, context, slice);
         });
      }
   };

   // free function
   template <typename Matrix, typename RHS>
   KKTSystem<Matrix, RHS> kkt_system(Matrix&& matrix, RHS&& rhs) {
      return {std::forward<Matrix>(matrix), std::forward<RHS>(rhs)};
   }

   // Condensation of a block system: the unknowns of the condensed blocks are eliminated. The blocks in the rows and
   // columns of a condensed block must be ScaledSelectionMatrix symbols (one nonzero per row and column), and:
   // - the diagonal block (the pivot D) is diagonal and nonsingular;
   // - each eliminated unknown v_j is coupled with at most one kept unknown u_k through a coefficient b_j (column of v_j, row k)
   //   and c_j (row of v_j, column k): kept row k: [...] + b_j Δv_j = r_k, own row j: c_j Δu_k + d_j Δv_j = r_j;
   // - two condensed blocks are not coupled with each other.
   // Eliminating Δv_j = (r_j - c_j Δu_k) / d_j yields a condensed system with a diagonal Schur complement:
   //    (k, k) += -b_j c_j / d_j,   r_k += -b_j r_j / d_j
   // Symbolic cancellation: if the right-hand side of a condensed block e is in residual form with the linear term of its
   // own pivot, r_e = r~_e - D v_e (e.g. μe - S z), the kept rows coupled with e must declare the linear term of their
   // coupling, r_k = r~_k - B v_e (e.g. +z in -∇L). The terms -B v_e and -b_j (-d_j v_j) / d_j = +B v_e cancel exactly
   // and are never evaluated: the condensed right-hand side is r~_k - B D^{-1} r~_e and Δv_j = (r~_j - c_j Δu_k) / d_j - v_j.
   // The eliminated unknowns are ordered by condensed block, then by index in the block.
   class Condensation {
   public:
      static constexpr size_t none = static_cast<size_t>(-1);

      Condensation() = default;

      template <typename Matrix, typename RHS>
      Condensation(const KKTSystem<Matrix, RHS>& system, const std::vector<size_t>& condensed_blocks) {
         const auto& matrix = system.get_matrix();
         const size_t number_blocks = matrix.number_blocks();
         this->condensed.assign(number_blocks, false);
         for (const size_t block_index: condensed_blocks) {
            if (number_blocks <= block_index) {
               throw std::invalid_argument("Condensation: condensed block out of range");
            }
            this->condensed[block_index] = true;
         }
         // offsets of the blocks in the kept and in the eliminated unknowns
         this->kept_offsets.assign(number_blocks, none);
         this->eliminated_offsets.assign(number_blocks, none);
         for (size_t block_index: Range(number_blocks)) {
            if (this->condensed[block_index]) {
               this->eliminated_offsets[block_index] = this->number_eliminated_unknowns;
               this->number_eliminated_unknowns += matrix.block_dimension(block_index);
            }
            else {
               this->kept_offsets[block_index] = this->number_kept_unknowns;
               this->number_kept_unknowns += matrix.block_dimension(block_index);
            }
         }
         const size_t number_eliminated = this->number_eliminated_unknowns;
         this->kept_row.assign(number_eliminated, none);
         this->kept_column.assign(number_eliminated, none);
         std::vector<bool> has_pivot(number_eliminated, false);

         // linear terms declared in the right-hand side: declared_terms[row][column]
         std::vector<std::vector<bool>> declared_terms(number_blocks, std::vector<bool>(number_blocks, false));
         system.get_rhs().for_each_component([&](auto index, const auto& component) {
            const size_t row = decltype(index)::value;
            for_each_linear_term(component, [&declared_terms, row](const auto& term) {
               declared_terms[row][std::decay_t<decltype(term)>::column] = true;
            });
         });
         this->cancelled.assign(number_blocks, false);
         for (size_t block_index: Range(number_blocks)) {
            this->cancelled[block_index] = this->condensed[block_index] && declared_terms[block_index][block_index];
         }

         // structure of the couplings
         matrix.for_each_block([&](const auto& block) {
            using B = std::decay_t<decltype(block)>;
            if (!this->condensed[B::row] && !this->condensed[B::column]) {
               return;
            }
            const auto& expression = block.get_expression();
            if constexpr (!is_scaled_selection_v<decltype(expression)>) {
               throw std::invalid_argument("Condensation: the blocks coupled with a condensed block must be selection or diagonal matrices");
            }
            else {
               // symbolic cancellation: the kept rows coupled with the block must declare the linear term of their coupling
               if (!this->condensed[B::row] && this->cancelled[B::column] && 0 < expression.number_nonzeros() &&
                     !declared_terms[B::row][B::column]) {
                  throw std::invalid_argument("Condensation: a kept row coupled with a condensed block in residual form must declare "
                     "the linear term of the coupling");
               }
               for (size_t entry: Range(expression.number_nonzeros())) {
                  const size_t row = expression.row_index(entry);
                  const size_t column = expression.column_index(entry);
                  if (this->condensed[B::row] && this->condensed[B::column]) { // pivot
                     if (B::row != B::column || row != column) {
                        throw std::invalid_argument("Condensation: the pivot of a condensed block must be diagonal");
                     }
                     const size_t eliminated_index = this->eliminated_offsets[B::row] + row;
                     if (has_pivot[eliminated_index]) {
                        throw std::invalid_argument("Condensation: duplicate pivot entry");
                     }
                     has_pivot[eliminated_index] = true;
                  }
                  else if (this->condensed[B::column]) { // b_j: kept row, eliminated column
                     const size_t eliminated_index = this->eliminated_offsets[B::column] + column;
                     if (this->kept_row[eliminated_index] != none) {
                        throw std::invalid_argument("Condensation: an eliminated unknown is coupled with several kept unknowns");
                     }
                     this->kept_row[eliminated_index] = this->kept_offsets[B::row] + row;
                  }
                  else { // c_j: eliminated row, kept column
                     const size_t eliminated_index = this->eliminated_offsets[B::row] + row;
                     if (this->kept_column[eliminated_index] != none) {
                        throw std::invalid_argument("Condensation: an eliminated unknown is coupled with several kept unknowns");
                     }
                     this->kept_column[eliminated_index] = this->kept_offsets[B::column] + column;
                  }
               }
            }
         });

         // positions of the diagonal Schur complement: distinct kept unknowns coupled with eliminated unknowns
         this->position.assign(number_eliminated, none);
         std::vector<size_t> position_of_kept(this->number_kept_unknowns, none);
         for (size_t eliminated_index: Range(number_eliminated)) {
            if (!has_pivot[eliminated_index]) {
               throw std::invalid_argument("Condensation: an eliminated unknown has no pivot");
            }
            const size_t row = this->kept_row[eliminated_index], column = this->kept_column[eliminated_index];
            if (row != none && column != none && row != column) {
               throw std::invalid_argument("Condensation: the Schur complement of a condensed block is not diagonal");
            }
            const size_t kept_index = (row != none) ? row : column;
            if (kept_index != none) {
               position_of_kept[kept_index] = 0; // mark
            }
         }
         for (size_t kept_index: Range(this->number_kept_unknowns)) {
            if (position_of_kept[kept_index] != none) {
               position_of_kept[kept_index] = this->kept_indices.size();
               this->kept_indices.push_back(kept_index);
            }
         }
         this->cancelled_unknown.assign(number_eliminated, false);
         for (size_t block_index: Range(number_blocks)) {
            if (this->cancelled[block_index]) {
               for (size_t index: Range(matrix.block_dimension(block_index))) {
                  this->cancelled_unknown[this->eliminated_offsets[block_index] + index] = true;
               }
            }
         }
         for (size_t eliminated_index: Range(number_eliminated)) {
            const size_t row = this->kept_row[eliminated_index], column = this->kept_column[eliminated_index];
            const size_t kept_index = (row != none) ? row : column;
            if (kept_index != none) {
               this->position[eliminated_index] = position_of_kept[kept_index];
            }
         }
         this->kept_coefficients.assign(number_eliminated, 0.);
         this->own_coefficients.assign(number_eliminated, 0.);
         this->pivots.assign(number_eliminated, 0.);
         this->eliminated_rhs.assign(number_eliminated, 0.);
         this->current_values.assign(number_eliminated, 0.);
         this->diagonal_values.assign(this->kept_indices.size(), 0.);
         this->rhs_corrections.assign(this->kept_indices.size(), 0.);
      }

      [[nodiscard]] bool is_condensed(size_t block_index) const { return this->condensed[block_index]; }
      [[nodiscard]] size_t number_kept() const { return this->number_kept_unknowns; }
      [[nodiscard]] size_t number_eliminated() const { return this->number_eliminated_unknowns; }
      [[nodiscard]] size_t kept_offset(size_t block_index) const { return this->kept_offsets[block_index]; }
      // diagonal Schur complement: one entry per kept unknown coupled with eliminated unknowns
      [[nodiscard]] size_t number_diagonal_entries() const { return this->kept_indices.size(); }
      [[nodiscard]] size_t diagonal_index(size_t position) const { return this->kept_indices[position]; }
      [[nodiscard]] double diagonal_value(size_t position) const { return this->diagonal_values[position]; }

      // evaluate the coupling and pivot coefficients at the current iterate and the diagonal Schur complement
      template <typename Matrix, typename RHS>
      void evaluate_matrix(const KKTSystem<Matrix, RHS>& system) {
         system.get_matrix().for_each_block([&](const auto& block) {
            using B = std::decay_t<decltype(block)>;
            if (!this->condensed[B::row] && !this->condensed[B::column]) {
               return;
            }
            const auto& expression = block.get_expression();
            if constexpr (is_scaled_selection_v<decltype(expression)>) {
               for (size_t entry: Range(expression.number_nonzeros())) {
                  if (this->condensed[B::row] && this->condensed[B::column]) {
                     this->pivots[this->eliminated_offsets[B::row] + expression.row_index(entry)] = expression.value(entry);
                  }
                  else if (this->condensed[B::column]) {
                     this->kept_coefficients[this->eliminated_offsets[B::column] + expression.column_index(entry)] = expression.value(entry);
                  }
                  else {
                     this->own_coefficients[this->eliminated_offsets[B::row] + expression.row_index(entry)] = expression.value(entry);
                  }
               }
            }
         });
         std::fill(this->diagonal_values.begin(), this->diagonal_values.end(), 0.);
         for (size_t eliminated_index: Range(this->number_eliminated_unknowns)) {
            if (this->position[eliminated_index] != none) {
               this->diagonal_values[this->position[eliminated_index]] += (-this->kept_coefficients[eliminated_index] *
                  this->own_coefficients[eliminated_index]) / this->pivots[eliminated_index];
            }
         }
      }

      // right-hand side of the condensed system: kept rows + corrections -b_j r_j / d_j (with symbolic cancellations).
      // Precondition: evaluate_matrix() was called at the same iterate
      template <typename Matrix, typename RHS>
      void evaluate_rhs(const KKTSystem<Matrix, RHS>& system, const NewtonSystemContext& context, View<double> rhs) {
         const auto& matrix = system.get_matrix();
         system.get_rhs().for_each_component([&](auto index, const auto& component) {
            this->template evaluate_rhs_row<decltype(index)::value>(matrix, component, context, rhs);
         });
         std::fill(this->rhs_corrections.begin(), this->rhs_corrections.end(), 0.);
         for (size_t eliminated_index: Range(this->number_eliminated_unknowns)) {
            if (this->position[eliminated_index] != none) {
               this->rhs_corrections[this->position[eliminated_index]] += (-this->kept_coefficients[eliminated_index] *
                  this->eliminated_rhs[eliminated_index]) / this->pivots[eliminated_index];
            }
         }
         for (size_t position: Range(this->kept_indices.size())) {
            rhs[this->kept_indices[position]] += this->rhs_corrections[position];
         }
      }

      // Δv_j = (r_j - c_j Δu_k) / d_j (- v_j if cancelled symbolically). Precondition: evaluate_rhs() was called at the current iterate
      void back_substitute(View<const double> kept_solution, View<double> eliminated_solution) const {
         for (size_t eliminated_index: Range(this->number_eliminated_unknowns)) {
            const size_t column = this->kept_column[eliminated_index];
            const double coupling = (column != none) ? this->own_coefficients[eliminated_index] * kept_solution[column] : 0.;
            eliminated_solution[eliminated_index] = (this->eliminated_rhs[eliminated_index] - coupling) / this->pivots[eliminated_index];
            if (this->cancelled_unknown[eliminated_index]) {
               eliminated_solution[eliminated_index] -= this->current_values[eliminated_index];
            }
         }
      }

   protected:
      std::vector<bool> condensed{};
      std::vector<bool> cancelled{}; // condensed blocks whose right-hand side is in residual form (own pivot term)
      std::vector<size_t> kept_offsets{};
      std::vector<size_t> eliminated_offsets{};
      size_t number_kept_unknowns{0};
      size_t number_eliminated_unknowns{0};
      // for each eliminated unknown
      std::vector<size_t> kept_row{}; // kept unknown in whose row it appears (b_j)
      std::vector<size_t> kept_column{}; // kept unknown that appears in its row (c_j)
      std::vector<size_t> position{}; // position of its kept unknown in the diagonal Schur complement
      std::vector<bool> cancelled_unknown{};
      std::vector<double> kept_coefficients{}; // b_j
      std::vector<double> own_coefficients{}; // c_j
      std::vector<double> pivots{}; // d_j
      std::vector<double> eliminated_rhs{}; // r_j (r~_j if cancelled symbolically)
      std::vector<double> current_values{}; // v_j (if cancelled symbolically)
      // diagonal Schur complement
      std::vector<size_t> kept_indices{};
      std::vector<double> diagonal_values{};
      std::vector<double> rhs_corrections{};

      // right-hand side of block row Row: kept rows into rhs, condensed rows into eliminated_rhs.
      // note: the block row is a template parameter (not a constexpr variable of an enclosing lambda) for portability (MSVC)
      template <size_t Row, typename Matrix, typename Component>
      void evaluate_rhs_row(const Matrix& matrix, const Component& component, const NewtonSystemContext& context, View<double> rhs) {
         const size_t dimension = matrix.block_dimension(Row);
         if (!this->condensed[Row]) {
            View<double> slice(rhs.data() + this->kept_offsets[Row], dimension);
            rhs_constant(component).evaluate(context, slice);
            for_each_linear_term(component, [&](const auto& term) {
               if (!this->cancelled[std::decay_t<decltype(term)>::column]) { // otherwise: cancelled symbolically
                  subtract_linear_term<Row>(matrix, term, context, slice);
               }
            });
         }
         else {
            View<double> slice(this->eliminated_rhs.data() + this->eliminated_offsets[Row], dimension);
            rhs_constant(component).evaluate(context, slice);
            for_each_linear_term(component, [&](const auto& term) {
               if (std::decay_t<decltype(term)>::column == Row) { // own pivot: cancelled symbolically, but its value is needed by the back substitution
                  for (size_t local_index: Range(dimension)) {
                     this->current_values[this->eliminated_offsets[Row] + local_index] = term.value(local_index);
                  }
               }
               else {
                  subtract_linear_term<Row>(matrix, term, context, slice);
               }
            });
         }
      }
   };

   // a block of eliminated unknowns, solved alone for a given kept vector u: v = D^{-1} (r - C u) (diagonal pivot D, coupling C
   // with one nonzero per row). The right-hand side is a component vector, or a NewtonRHS whose linear term is the own pivot
   // term (r = r~ - D v): then v = D^{-1} (r~ - C u) - v
   template <typename Coupling, typename Pivot, typename RHS>
   void back_substitute_block(const Coupling& coupling, const Pivot& pivot, const RHS& rhs, View<const double> kept_solution,
         View<double> solution) {
      static_assert(is_scaled_selection_v<Coupling> && is_scaled_selection_v<Pivot> && is_component_vector_v<decltype(rhs_constant(rhs))>);
      const auto& constant = rhs_constant(rhs);
      std::vector<double> coupling_terms(constant.size(), 0.);
      for (size_t entry: Range(coupling.number_nonzeros())) {
         coupling_terms[coupling.row_index(entry)] = coupling.value(entry) * kept_solution[coupling.column_index(entry)];
      }
      for (size_t entry: Range(pivot.number_nonzeros())) {
         const size_t index = pivot.row_index(entry);
         solution[index] = (constant.entry(index) - coupling_terms[index]) / pivot.value(entry);
      }
      for_each_linear_term(rhs, [&](const auto& term) {
         for (size_t index: Range(constant.size())) {
            solution[index] -= term.value(index);
         }
      });
   }

   // diagonal Schur complement of the condensed blocks, as a symbol of the condensed matrix
   class CondensedDiagonal: public SymbolicExpression {
   public:
      using value_type = double;

      CondensedDiagonal(size_t dimension, const Condensation& condensation): dimension(dimension), condensation(condensation) { }

      [[nodiscard]] size_t number_rows() const { return this->dimension; }
      [[nodiscard]] size_t number_columns() const { return this->dimension; }
      [[nodiscard]] size_t number_nonzeros() const { return this->condensation.number_diagonal_entries(); }

      void compute_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, size_t row_offset, size_t column_offset,
            uno_int solver_indexing) const {
         for (size_t position: Range(this->number_nonzeros())) {
            const size_t index = this->condensation.diagonal_index(position);
            row_indices[position] = static_cast<uno_int>(row_offset + index) + solver_indexing;
            column_indices[position] = static_cast<uno_int>(column_offset + index) + solver_indexing;
         }
      }

      // precondition: Condensation::evaluate_matrix() was called at the current iterate
      void evaluate(const NewtonSystemContext& /*context*/, View<double> values) const {
         for (size_t position: Range(this->number_nonzeros())) {
            values[position] = this->condensation.diagonal_value(position);
         }
      }

      void add_matrix_vector_product(const NewtonSystemContext& /*context*/, View<const double> x, View<double> y) const {
         for (size_t position: Range(this->number_nonzeros())) {
            const size_t index = this->condensation.diagonal_index(position);
            y[index] += this->condensation.diagonal_value(position) * x[index];
         }
      }

   protected:
      const size_t dimension;
      const Condensation& condensation;
   };

   // condensed symmetric matrix for the linear solver. The kept blocks are the primal (0) and dual (1) unknowns:
   //    [A00 + Schur   *  ]   with A00 = curvature + primal regularization, A10 = constraint Jacobian and
   //    [A10          A11 ]   A11 = dual regularization (the block (0, 1) is the transpose of A10)
   // The expression refers to the system and the condensation, which must outlive it
   template <typename Matrix, typename RHS>
   auto condensed_matrix(const KKTSystem<Matrix, RHS>& system, const Condensation& condensation) {
      const auto& matrix = system.get_matrix();
      const auto& a00 = matrix.template get_block<0, 0>();
      return symmetric_block_matrix(a00.get_left() + CondensedDiagonal(matrix.block_dimension(0), condensation) + a00.get_right(),
         matrix.template get_block<1, 0>(), matrix.template get_block<1, 1>());
   }
} // namespace

#endif // UNO_BLOCKSYSTEM_H
