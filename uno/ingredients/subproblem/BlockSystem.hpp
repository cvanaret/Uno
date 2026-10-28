// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_BLOCKSYSTEM_H
#define UNO_BLOCKSYSTEM_H

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
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
      using expression_type = Expression;

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

      // whether the block (Row, Column) is present (structurally nonzero)
      template <size_t Row, size_t Column>
      static constexpr bool has_block() { return find_block<Row, Column>() < sizeof...(Blocks); }

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

   // structurally zero block of a block matrix declared row by row
   struct Zero { };
   inline constexpr Zero zero{};

   // row of a block matrix declared row by row: one entry (a symbol or zero) per block column
   template <typename... Entries>
   class BlockRow {
   public:
      static constexpr size_t number_entries = sizeof...(Entries);
      // declared type of the I-th entry: a value type (owned) or an lvalue reference
      template <size_t I>
      using entry_type = std::tuple_element_t<I, std::tuple<Entries...>>;

      explicit BlockRow(Entries&&... entries): entries(std::forward<Entries>(entries)...) { }

      // the I-th entry, forwarded with its declared value category (lvalue references stay lvalue references)
      template <size_t I>
      decltype(auto) forward_entry() {
         return std::forward<entry_type<I>>(std::get<I>(this->entries));
      }

   protected:
      std::tuple<storage_t<Entries>...> entries;
   };

   template <typename... Entries>
   BlockRow<Entries...> row(Entries&&... entries) {
      return BlockRow<Entries...>(std::forward<Entries>(entries)...);
   }

   template <typename T>
   struct is_block_row: std::false_type { };
   template <typename... Entries>
   struct is_block_row<BlockRow<Entries...>>: std::true_type { };
   template <typename T>
   inline constexpr bool is_block_row_v = is_block_row<std::decay_t<T>>::value;

   namespace detail {
      // entry (Row, Column) of a block row -> tuple with the block, or empty tuple if the entry is zero
      template <size_t Row, size_t Column, typename BlockRowType>
      auto lift_entry(BlockRowType& block_row) {
         using Entry = typename BlockRowType::template entry_type<Column>;
         if constexpr (std::is_same_v<std::decay_t<Entry>, Zero>) {
            return std::tuple<>{};
         }
         else {
            return std::tuple<Block<Row, Column, Entry>>(Block<Row, Column, Entry>(block_row.template forward_entry<Column>()));
         }
      }

      template <size_t Row, typename BlockRowType, size_t... Columns>
      auto lift_row(BlockRowType& block_row, std::index_sequence<Columns...>) {
         return std::tuple_cat(lift_entry<Row, Columns>(block_row)...);
      }

      template <typename RowTuple, size_t... Rows>
      auto lift_rows(RowTuple& rows, std::index_sequence<Rows...>) {
         return std::tuple_cat(lift_row<Rows>(std::get<Rows>(rows),
            std::make_index_sequence<std::tuple_element_t<Rows, RowTuple>::number_entries>{})...);
      }
   } // namespace detail

   // free function. Two declarations:
   // - sparse: block_matrix(dimensions, block<0, 0>(A), block<1, 0>(B), ...) with the nonzero blocks in any order;
   // - dense, row by row: block_matrix(dimensions, row(A, B, zero), row(C, zero, D), ...) with explicit zero blocks.
   //   It is lowered at compile time to the sparse declaration (the zero blocks are dropped, the blocks are ordered row by row)
   template <typename... Arguments>
   auto block_matrix(std::vector<size_t> block_dimensions, Arguments&&... arguments) {
      if constexpr (0 < sizeof...(Arguments) && (is_block_row_v<Arguments> && ...)) {
         constexpr size_t number_rows = sizeof...(Arguments);
         static_assert(((std::decay_t<Arguments>::number_entries == number_rows) && ...),
            "block_matrix: each row must have one entry per block row (the block matrix is square)");
         if (block_dimensions.size() != number_rows) {
            throw std::invalid_argument("block_matrix: the number of block dimensions differs from the number of rows");
         }
         std::tuple<std::decay_t<Arguments>...> rows(std::forward<Arguments>(arguments)...);
         auto blocks = detail::lift_rows(rows, std::make_index_sequence<number_rows>{});
         return std::apply([&](auto&&... lifted_blocks) {
            return BlockMatrix<std::decay_t<decltype(lifted_blocks)>...>(std::move(block_dimensions),
               std::move(lifted_blocks)...);
         }, std::move(blocks));
      }
      else {
         static_assert(!(is_block_row_v<Arguments> || ...), "block_matrix: rows and sparse blocks cannot be mixed");
         return BlockMatrix<Arguments...>(std::move(block_dimensions), std::forward<Arguments>(arguments)...);
      }
   }


   // block vector: one component (a vector symbol) per block of unknowns
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

   // result -= K(Row, Block) v for the current value v of block Block, if the block (Row, Block) is present. The buffer
   // holds -v: its capacity is reused across calls (no allocation once it has reached the largest block dimension)
   template <size_t Row, typename Matrix, typename Current>
   void subtract_current_term(const Matrix& matrix, const Current& current, const NewtonSystemContext& context, View<double> result,
         std::vector<double>& buffer) {
      constexpr size_t column = std::decay_t<Current>::block;
      if constexpr (std::decay_t<Matrix>::template has_block<Row, column>()) {
         const size_t number_columns = matrix.block_dimension(column);
         buffer.resize(number_columns);
         for (size_t index: Range(number_columns)) {
            buffer[index] = -current.value(index);
         }
         add_matrix_vector_product(matrix.template get_block<Row, column>(), context, View<const double>(buffer.data(), number_columns),
            result);
      }
   }

   // largest block dimension of a block matrix
   template <typename Matrix>
   size_t largest_block_dimension(const Matrix& matrix) {
      size_t dimension = 0;
      for (size_t block_index: Range(matrix.number_blocks())) {
         dimension = std::max(dimension, matrix.block_dimension(block_index));
      }
      return dimension;
   }

   // uncondensed Newton system: block matrix K, right-hand side components r0 and (optionally) the current values v of some
   // blocks of unknowns. The right-hand side is in residual form: r = r0 - sum_j K(:, j) v_j over the blocks j with a
   // current value (e.g. r0 = μe and r = μe - S z in the complementarity rows of an interior-point method)
   template <typename Matrix, typename RHS, typename Currents = CurrentValues<>>
   class KKTSystem {
   public:
      using matrix_type = std::decay_t<Matrix>;
      using current_values_type = Currents;

      KKTSystem(Matrix&& matrix, RHS&& rhs, Currents current_values = Currents{}): matrix(std::forward<Matrix>(matrix)),
            rhs(std::forward<RHS>(rhs)), current_values(std::move(current_values)) { }

      [[nodiscard]] const auto& get_matrix() const { return this->matrix; }
      [[nodiscard]] const auto& get_rhs() const { return this->rhs; }
      [[nodiscard]] const Currents& get_current_values() const { return this->current_values; }

      // uncondensed right-hand side
      void evaluate_rhs(const NewtonSystemContext& context, View<double> result) const {
         std::vector<double> buffer;
         buffer.reserve(largest_block_dimension(this->matrix));
         this->rhs.for_each_component([&](auto index, const auto& component) {
            this->template evaluate_rhs_row<decltype(index)::value>(component, context, result, buffer);
         });
      }

   protected:
      storage_t<Matrix> matrix;
      storage_t<RHS> rhs;
      Currents current_values;

      // note: the block row is a template parameter (not a constexpr variable of an enclosing lambda) for portability (MSVC)
      template <size_t Row, typename Component>
      void evaluate_rhs_row(const Component& component, const NewtonSystemContext& context, View<double> result,
            std::vector<double>& buffer) const {
         View<double> slice(result.data() + this->matrix.block_offset(Row), this->matrix.block_dimension(Row));
         component.evaluate(context, slice);
         this->current_values.for_each([&](const auto& current) {
            subtract_current_term<Row>(this->matrix, current, context, slice, buffer);
         });
      }
   };

   // free function (without current values, the right-hand side is not in residual form)
   template <typename Matrix, typename RHS, typename Currents = CurrentValues<>>
   KKTSystem<Matrix, RHS, Currents> kkt_system(Matrix&& matrix, RHS&& rhs, Currents current_values = Currents{}) {
      return {std::forward<Matrix>(matrix), std::forward<RHS>(rhs), std::move(current_values)};
   }

   // list of the blocks of unknowns to condense, known at compile time: condense<2, 3>
   template <size_t... Indices>
   struct CondensedBlocks {
      static_assert(all_distinct({Indices...}), "condense: duplicate block");

      static constexpr bool contains(size_t block_index) { return ((block_index == Indices) || ...); }
      static std::vector<size_t> list() { return {Indices...}; }
   };

   template <size_t... Indices>
   inline constexpr CondensedBlocks<Indices...> condense{};

   namespace detail {
      // Blocks: the (decayed) Block types of the matrix
      template <typename Condensed, typename... Blocks>
      struct CondensationStructure;

      template <size_t... Indices, typename... Blocks>
      struct CondensationStructure<CondensedBlocks<Indices...>, Blocks...> {
         using Condensed = CondensedBlocks<Indices...>;
         static constexpr size_t no_block = static_cast<size_t>(-1);

         // the block lies in the row or the column of a condensed block
         template <typename B>
         static constexpr bool touches_condensed_block() { return Condensed::contains(B::row) || Condensed::contains(B::column); }

         // the kept block coupled with a condensed block by the block (no_block if the block does not couple them)
         template <typename B>
         static constexpr size_t kept_block_coupled() {
            if (Condensed::contains(B::row) == Condensed::contains(B::column)) {
               return no_block;
            }
            return Condensed::contains(B::row) ? B::column : B::row;
         }

         static constexpr bool has_pivot(size_t block_index) {
            return ((Blocks::row == block_index && Blocks::column == block_index) || ...);
         }

         // the blocks in the rows and columns of a condensed block are scaled partial permutations
         static constexpr bool couplings_are_partial_permutations = ((!touches_condensed_block<Blocks>() ||
            is_scaled_partial_permutation_v<typename Blocks::expression_type>) && ...);
         // two condensed blocks are not coupled with each other
         static constexpr bool condensed_blocks_are_uncoupled = ((!(Condensed::contains(Blocks::row) &&
            Condensed::contains(Blocks::column)) || Blocks::row == Blocks::column) && ...);
         // each condensed block has a pivot (its diagonal block)
         static constexpr bool pivots_are_present = (has_pivot(Indices) && ...);

         // the condensed blocks are only coupled with the kept block KeptBlock (where the Schur complement lands)
         template <size_t KeptBlock>
         static constexpr bool coupled_only_with() {
            return ((kept_block_coupled<Blocks>() == no_block || kept_block_coupled<Blocks>() == KeptBlock) && ...);
         }
      };
   } // namespace detail

   // structural (type-level) requirements of a condensation, checked at compile time
   template <typename Matrix, typename Condensed>
   struct CondensationStructure;

   template <typename... Blocks, typename Condensed>
   struct CondensationStructure<BlockMatrix<Blocks...>, Condensed>: detail::CondensationStructure<Condensed, std::decay_t<Blocks>...> { };

   // Condensation of a block system: the unknowns of the condensed blocks are eliminated. Structure (checked at compile
   // time, see CondensationStructure): the blocks in the rows and columns of a condensed block are scaled partial
   // permutations (at most one nonzero per row and column), each condensed block has a diagonal pivot block D, and two
   // condensed blocks are not coupled with each other. Structure checked at construction: the pivot blocks are diagonal,
   // and each eliminated unknown v_j is coupled with at most one kept unknown u_k through a coefficient b_j (column of v_j,
   // row k) and c_j (row of v_j, column k): kept row k: [...] + b_j Δv_j = r_k, own row j: c_j Δu_k + d_j Δv_j = r_j.
   // Eliminating Δv_j = (r_j - c_j Δu_k) / d_j yields a condensed system with a diagonal Schur complement:
   //    (k, k) += -b_j c_j / d_j,   r_k += -b_j r_j / d_j
   // (a missing coefficient is 0: one-sided and uncoupled eliminated unknowns are particular cases).
   // Symbolic cancellation: if a condensed block e has a current value v_e (residual form, see KKTSystem), its own row
   // reads r_e = r0_e - D v_e (e.g. μe - S z) and each kept row k coupled with e contains -B v_e (e.g. +z in -∇L). The terms
   // -B v_e and -b_j (-d_j v_j) / d_j = +B v_e cancel exactly and are never evaluated: the condensed right-hand side is
   // r0_k - B D^{-1} r0_e and Δv_j = (r0_j - c_j Δu_k) / d_j - v_j.
   // Numerical precondition: the pivots d_j are nonzero and finite. It is checked when the pivots are evaluated
   // (evaluate_matrix), unless the pivot block carries a sign precondition (EntrySign::positive/negative, e.g. the slacks
   // of an interior-point method), in which case the sign is only asserted in debug mode.
   // Call protocol (checked): evaluate_matrix() -> evaluate_rhs() -> back_substitute(), all at the same iterate. A new
   // evaluate_matrix() invalidates the right-hand side.
   // The eliminated unknowns are ordered by condensed block, then by index in the block. Uncondensed vectors (iterative
   // refinement on the uncondensed system) are laid out as [kept unknowns | eliminated unknowns].
   class Condensation {
   public:
      static constexpr size_t none = static_cast<size_t>(-1);

      Condensation() = default;

      template <typename Matrix, typename RHS, typename Currents, size_t... Indices>
      Condensation(const KKTSystem<Matrix, RHS, Currents>& system, CondensedBlocks<Indices...> condensed_blocks) {
         using Structure = CondensationStructure<std::decay_t<Matrix>, decltype(condensed_blocks)>;
         static_assert(Structure::couplings_are_partial_permutations,
            "Condensation: the blocks coupled with a condensed block must be scaled partial permutations (e.g. diagonal or selection matrices)");
         static_assert(Structure::condensed_blocks_are_uncoupled, "Condensation: two condensed blocks cannot be coupled");
         static_assert(Structure::pivots_are_present, "Condensation: each condensed block needs a pivot (its diagonal block)");

         // blocks: condensed, cancelled (condensed with a current value), offsets in the kept and eliminated unknowns
         const auto& matrix = system.get_matrix();
         const size_t number_blocks = matrix.number_blocks();
         this->condensed.assign(number_blocks, false);
         for (const size_t block_index: decltype(condensed_blocks)::list()) {
            if (number_blocks <= block_index) {
               throw std::invalid_argument("Condensation: condensed block out of range");
            }
            this->condensed[block_index] = true;
         }
         this->cancelled.assign(number_blocks, false);
         system.get_current_values().for_each([&](const auto& current) {
            const size_t block_index = std::decay_t<decltype(current)>::block;
            if (number_blocks <= block_index) {
               throw std::invalid_argument("Condensation: current value of a block out of range");
            }
            this->cancelled[block_index] = this->condensed[block_index];
         });
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

         // eliminated unknowns: pivot and (at most one) coupled kept unknown
         const size_t number_eliminated = this->number_eliminated_unknowns;
         this->kept_unknown.assign(number_eliminated, none);
         std::vector<bool> has_pivot(number_eliminated, false);
         const auto couple = [&](size_t eliminated_index, size_t kept_index) {
            if (this->kept_unknown[eliminated_index] != none && this->kept_unknown[eliminated_index] != kept_index) {
               throw std::invalid_argument("Condensation: an eliminated unknown is coupled with several kept unknowns (the Schur "
                  "complement would not be diagonal)");
            }
            this->kept_unknown[eliminated_index] = kept_index;
         };
         matrix.for_each_block([&](const auto& block) {
            using B = std::decay_t<decltype(block)>;
            if constexpr (Structure::template touches_condensed_block<B>()) {
               const auto& expression = block.get_expression();
               for (size_t entry: Range(expression.number_nonzeros())) {
                  const size_t row = expression.row_index(entry);
                  const size_t column = expression.column_index(entry);
                  if (this->condensed[B::row] && this->condensed[B::column]) { // pivot d_j
                     const size_t eliminated_index = this->eliminated_offsets[B::row] + row;
                     if (row != column || has_pivot[eliminated_index]) {
                        throw std::invalid_argument("Condensation: the pivot of a condensed block must be diagonal (one entry per unknown)");
                     }
                     has_pivot[eliminated_index] = true;
                  }
                  else if (this->condensed[B::column]) { // b_j: kept row, eliminated column
                     couple(this->eliminated_offsets[B::column] + column, this->kept_offsets[B::row] + row);
                  }
                  else { // c_j: eliminated row, kept column
                     couple(this->eliminated_offsets[B::row] + row, this->kept_offsets[B::column] + column);
                  }
               }
            }
         });
         if (std::find(has_pivot.begin(), has_pivot.end(), false) != has_pivot.end()) {
            throw std::invalid_argument("Condensation: an eliminated unknown has no pivot");
         }

         // diagonal Schur complement: one position per distinct coupled kept unknown, in increasing order
         std::vector<size_t> position_of_kept(this->number_kept_unknowns, none);
         for (const size_t kept_index: this->kept_unknown) {
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
         this->position.assign(number_eliminated, none);
         for (size_t eliminated_index: Range(number_eliminated)) {
            if (this->kept_unknown[eliminated_index] != none) {
               this->position[eliminated_index] = position_of_kept[this->kept_unknown[eliminated_index]];
            }
         }

         this->kept_coefficients.assign(number_eliminated, 0.);
         this->own_coefficients.assign(number_eliminated, 0.);
         this->pivots.assign(number_eliminated, 0.);
         this->eliminated_rhs.assign(number_eliminated, 0.);
         this->current_values.assign(number_eliminated, 0.);
         this->diagonal_values.assign(this->kept_indices.size(), 0.);
         this->rhs_corrections.assign(this->kept_indices.size(), 0.);
         this->term_buffer.reserve(largest_block_dimension(matrix));
      }

      [[nodiscard]] bool is_condensed(size_t block_index) const { return this->condensed[block_index]; }
      [[nodiscard]] size_t number_kept() const { return this->number_kept_unknowns; }
      [[nodiscard]] size_t number_eliminated() const { return this->number_eliminated_unknowns; }
      [[nodiscard]] size_t uncondensed_dimension() const { return this->number_kept_unknowns + this->number_eliminated_unknowns; }
      // diagonal Schur complement: one entry per kept unknown coupled with eliminated unknowns
      [[nodiscard]] size_t number_diagonal_entries() const { return this->kept_indices.size(); }
      [[nodiscard]] size_t diagonal_index(size_t position) const { return this->kept_indices[position]; }
      // precondition: evaluate_matrix() was called at the current iterate
      [[nodiscard]] double diagonal_value(size_t position) const { return this->diagonal_values[position]; }

      // evaluate the coupling and pivot coefficients at the current iterate and the diagonal Schur complement
      template <typename Matrix, typename RHS, typename Currents>
      void evaluate_matrix(const KKTSystem<Matrix, RHS, Currents>& system) {
         this->stage = Stage::structure; // invalid until the evaluation completes
         system.get_matrix().for_each_block([&](const auto& block) {
            using B = std::decay_t<decltype(block)>;
            if (!this->condensed[B::row] && !this->condensed[B::column]) {
               return;
            }
            const auto& expression = block.get_expression();
            if constexpr (is_scaled_partial_permutation_v<decltype(expression)>) {
               for (size_t entry: Range(expression.number_nonzeros())) {
                  if (this->condensed[B::row] && this->condensed[B::column]) {
                     // the pivots are the only divisors of the condensation (Schur complement, RHS, back substitution)
                     const double pivot = expression.value(entry);
                     check_pivot<std::decay_t<decltype(expression)>::sign>(pivot);
                     this->pivots[this->eliminated_offsets[B::row] + expression.row_index(entry)] = pivot;
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
         this->stage = Stage::matrix_evaluated;
      }

      // right-hand side of the condensed system: kept rows + corrections -b_j r_j / d_j (with symbolic cancellations).
      // Precondition: evaluate_matrix() was called at the same iterate
      template <typename Matrix, typename RHS, typename Currents>
      void evaluate_rhs(const KKTSystem<Matrix, RHS, Currents>& system, const NewtonSystemContext& context, View<double> rhs) {
         this->require(Stage::matrix_evaluated, "evaluate_rhs");
         system.get_rhs().for_each_component([&](auto index, const auto& component) {
            this->template evaluate_rhs_row<decltype(index)::value>(system.get_matrix(), component, system.get_current_values(),
               context, rhs);
         });
         this->compute_corrections(this->eliminated_rhs.data(), this->rhs_corrections);
         for (size_t position: Range(this->kept_indices.size())) {
            rhs[this->kept_indices[position]] += this->rhs_corrections[position];
         }
         this->stage = Stage::rhs_evaluated;
      }

      // Δv_j = (r_j - c_j Δu_k) / d_j - v_j (v_j = 0 unless cancelled symbolically). Precondition: evaluate_rhs() was called
      // at the current iterate (the kept solution may be that of a modified right-hand side whose condensed rows are
      // unchanged, e.g. a second-order correction)
      void back_substitute(View<const double> kept_solution, View<double> eliminated_solution) const {
         this->require(Stage::rhs_evaluated, "back_substitute");
         for (size_t eliminated_index: Range(this->number_eliminated_unknowns)) {
            eliminated_solution[eliminated_index] = this->substitute(eliminated_index, this->eliminated_rhs[eliminated_index],
               kept_solution) - this->current_values[eliminated_index];
         }
      }

      // operations on the uncondensed system (no residual form: the current-value terms are evaluated numerically)

      // uncondensed right-hand side from a condensed right-hand side whose eliminated rows are those of the last
      // evaluate_rhs() (e.g. the right-hand side of the system, or of a second-order correction):
      // kept rows: condensed - corrections - b_j v_j, eliminated rows: r0_j - d_j v_j
      void uncondense_rhs(View<const double> condensed_rhs, View<double> uncondensed_rhs) const {
         this->require(Stage::rhs_evaluated, "uncondense_rhs");
         const size_t number_kept = this->number_kept_unknowns;
         for (size_t kept_index: Range(number_kept)) {
            uncondensed_rhs[kept_index] = condensed_rhs[kept_index];
         }
         for (size_t position: Range(this->kept_indices.size())) {
            uncondensed_rhs[this->kept_indices[position]] -= this->rhs_corrections[position];
         }
         for (size_t eliminated_index: Range(this->number_eliminated_unknowns)) {
            const double current_value = this->current_values[eliminated_index];
            if (this->kept_unknown[eliminated_index] != none) {
               uncondensed_rhs[this->kept_unknown[eliminated_index]] -= this->kept_coefficients[eliminated_index] * current_value;
            }
            uncondensed_rhs[number_kept + eliminated_index] = this->eliminated_rhs[eliminated_index] -
               this->pivots[eliminated_index] * current_value;
         }
      }

      // condensed right-hand side of an uncondensed right-hand side: kept rows + corrections -b_j ρ_j / d_j
      void condense_rhs(View<const double> uncondensed_rhs, View<double> condensed_rhs) {
         this->require(Stage::matrix_evaluated, "condense_rhs");
         for (size_t kept_index: Range(this->number_kept_unknowns)) {
            condensed_rhs[kept_index] = uncondensed_rhs[kept_index];
         }
         this->compute_corrections(uncondensed_rhs.data() + this->number_kept_unknowns, this->correction_buffer);
         for (size_t position: Range(this->kept_indices.size())) {
            condensed_rhs[this->kept_indices[position]] += this->correction_buffer[position];
         }
      }

      // uncondensed solution of the system with an uncondensed right-hand side ρ, given the condensed solution:
      // kept unknowns: condensed solution, eliminated unknowns: Δv_j = (ρ_j - c_j Δu_k) / d_j
      void expand_solution(View<const double> uncondensed_rhs, View<const double> condensed_solution,
            View<double> uncondensed_solution) const {
         this->require(Stage::matrix_evaluated, "expand_solution");
         const size_t number_kept = this->number_kept_unknowns;
         for (size_t kept_index: Range(number_kept)) {
            uncondensed_solution[kept_index] = condensed_solution[kept_index];
         }
         for (size_t eliminated_index: Range(this->number_eliminated_unknowns)) {
            uncondensed_solution[number_kept + eliminated_index] = this->substitute(eliminated_index,
               uncondensed_rhs[number_kept + eliminated_index], condensed_solution);
         }
      }

      // product with the uncondensed matrix. On entry, the kept entries of product contain the product of the kept block
      // (as factorized, without the Schur complement) with the kept unknowns. On exit, product contains the uncondensed
      // product: kept rows + b_j u_j, eliminated rows c_j u_k + d_j u_j
      void complete_uncondensed_product(View<const double> vector, View<double> product) const {
         this->require(Stage::matrix_evaluated, "complete_uncondensed_product");
         const size_t number_kept = this->number_kept_unknowns;
         for (size_t eliminated_index: Range(this->number_eliminated_unknowns)) {
            const size_t kept_index = this->kept_unknown[eliminated_index];
            const double eliminated_value = vector[number_kept + eliminated_index];
            double& eliminated_row = product[number_kept + eliminated_index];
            eliminated_row = this->pivots[eliminated_index] * eliminated_value;
            if (kept_index != none) {
               product[kept_index] += this->kept_coefficients[eliminated_index] * eliminated_value;
               eliminated_row += this->own_coefficients[eliminated_index] * vector[kept_index];
            }
         }
      }

   protected:
      enum class Stage { structure, matrix_evaluated, rhs_evaluated };
      Stage stage{Stage::structure};
      // for each block
      std::vector<bool> condensed{};
      std::vector<bool> cancelled{}; // condensed blocks with a current value (residual form)
      std::vector<size_t> kept_offsets{};
      std::vector<size_t> eliminated_offsets{};
      size_t number_kept_unknowns{0};
      size_t number_eliminated_unknowns{0};
      // for each eliminated unknown j (a missing coupling coefficient is 0)
      std::vector<size_t> kept_unknown{}; // kept unknown k coupled with j (none if uncoupled)
      std::vector<size_t> position{}; // position of k in the diagonal Schur complement (none if uncoupled)
      std::vector<double> kept_coefficients{}; // b_j
      std::vector<double> own_coefficients{}; // c_j
      std::vector<double> pivots{}; // d_j
      std::vector<double> eliminated_rhs{}; // r_j (r0_j if cancelled symbolically)
      std::vector<double> current_values{}; // v_j if cancelled symbolically, 0 otherwise
      // diagonal Schur complement
      std::vector<size_t> kept_indices{};
      std::vector<double> diagonal_values{};
      std::vector<double> rhs_corrections{};
      std::vector<double> term_buffer{}; // current-value terms evaluated numerically (not cancelled)
      std::vector<double> correction_buffer{}; // corrections of condense_rhs()

      void require(Stage minimum_stage, const char* function_name) const {
         if (this->stage < minimum_stage) {
            throw std::logic_error(std::string("Condensation::") + function_name + ": called before " +
               (minimum_stage == Stage::matrix_evaluated ? "evaluate_matrix" : "evaluate_rhs") + " at the current matrix");
         }
      }

      // corrections -b_j e_j / d_j of the eliminated values e_j, accumulated per position of the Schur complement
      void compute_corrections(const double* eliminated_values, std::vector<double>& corrections) const {
         corrections.assign(this->kept_indices.size(), 0.);
         for (size_t eliminated_index: Range(this->number_eliminated_unknowns)) {
            if (this->position[eliminated_index] != none) {
               corrections[this->position[eliminated_index]] += (-this->kept_coefficients[eliminated_index] *
                  eliminated_values[eliminated_index]) / this->pivots[eliminated_index];
            }
         }
      }

      // (e_j - c_j u_k) / d_j
      [[nodiscard]] double substitute(size_t eliminated_index, double eliminated_value, View<const double> kept_solution) const {
         const size_t kept_index = this->kept_unknown[eliminated_index];
         const double coupling = (kept_index != none) ? this->own_coefficients[eliminated_index] * kept_solution[kept_index] : 0.;
         return (eliminated_value - coupling) / this->pivots[eliminated_index];
      }

      // right-hand side of block row Row (kept rows into rhs, condensed rows into eliminated_rhs): constant part, then the
      // current-value terms -K(Row, j) v_j: those of a cancelled block j are skipped (the own term of a condensed row is
      // kept for the back substitution), the others are evaluated.
      // note: the block row is a template parameter (not a constexpr variable of an enclosing lambda) for portability (MSVC)
      template <size_t Row, typename Matrix, typename Component, typename Currents>
      void evaluate_rhs_row(const Matrix& matrix, const Component& component, const Currents& currents, const NewtonSystemContext& context,
            View<double> rhs) {
         View<double> slice = this->condensed[Row] ?
            View<double>(this->eliminated_rhs.data() + this->eliminated_offsets[Row], matrix.block_dimension(Row)) :
            View<double>(rhs.data() + this->kept_offsets[Row], matrix.block_dimension(Row));
         component.evaluate(context, slice);
         currents.for_each([&](const auto& current) {
            this->template current_term<Row>(matrix, current, context, slice);
         });
      }

      template <size_t Row, typename Matrix, typename Current>
      void current_term(const Matrix& matrix, const Current& current, const NewtonSystemContext& context, View<double> slice) {
         constexpr size_t block = std::decay_t<Current>::block;
         if (!this->cancelled[block]) {
            subtract_current_term<Row>(matrix, current, context, slice, this->term_buffer);
         }
         else if constexpr (block == Row) {
            for (size_t local_index: Range(matrix.block_dimension(Row))) {
               this->current_values[this->eliminated_offsets[Row] + local_index] = current.value(local_index);
            }
         }
      }
   };

   // a block of eliminated unknowns, solved alone for a given kept vector u: v = D^{-1} (r0 - C u) (diagonal pivot D,
   // coupling C with one nonzero per row, right-hand side r0), or v = D^{-1} (r0 - C u) - v_current if the block has a
   // current value (residual form r = r0 - D v_current, see Condensation)
   template <typename Coupling, typename Pivot, typename Constant, typename Current = std::nullptr_t>
   void back_substitute_block(const Coupling& coupling, const Pivot& pivot, const Constant& constant, View<const double> kept_solution,
         View<double> solution, const Current& current = nullptr) {
      static_assert(is_scaled_partial_permutation_v<Coupling> && is_scaled_partial_permutation_v<Pivot> && is_component_vector_v<Constant>);
      std::vector<double> coupling_terms(constant.size(), 0.);
      for (size_t entry: Range(coupling.number_nonzeros())) {
         coupling_terms[coupling.row_index(entry)] = coupling.value(entry) * kept_solution[coupling.column_index(entry)];
      }
      for (size_t entry: Range(pivot.number_nonzeros())) {
         const size_t index = pivot.row_index(entry);
         const double pivot_value = pivot.value(entry);
         check_pivot<std::decay_t<Pivot>::sign>(pivot_value);
         solution[index] = (constant.entry(index) - coupling_terms[index]) / pivot_value;
      }
      if constexpr (!std::is_same_v<Current, std::nullptr_t>) {
         for (size_t index: Range(constant.size())) {
            solution[index] -= current(index);
         }
      }
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

   // roles of the blocks of a primal-dual Newton system, checked at compile time
   template <typename T>
   struct curvature_plus_primal_regularization: std::false_type { };
   template <typename Curvature, int Enable>
   struct curvature_plus_primal_regularization<Sum<Curvature, PrimalRegularization, Enable>>: std::true_type { };

   // condensed symmetric matrix for the linear solver. The kept blocks are the primal (0) and dual (1) unknowns:
   //    [A00 + Schur   *  ]   with A00 = curvature + primal regularization, A10 = constraint Jacobian and
   //    [A10          A11 ]   A11 = dual regularization (the block (0, 1) is the transpose of A10)
   // Value layout: curvature, Schur complement, primal regularization, Jacobian, dual regularization.
   // The expression refers to the system and the condensation, which must outlive it
   template <typename Matrix, typename RHS, typename Currents>
   auto condensed_matrix(const KKTSystem<Matrix, RHS, Currents>& system, const Condensation& condensation) {
      const auto& matrix = system.get_matrix();
      const auto& a00 = matrix.template get_block<0, 0>();
      static_assert(curvature_plus_primal_regularization<std::decay_t<decltype(a00)>>::value,
         "Newton system: the block (0, 0) must be declared as curvature + PrimalRegularization");
      static_assert(std::is_same_v<std::decay_t<decltype(matrix.template get_block<1, 0>())>, ConstraintJacobian>,
         "Newton system: the block (1, 0) must be the ConstraintJacobian");
      static_assert(std::is_same_v<std::decay_t<decltype(matrix.template get_block<1, 1>())>, DualRegularization>,
         "Newton system: the block (1, 1) must be the DualRegularization");
      return symmetric_block_matrix(a00.get_left() + CondensedDiagonal(matrix.block_dimension(0), condensation) + a00.get_right(),
         matrix.template get_block<1, 0>(), matrix.template get_block<1, 1>());
   }
} // namespace

#endif // UNO_BLOCKSYSTEM_H
