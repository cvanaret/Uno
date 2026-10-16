// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NEWTONSYSTEMSYMBOLS_H
#define UNO_NEWTONSYSTEMSYMBOLS_H

#include <cstddef>
#include <tuple>
#include <type_traits>
#include <utility>
#include "linear_algebra/DiagonalMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/View.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/Sum.hpp"
#include "symbolic/Transpose.hpp"
#include "symbolic/symbolic_traits.hpp"
#include "../interfaces/C/uno_int.h"

namespace uno {
   // forward declarations
   class Evaluations;
   class HessianModel;
   class Statistics;

   // what the symbols of a Newton system may use when they are evaluated at the current iterate
   struct NewtonSystemContext {
      Statistics& statistics;
      const Iterate& iterate;
      Evaluations& evaluations;
   };

   // Symbols of the Newton system. A symbol is a sparse matrix whose entries are stored in a contiguous slice of a COO
   // array (in the order given by compute_sparsity). It exposes:
   // - number_rows(), number_columns(), number_nonzeros()
   // - compute_sparsity(row_indices, column_indices, row_offset, column_offset, solver_indexing): writes its pattern (lower
   //   triangle for symmetric matrices) into views of size number_nonzeros()
   // - evaluate(context, values): writes its values into a view of size number_nonzeros()
   // - add_matrix_vector_product(context, x, y): y += A x (operator form, if available)

   // Lagrangian Hessian of an optimization problem, embedded in the top-left corner of a (possibly larger) square matrix
   class LagrangianHessian: public SymbolicExpression {
   public:
      using value_type = double;

      LagrangianHessian(const OptimizationProblem& problem, HessianModel& hessian_model, size_t dimension):
            problem(problem), hessian_model(hessian_model), dimension(dimension) { }

      [[nodiscard]] size_t number_rows() const { return this->dimension; }
      [[nodiscard]] size_t number_columns() const { return this->dimension; }
      [[nodiscard]] size_t number_nonzeros() const { return this->problem.number_hessian_nonzeros(this->hessian_model); }

      void compute_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, size_t row_offset, size_t column_offset,
            uno_int solver_indexing) const {
         this->problem.compute_hessian_sparsity(this->hessian_model, row_indices, column_indices, solver_indexing);
         for (size_t nonzero_index: Range(this->number_nonzeros())) {
            row_indices[nonzero_index] += static_cast<uno_int>(row_offset);
            column_indices[nonzero_index] += static_cast<uno_int>(column_offset);
         }
      }

      void evaluate(const NewtonSystemContext& context, View<double> values) const {
         this->problem.evaluate_lagrangian_hessian(context.statistics, this->hessian_model, context.iterate.primals,
            context.iterate.multipliers, values);
      }

      // y += ∇²L x (Hessian-vector product of the problem)
      void add_matrix_vector_product(const NewtonSystemContext& context, View<const double> x, View<double> y) const {
         const size_t number_variables = this->problem.number_variables;
         this->product.resize(number_variables);
         this->problem.compute_hessian_vector_product(this->hessian_model, view(context.iterate.primals.data(), number_variables),
            view(x.data(), number_variables), context.iterate.multipliers, this->product.view());
         for (size_t variable_index: Range(number_variables)) {
            y[variable_index] += this->product[variable_index];
         }
      }

   protected:
      const OptimizationProblem& problem;
      HessianModel& hessian_model;
      const size_t dimension;
      mutable Vector<double> product{};
   };

   // constraint Jacobian of an optimization problem
   class ConstraintJacobian: public SymbolicExpression {
   public:
      using value_type = double;

      explicit ConstraintJacobian(const OptimizationProblem& problem): problem(problem) { }

      [[nodiscard]] size_t number_rows() const { return this->problem.number_constraints; }
      [[nodiscard]] size_t number_columns() const { return this->problem.number_variables; }
      [[nodiscard]] size_t number_nonzeros() const { return this->problem.number_jacobian_nonzeros(); }

      void compute_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, size_t row_offset, size_t column_offset,
            uno_int solver_indexing) const {
         const auto jacobian_row_indices = this->problem.get_jacobian_row_indices();
         const auto jacobian_column_indices = this->problem.get_jacobian_column_indices();
         for (size_t nonzero_index: Range(this->number_nonzeros())) {
            row_indices[nonzero_index] = jacobian_row_indices[nonzero_index] + static_cast<uno_int>(row_offset) + solver_indexing;
            column_indices[nonzero_index] = jacobian_column_indices[nonzero_index] + static_cast<uno_int>(column_offset) + solver_indexing;
         }
      }

      void evaluate(const NewtonSystemContext& context, View<double> values) const {
         this->problem.evaluate_jacobian(context.iterate.primals, values, context.evaluations);
      }

      // y += J x
      void add_matrix_vector_product(const NewtonSystemContext& context, View<const double> x, View<double> y) const {
         this->problem.add_jacobian_vector_product(x, y, context.evaluations);
      }

      // y += J^T x
      void add_transposed_matrix_vector_product(const NewtonSystemContext& context, View<const double> x, View<double> y) const {
         this->problem.add_jacobian_transposed_vector_product(x, y, context.evaluations);
      }

   protected:
      const OptimizationProblem& problem;
   };

   // diagonal block whose values are set outside of the evaluation (by the inertia correction strategy)
   class Regularization: public SymbolicExpression {
   public:
      using value_type = double;

      explicit Regularization(size_t dimension): dimension(dimension) { }

      [[nodiscard]] size_t number_rows() const { return this->dimension; }
      [[nodiscard]] size_t number_columns() const { return this->dimension; }
      [[nodiscard]] size_t number_nonzeros() const { return this->dimension; }

      void compute_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, size_t row_offset, size_t column_offset,
            uno_int solver_indexing) const {
         for (size_t index: Range(this->dimension)) {
            row_indices[index] = static_cast<uno_int>(row_offset + index) + solver_indexing;
            column_indices[index] = static_cast<uno_int>(column_offset + index) + solver_indexing;
         }
      }

      void evaluate(const NewtonSystemContext& /*context*/, View<double> /*values*/) const {
         // the values are set by the inertia correction strategy
      }

      // the values are owned by the linear system (inertia correction): the operator is applied with δ = 0
      void add_matrix_vector_product(const NewtonSystemContext& /*context*/, View<const double> /*x*/, View<double> /*y*/) const {
      }

   protected:
      const size_t dimension;
   };

   // index map k -> k
   struct IdentityIndex {
      [[nodiscard]] constexpr size_t operator()(size_t index) const noexcept { return index; }
   };

   // index map k -> constant index
   struct ConstantIndex {
      size_t index;
      [[nodiscard]] constexpr size_t operator()(size_t /*entry*/) const noexcept { return this->index; }
   };

   // constant value
   struct ConstantValue {
      double value;
      [[nodiscard]] constexpr double operator()(size_t /*index*/) const noexcept { return this->value; }
   };

   // sparse matrix with at most one nonzero per row and per column, given as dematerialized entries: entry k (in
   // [0, number_entries)) is (row_index(k), column_index(k), value(k)). Diagonal matrices and (scaled) selection matrices
   // are particular cases. Its entries are computed when the matrix is evaluated or applied
   template <typename RowIndex, typename ColumnIndex, typename Value>
   class ScaledSelectionMatrix: public SymbolicExpression {
   public:
      using value_type = double;

      ScaledSelectionMatrix(size_t number_rows, size_t number_columns, size_t number_entries, RowIndex row_index,
            ColumnIndex column_index, Value value): rows(number_rows), columns(number_columns), entries(number_entries),
            row_index_function(std::move(row_index)), column_index_function(std::move(column_index)), value_function(std::move(value)) { }

      [[nodiscard]] size_t number_rows() const { return this->rows; }
      [[nodiscard]] size_t number_columns() const { return this->columns; }
      [[nodiscard]] size_t number_nonzeros() const { return this->entries; }
      [[nodiscard]] size_t row_index(size_t entry) const { return this->row_index_function(entry); }
      [[nodiscard]] size_t column_index(size_t entry) const { return this->column_index_function(entry); }
      [[nodiscard]] double value(size_t entry) const { return this->value_function(entry); }

      void compute_sparsity(View<uno_int> row_indices, View<uno_int> column_indices, size_t row_offset, size_t column_offset,
            uno_int solver_indexing) const {
         for (size_t entry: Range(this->entries)) {
            row_indices[entry] = static_cast<uno_int>(row_offset + this->row_index(entry)) + solver_indexing;
            column_indices[entry] = static_cast<uno_int>(column_offset + this->column_index(entry)) + solver_indexing;
         }
      }

      void evaluate(const NewtonSystemContext& /*context*/, View<double> values) const {
         for (size_t entry: Range(this->entries)) {
            values[entry] = this->value(entry);
         }
      }

      // y += A x
      void add_matrix_vector_product(const NewtonSystemContext& /*context*/, View<const double> x, View<double> y) const {
         for (size_t entry: Range(this->entries)) {
            y[this->row_index(entry)] += this->value(entry) * x[this->column_index(entry)];
         }
      }

      // y += A^T x
      void add_transposed_matrix_vector_product(const NewtonSystemContext& /*context*/, View<const double> x, View<double> y) const {
         for (size_t entry: Range(this->entries)) {
            y[this->column_index(entry)] += this->value(entry) * x[this->row_index(entry)];
         }
      }

   protected:
      const size_t rows;
      const size_t columns;
      const size_t entries;
      RowIndex row_index_function;
      ColumnIndex column_index_function;
      Value value_function;
   };

   template <typename T>
   struct is_scaled_selection: std::false_type { };
   template <typename RowIndex, typename ColumnIndex, typename Value>
   struct is_scaled_selection<ScaledSelectionMatrix<RowIndex, ColumnIndex, Value>>: std::true_type { };
   template <typename T>
   inline constexpr bool is_scaled_selection_v = is_scaled_selection<std::decay_t<T>>::value;

   // free functions
   template <typename RowIndex, typename ColumnIndex, typename Value>
   auto scaled_selection(size_t number_rows, size_t number_columns, size_t number_entries, RowIndex&& row_index,
         ColumnIndex&& column_index, Value&& value) {
      return ScaledSelectionMatrix<std::decay_t<RowIndex>, std::decay_t<ColumnIndex>, std::decay_t<Value>>(number_rows,
         number_columns, number_entries, std::forward<RowIndex>(row_index), std::forward<ColumnIndex>(column_index),
         std::forward<Value>(value));
   }

   // square diagonal matrix with entries value(k)
   template <typename Value>
   auto diagonal(size_t dimension, Value&& value) {
      return scaled_selection(dimension, dimension, dimension, IdentityIndex{}, IdentityIndex{}, std::forward<Value>(value));
   }

   // right-hand side components (vectors evaluated at the current iterate)
   // vector whose entry k is value(k)
   template <typename Value>
   class ComponentVector {
   public:
      ComponentVector(size_t dimension, Value value): dimension(dimension), value_function(std::move(value)) { }

      [[nodiscard]] size_t size() const { return this->dimension; }
      [[nodiscard]] double entry(size_t index) const { return this->value_function(index); }

      void evaluate(const NewtonSystemContext& /*context*/, View<double> result) const {
         for (size_t index: Range(this->dimension)) {
            result[index] = this->entry(index);
         }
      }

   protected:
      const size_t dimension;
      Value value_function;
   };

   template <typename T>
   struct is_component_vector: std::false_type { };
   template <typename Value>
   struct is_component_vector<ComponentVector<Value>>: std::true_type { };
   template <typename T>
   inline constexpr bool is_component_vector_v = is_component_vector<std::decay_t<T>>::value;

   template <typename Value>
   ComponentVector<std::decay_t<Value>> component_vector(size_t dimension, Value&& value) {
      return {dimension, std::forward<Value>(value)};
   }

   // -g(x) + J(x)^T y, where g is the gradient of the objective (given as a function). The terms of the bound duals in
   // -∇_x L(x, y, z) are linear terms of the Newton right-hand side (see NewtonRHS)
   template <typename ObjectiveGradient>
   class NegativeLagrangianGradient {
   public:
      NegativeLagrangianGradient(const OptimizationProblem& problem, ObjectiveGradient objective_gradient):
            problem(problem), objective_gradient(std::move(objective_gradient)) { }

      [[nodiscard]] size_t size() const { return this->problem.number_variables; }

      void evaluate(const NewtonSystemContext& context, View<double> result) const {
         const Iterate& iterate = context.iterate;
         this->objective_gradient(iterate, result, context.evaluations);
         result.scale(-1.);
         this->problem.add_jacobian_transposed_vector_product(iterate.multipliers.constraints.view(), result, context.evaluations);
      }

   protected:
      const OptimizationProblem& problem;
      ObjectiveGradient objective_gradient;
   };

   template <typename ObjectiveGradient>
   NegativeLagrangianGradient<std::decay_t<ObjectiveGradient>> negative_lagrangian_gradient(const OptimizationProblem& problem,
         ObjectiveGradient&& objective_gradient) {
      return {problem, std::forward<ObjectiveGradient>(objective_gradient)};
   }

   // Newton right-hand side of a block row in residual form: the row residual is F = F0 + sum_j K(row, j) v_j, linear in
   // the current values v_j of some blocks of unknowns, with the blocks K(row, j) of the Newton matrix as coefficients.
   // The right-hand side is -F = constant - sum_j K(row, j) v_j with constant = -F0.
   // Declaring the linear terms lets the condensation cancel them symbolically instead of numerically

   // linear term K(row, Column) v_Column, where value(k) is the current value of the k-th unknown of block Column
   template <size_t Column, typename Value>
   class LinearTerm {
   public:
      static constexpr size_t column = Column;

      explicit LinearTerm(Value value): value_function(std::move(value)) { }

      [[nodiscard]] double value(size_t index) const { return this->value_function(index); }

   protected:
      Value value_function;
   };

   template <size_t Column, typename Value>
   LinearTerm<Column, std::decay_t<Value>> linear_term(Value&& value) {
      return LinearTerm<Column, std::decay_t<Value>>(std::forward<Value>(value));
   }

   template <typename Constant, typename... Terms>
   class NewtonRHS {
   public:
      explicit NewtonRHS(Constant constant, Terms... terms): constant(std::move(constant)), terms(std::move(terms)...) { }

      [[nodiscard]] size_t size() const { return this->constant.size(); }
      [[nodiscard]] const Constant& get_constant() const { return this->constant; }

      template <typename Function>
      void for_each_term(Function&& function) const {
         std::apply([&](const auto&... term) { (function(term), ...); }, this->terms);
      }

   protected:
      Constant constant;
      std::tuple<Terms...> terms;
   };

   template <typename Constant, typename... Terms>
   NewtonRHS<std::decay_t<Constant>, std::decay_t<Terms>...> newton_rhs(Constant&& constant, Terms&&... terms) {
      return NewtonRHS<std::decay_t<Constant>, std::decay_t<Terms>...>(std::forward<Constant>(constant), std::forward<Terms>(terms)...);
   }

   template <typename T>
   struct is_newton_rhs: std::false_type { };
   template <typename Constant, typename... Terms>
   struct is_newton_rhs<NewtonRHS<Constant, Terms...>>: std::true_type { };
   template <typename T>
   inline constexpr bool is_newton_rhs_v = is_newton_rhs<std::decay_t<T>>::value;

   // a plain right-hand side component (a vector) is a NewtonRHS without linear terms
   template <typename Component>
   const auto& rhs_constant(const Component& component) {
      if constexpr (is_newton_rhs_v<Component>) {
         return component.get_constant();
      }
      else {
         return component;
      }
   }

   template <typename Component, typename Function>
   void for_each_linear_term(const Component& component, Function&& function) {
      if constexpr (is_newton_rhs_v<Component>) {
         component.for_each_term(std::forward<Function>(function));
      }
   }

   // -(c(x) - c_L)
   class NegativeConstraintResidual {
   public:
      explicit NegativeConstraintResidual(const OptimizationProblem& problem): problem(problem) { }

      [[nodiscard]] size_t size() const { return this->problem.number_constraints; }

      void evaluate(const NewtonSystemContext& context, View<double> result) const {
         this->problem.evaluate_constraints(context.iterate, result, context.evaluations);
         const auto& constraints_lower_bounds = this->problem.get_constraints_lower_bounds();
         for (size_t constraint_index: Range(this->problem.number_constraints)) {
            result[constraint_index] = -(result[constraint_index] - constraints_lower_bounds[constraint_index]);
         }
      }

   protected:
      const OptimizationProblem& problem;
   };

   // symmetric 2x2 block matrix [A11 *; A21 A22]: the (1, 2) block is the transpose of A21 and is not stored
   template <typename A11, typename A21, typename A22>
   class SymmetricBlockMatrix: public SymbolicExpression {
   public:
      using value_type = double;

      SymmetricBlockMatrix(A11&& a11, A21&& a21, A22&& a22): a11(std::forward<A11>(a11)), a21(std::forward<A21>(a21)),
            a22(std::forward<A22>(a22)) { }

      UNO_FORWARD_ACCESSOR(get_a11, this->a11)
      UNO_FORWARD_ACCESSOR(get_a21, this->a21)
      UNO_FORWARD_ACCESSOR(get_a22, this->a22)

   protected:
      storage_t<A11> a11;
      storage_t<A21> a21;
      storage_t<A22> a22;
   };

   template <typename A11, typename A21, typename A22>
   SymmetricBlockMatrix<A11, A21, A22> symmetric_block_matrix(A11&& a11, A21&& a21, A22&& a22) {
      return {std::forward<A11>(a11), std::forward<A21>(a21), std::forward<A22>(a22)};
   }

   //*****************************************************************
   // operations dispatched on the type of the expression
   //*****************************************************************

   // symbols (member functions)
   template <typename Symbol>
   size_t number_rows(const Symbol& symbol) { return symbol.number_rows(); }

   template <typename Symbol>
   size_t number_nonzeros(const Symbol& symbol) { return symbol.number_nonzeros(); }

   template <typename Symbol>
   void compute_sparsity(const Symbol& symbol, View<uno_int> row_indices, View<uno_int> column_indices, size_t row_offset,
         size_t column_offset, uno_int solver_indexing) {
      symbol.compute_sparsity(row_indices, column_indices, row_offset, column_offset, solver_indexing);
   }

   template <typename Symbol>
   void evaluate(const Symbol& symbol, const NewtonSystemContext& context, View<double> values) {
      symbol.evaluate(context, values);
   }

   template <typename Symbol>
   void add_matrix_vector_product(const Symbol& symbol, const NewtonSystemContext& context, View<const double> x, View<double> y) {
      symbol.add_matrix_vector_product(context, x, y);
   }

   template <typename Symbol>
   void add_transposed_matrix_vector_product(const Symbol& symbol, const NewtonSystemContext& context, View<const double> x,
         View<double> y) {
      symbol.add_transposed_matrix_vector_product(context, x, y);
   }

   // transpose: y += A^T x
   template <typename Matrix>
   void add_matrix_vector_product(const Transpose<Matrix>& transposed, const NewtonSystemContext& context,
         View<const double> x, View<double> y) {
      add_transposed_matrix_vector_product(transposed.get_matrix(), context, x, y);
   }

   // diagonal matrix: one entry per index
   template <typename Indices, typename ComponentFunction>
   size_t number_nonzeros(const DiagonalMatrix<Indices, ComponentFunction>& matrix) {
      return matrix.get_indices().size();
   }

   template <typename Indices, typename ComponentFunction>
   void compute_sparsity(const DiagonalMatrix<Indices, ComponentFunction>& matrix, View<uno_int> row_indices,
         View<uno_int> column_indices, size_t row_offset, size_t column_offset, uno_int solver_indexing) {
      size_t nonzero_index = 0;
      for (const size_t index: matrix.get_indices()) {
         row_indices[nonzero_index] = static_cast<uno_int>(row_offset + index) + solver_indexing;
         column_indices[nonzero_index] = static_cast<uno_int>(column_offset + index) + solver_indexing;
         ++nonzero_index;
      }
   }

   template <typename Indices, typename ComponentFunction>
   void evaluate(const DiagonalMatrix<Indices, ComponentFunction>& matrix, const NewtonSystemContext& /*context*/,
         View<double> values) {
      size_t nonzero_index = 0;
      for (const size_t index: matrix.get_indices()) {
         values[nonzero_index] = matrix.entry(index);
         ++nonzero_index;
      }
   }

   template <typename Indices, typename ComponentFunction>
   void add_matrix_vector_product(const DiagonalMatrix<Indices, ComponentFunction>& matrix, const NewtonSystemContext& /*context*/,
         View<const double> x, View<double> y) {
      for (const size_t index: matrix.get_indices()) {
         y[index] += matrix.entry(index) * x[index];
      }
   }

   // sum: the entries of both operands are stored one after the other (duplicate positions are summed by the solvers)
   template <typename L, typename R>
   size_t number_rows(const Sum<L, R>& sum) { return number_rows(sum.get_left()); }

   template <typename L, typename R>
   size_t number_nonzeros(const Sum<L, R>& sum) {
      return number_nonzeros(sum.get_left()) + number_nonzeros(sum.get_right());
   }

   template <typename L, typename R>
   void compute_sparsity(const Sum<L, R>& sum, View<uno_int> row_indices, View<uno_int> column_indices, size_t row_offset,
         size_t column_offset, uno_int solver_indexing) {
      const size_t left_nonzeros = number_nonzeros(sum.get_left());
      const size_t right_nonzeros = number_nonzeros(sum.get_right());
      compute_sparsity(sum.get_left(), view(row_indices.data(), left_nonzeros), view(column_indices.data(), left_nonzeros),
         row_offset, column_offset, solver_indexing);
      compute_sparsity(sum.get_right(), view(row_indices.data() + left_nonzeros, right_nonzeros),
         view(column_indices.data() + left_nonzeros, right_nonzeros), row_offset, column_offset, solver_indexing);
   }

   template <typename L, typename R>
   void evaluate(const Sum<L, R>& sum, const NewtonSystemContext& context, View<double> values) {
      const size_t left_nonzeros = number_nonzeros(sum.get_left());
      evaluate(sum.get_left(), context, view(values.data(), left_nonzeros));
      evaluate(sum.get_right(), context, view(values.data() + left_nonzeros, number_nonzeros(sum.get_right())));
   }

   template <typename L, typename R>
   void add_matrix_vector_product(const Sum<L, R>& sum, const NewtonSystemContext& context, View<const double> x, View<double> y) {
      add_matrix_vector_product(sum.get_left(), context, x, y);
      add_matrix_vector_product(sum.get_right(), context, x, y);
   }

   // symmetric block matrix [A11 *; A21 A22]: entries of A11, then A21, then A22
   // dimension: rows of A11 + rows of A21 (A22 may be empty, e.g. without dual regularization)
   template <typename A11, typename A21, typename A22>
   size_t number_rows(const SymmetricBlockMatrix<A11, A21, A22>& matrix) {
      return number_rows(matrix.get_a11()) + number_rows(matrix.get_a21());
   }

   template <typename A11, typename A21, typename A22>
   size_t number_nonzeros(const SymmetricBlockMatrix<A11, A21, A22>& matrix) {
      return number_nonzeros(matrix.get_a11()) + number_nonzeros(matrix.get_a21()) + number_nonzeros(matrix.get_a22());
   }

   template <typename A11, typename A21, typename A22>
   void compute_sparsity(const SymmetricBlockMatrix<A11, A21, A22>& matrix, View<uno_int> row_indices, View<uno_int> column_indices,
         size_t row_offset, size_t column_offset, uno_int solver_indexing) {
      const size_t dimension1 = number_rows(matrix.get_a11());
      size_t position = 0;
      const auto block = [&](const auto& a, size_t block_row_offset, size_t block_column_offset) {
         const size_t block_nonzeros = number_nonzeros(a);
         compute_sparsity(a, view(row_indices.data() + position, block_nonzeros), view(column_indices.data() + position, block_nonzeros),
            row_offset + block_row_offset, column_offset + block_column_offset, solver_indexing);
         position += block_nonzeros;
      };
      block(matrix.get_a11(), 0, 0);
      block(matrix.get_a21(), dimension1, 0);
      block(matrix.get_a22(), dimension1, dimension1);
   }

   template <typename A11, typename A21, typename A22>
   void evaluate(const SymmetricBlockMatrix<A11, A21, A22>& matrix, const NewtonSystemContext& context, View<double> values) {
      size_t position = 0;
      const auto block = [&](const auto& a) {
         const size_t block_nonzeros = number_nonzeros(a);
         evaluate(a, context, view(values.data() + position, block_nonzeros));
         position += block_nonzeros;
      };
      block(matrix.get_a11());
      block(matrix.get_a21());
      block(matrix.get_a22());
   }
} // namespace

#endif // UNO_NEWTONSYSTEMSYMBOLS_H
