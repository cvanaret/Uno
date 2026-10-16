// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_NEWTONSYSTEMSYMBOLS_H
#define UNO_NEWTONSYSTEMSYMBOLS_H

#include <cstddef>
#include <utility>
#include "linear_algebra/DiagonalMatrix.hpp"
#include "linear_algebra/Vector.hpp"
#include "linear_algebra/View.hpp"
#include "optimization/Iterate.hpp"
#include "optimization/OptimizationProblem.hpp"
#include "symbolic/Range.hpp"
#include "symbolic/Sum.hpp"
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

   protected:
      const size_t dimension;
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
