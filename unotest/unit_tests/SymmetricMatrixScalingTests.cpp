// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>
#include <gtest/gtest.h>
#include "ingredients/subproblem_solvers/SymmetricMatrixScaling.hpp"

using namespace uno;

namespace {
   struct COOMatrix {
      size_t dimension;
      std::vector<uno_int> rows{};
      std::vector<uno_int> columns{};
      std::vector<double> values{};

      void add(size_t row, size_t column, double value) {
         this->rows.push_back(static_cast<uno_int>(row));
         this->columns.push_back(static_cast<uno_int>(column));
         this->values.push_back(value);
      }

      // dense symmetric matrix (duplicates summed, one triangle stored)
      [[nodiscard]] std::vector<std::vector<double>> dense() const {
         std::vector<std::vector<double>> matrix(this->dimension, std::vector<double>(this->dimension, 0.));
         for (size_t index = 0; index < this->values.size(); ++index) {
            const auto row = static_cast<size_t>(this->rows[index]), column = static_cast<size_t>(this->columns[index]);
            matrix[row][column] += this->values[index];
            if (row != column) {
               matrix[column][row] += this->values[index];
            }
         }
         return matrix;
      }
   };

   std::vector<double> compute_scaling(const COOMatrix& matrix, LinearSystemScalingMethod method, SymmetricMatrixScaling& scaling_method,
         int indexing = 0) {
      std::vector<uno_int> rows = matrix.rows, columns = matrix.columns;
      for (size_t index = 0; index < rows.size(); ++index) {
         rows[index] += indexing;
         columns[index] += indexing;
      }
      scaling_method.initialize(matrix.dimension, rows.data(), columns.data(), matrix.values.size(), indexing);
      std::vector<double> scaling;
      scaling_method.compute(method, matrix.values.data(), scaling);
      return scaling;
   }

   std::vector<double> compute_scaling(const COOMatrix& matrix, LinearSystemScalingMethod method) {
      SymmetricMatrixScaling scaling_method;
      return compute_scaling(matrix, method, scaling_method);
   }

   double max_scaled_entry(const COOMatrix& matrix, const std::vector<double>& scaling) {
      const auto dense = matrix.dense();
      double maximum = 0.;
      for (size_t i = 0; i < matrix.dimension; ++i) {
         for (size_t j = 0; j < matrix.dimension; ++j) {
            maximum = std::max(maximum, std::abs(scaling[i] * dense[i][j] * scaling[j]));
         }
      }
      return maximum;
   }

   // ratio between the largest and smallest infinity norms of the rows of S A S
   double row_norm_spread(const COOMatrix& matrix, const std::vector<double>& scaling) {
      const auto dense = matrix.dense();
      double largest = 0., smallest = std::numeric_limits<double>::infinity();
      for (size_t i = 0; i < matrix.dimension; ++i) {
         double row_norm = 0.;
         for (size_t j = 0; j < matrix.dimension; ++j) {
            row_norm = std::max(row_norm, std::abs(scaling[i] * dense[i][j] * scaling[j]));
         }
         largest = std::max(largest, row_norm);
         smallest = std::min(smallest, row_norm);
      }
      return largest / smallest;
   }

   // badly scaled KKT matrix [H J^T; J -delta I]: entries spread over 16 orders of magnitude
   COOMatrix badly_scaled_kkt() {
      COOMatrix matrix{5};
      matrix.add(0, 0, 1e8);
      matrix.add(1, 1, 2e-6);
      matrix.add(2, 2, 3.);
      matrix.add(1, 0, 1e2);
      matrix.add(3, 0, 1e4); // J
      matrix.add(3, 1, 1e-3);
      matrix.add(4, 1, 5e-4);
      matrix.add(4, 2, 7e3);
      matrix.add(3, 3, 0.); // dual regularization (zero)
      matrix.add(4, 4, 0.);
      return matrix;
   }

   std::vector<double> unit_scaling(size_t dimension) {
      return std::vector<double>(dimension, 1.);
   }
} // namespace

TEST(SymmetricMatrixScaling, MethodFromString) {
   EXPECT_EQ(linear_system_scaling_method_from_string("none"), LinearSystemScalingMethod::NONE);
   EXPECT_EQ(linear_system_scaling_method_from_string("ruiz"), LinearSystemScalingMethod::RUIZ);
   EXPECT_THROW((void)linear_system_scaling_method_from_string("mc64"), std::invalid_argument);
}

TEST(SymmetricMatrixScaling, ComputeRequiresInitialization) {
   SymmetricMatrixScaling scaling_method;
   std::vector<double> scaling;
   const double values[1] = {1.};
   EXPECT_THROW(scaling_method.compute(LinearSystemScalingMethod::RUIZ, values, scaling), std::logic_error);
}

// diagonal matrix: both methods scale the diagonal to +-1
TEST(SymmetricMatrixScaling, DiagonalMatrix) {
   COOMatrix matrix{3};
   matrix.add(0, 0, 4.);
   matrix.add(1, 1, -1e-6);
   matrix.add(2, 2, 1e10);
   const std::vector<double> scaling = compute_scaling(matrix, LinearSystemScalingMethod::RUIZ);
   EXPECT_NEAR(scaling[0], 0.5, 1e-12);
   EXPECT_NEAR(scaling[1] / 1e3, 1., 1e-12);
   EXPECT_NEAR(scaling[2] / 1e-5, 1., 1e-12);
}

TEST(SymmetricMatrixScaling, NoneIsIdentity) {
   const COOMatrix matrix = badly_scaled_kkt();
   EXPECT_EQ(compute_scaling(matrix, LinearSystemScalingMethod::NONE), unit_scaling(matrix.dimension));
}

// Ruiz: the row norms of S A S are close to 1
TEST(SymmetricMatrixScaling, RuizEquilibratesRows) {
   const COOMatrix matrix = badly_scaled_kkt();
   const std::vector<double> scaling = compute_scaling(matrix, LinearSystemScalingMethod::RUIZ);
   EXPECT_LE(max_scaled_entry(matrix, scaling), 1. + 0.1);
   EXPECT_LT(row_norm_spread(matrix, scaling), 1.5);
}

// the duplicates are summed, and either triangle may be given
TEST(SymmetricMatrixScaling, DuplicatesAndTriangles) {
   const COOMatrix reference = badly_scaled_kkt();
   COOMatrix split{reference.dimension};
   for (size_t index = 0; index < reference.values.size(); ++index) {
      const auto row = static_cast<size_t>(reference.rows[index]), column = static_cast<size_t>(reference.columns[index]);
      // upper triangle, value split into two duplicates
      split.add(column, row, 0.25 * reference.values[index]);
      split.add(row, column, 0.75 * reference.values[index]);
   }
   const std::vector<double> expected = compute_scaling(reference, LinearSystemScalingMethod::RUIZ);
   const std::vector<double> scaling = compute_scaling(split, LinearSystemScalingMethod::RUIZ);
   for (size_t index = 0; index < reference.dimension; ++index) {
      EXPECT_NEAR(scaling[index] / expected[index], 1., 1e-12);
   }
}

TEST(SymmetricMatrixScaling, FortranIndexing) {
   const COOMatrix matrix = badly_scaled_kkt();
   SymmetricMatrixScaling c_scaling, fortran_scaling;
   EXPECT_EQ(compute_scaling(matrix, LinearSystemScalingMethod::RUIZ, c_scaling, 0),
             compute_scaling(matrix, LinearSystemScalingMethod::RUIZ, fortran_scaling, 1));
}

TEST(SymmetricMatrixScaling, IndexOutOfRange) {
   SymmetricMatrixScaling scaling_method;
   const uno_int rows[1] = {2}, columns[1] = {0};
   EXPECT_THROW(scaling_method.initialize(2, rows, columns, 1, 0), std::out_of_range);
}

TEST(SymmetricMatrixScaling, NonFiniteValuesLeaveMatrixUnscaled) {
   COOMatrix matrix = badly_scaled_kkt();
   matrix.values[2] = std::numeric_limits<double>::quiet_NaN();
   EXPECT_EQ(compute_scaling(matrix, LinearSystemScalingMethod::RUIZ), unit_scaling(matrix.dimension));
}