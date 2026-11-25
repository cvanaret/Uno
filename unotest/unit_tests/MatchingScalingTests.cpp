// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

// symmetric matching-based scaling (Duff & Koster 2001, MC64 job 5; Duff & Pralet 2005): S A S has all its entries at
// most 1 in absolute value, and a (maximum product) matching of entries equal to 1 in absolute value. These properties
// hold for any optimal dual solution of the matching problem, so they do not depend on the implementation

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <random>
#include <vector>
#include <gtest/gtest.h>
#include "ingredients/subproblem_solvers/SymmetricMatrixScaling.hpp"
#include "options/DefaultOptions.hpp"
#include "options/Options.hpp"

using namespace uno;

namespace {
   struct SymmetricMatrix {
      size_t dimension;
      std::vector<uno_int> rows{};
      std::vector<uno_int> columns{};
      std::vector<double> values{};

      void add(size_t row, size_t column, double value) {
         this->rows.push_back(static_cast<uno_int>(row));
         this->columns.push_back(static_cast<uno_int>(column));
         this->values.push_back(value);
      }

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

   std::vector<double> matching_scaling(const SymmetricMatrix& matrix) {
      SymmetricMatrixScaling scaling_method;
      scaling_method.initialize(matrix.dimension, matrix.rows.data(), matrix.columns.data(), matrix.values.size(), 0);
      std::vector<double> scaling;
      scaling_method.compute(linear_system_scaling_method_from_string("matching"), matrix.values.data(), scaling);
      return scaling;
   }

   // S A S: all entries at most 1, and a perfect matching of entries equal to 1 (checked with a bipartite matching of the
   // unit entries, by augmenting paths)
   void expect_matching_scaling_properties(const SymmetricMatrix& matrix, const std::vector<double>& scaling) {
      ASSERT_EQ(scaling.size(), matrix.dimension);
      const auto dense = matrix.dense();
      const size_t n = matrix.dimension;
      std::vector<std::vector<size_t>> unit_entries(n);
      for (size_t i = 0; i < n; ++i) {
         EXPECT_TRUE(std::isfinite(scaling[i]) && 0. < scaling[i]) << "scaling " << i;
         for (size_t j = 0; j < n; ++j) {
            const double entry = std::abs(scaling[i] * dense[i][j] * scaling[j]);
            EXPECT_LE(entry, 1. + 1e-10) << "entry (" << i << ", " << j << ")";
            if (dense[i][j] != 0. && 1. - 1e-10 <= entry) {
               unit_entries[i].push_back(j);
            }
         }
      }
      std::vector<size_t> column_match(n, n);
      std::vector<bool> visited;
      const std::function<bool(size_t)> augment = [&](size_t row) {
         for (size_t column: unit_entries[row]) {
            if (!visited[column]) {
               visited[column] = true;
               if (column_match[column] == n || augment(column_match[column])) {
                  column_match[column] = row;
                  return true;
               }
            }
         }
         return false;
      };
      size_t matching_size = 0;
      for (size_t row = 0; row < n; ++row) {
         visited.assign(n, false);
         matching_size += augment(row) ? 1 : 0;
      }
      EXPECT_EQ(matching_size, n) << "the unit entries of S A S do not contain a perfect matching";
   }

   // badly scaled KKT matrix [H J^T; J 0] with explicit zeros in the (2, 2) block
   SymmetricMatrix badly_scaled_kkt() {
      SymmetricMatrix matrix{5};
      matrix.add(0, 0, 1e8);
      matrix.add(1, 1, 2e-6);
      matrix.add(2, 2, 3.);
      matrix.add(1, 0, 1e2);
      matrix.add(3, 0, 1e4);
      matrix.add(3, 1, 1e-3);
      matrix.add(4, 1, 5e-4);
      matrix.add(4, 2, 7e3);
      matrix.add(3, 3, 0.);
      matrix.add(4, 4, 0.);
      return matrix;
   }

   // properties of the scaling of any matrix (Duff & Pralet 2005 for structurally singular matrices): finite positive
   // scaling, |s_i a_ij s_j| <= 1, and the largest entry of each nonempty row is 1
   void expect_scaling_properties(const SymmetricMatrix& matrix, const std::vector<double>& scaling) {
      ASSERT_EQ(scaling.size(), matrix.dimension);
      const auto dense = matrix.dense();
      for (size_t i = 0; i < matrix.dimension; ++i) {
         EXPECT_TRUE(std::isfinite(scaling[i]) && 0. < scaling[i]);
         double row_maximum = 0.;
         bool empty_row = true;
         for (size_t j = 0; j < matrix.dimension; ++j) {
            const double entry = std::abs(scaling[i] * dense[i][j] * scaling[j]);
            EXPECT_LE(entry, 1. + 1e-10);
            row_maximum = std::max(row_maximum, entry);
            empty_row = empty_row && dense[i][j] == 0.;
         }
         if (!empty_row) {
            EXPECT_NEAR(row_maximum, 1., 1e-10) << "row " << i;
         }
      }
   }
} // namespace

TEST(MatchingScaling, MethodFromString) {
   EXPECT_NO_THROW((void)linear_system_scaling_method_from_string("matching"));
}

// diagonal matrix: s_i = 1 / sqrt(|d_i|)
TEST(MatchingScaling, DiagonalMatrix) {
   SymmetricMatrix matrix{3};
   matrix.add(0, 0, 4.);
   matrix.add(1, 1, -1e-6);
   matrix.add(2, 2, 1e10);
   const std::vector<double> scaling = matching_scaling(matrix);
   ASSERT_EQ(scaling.size(), 3u);
   EXPECT_NEAR(scaling[0], 0.5, 1e-12);
   EXPECT_NEAR(scaling[1] / 1e3, 1., 1e-12);
   EXPECT_NEAR(scaling[2] / 1e-5, 1., 1e-12);
}

// [eps a; a 0]: the matching is the off-diagonal pair, which is scaled to 1 (the tiny diagonal entry to at most 1)
TEST(MatchingScaling, SaddlePointBlock) {
   SymmetricMatrix matrix{2};
   matrix.add(0, 0, 1e-9);
   matrix.add(1, 0, 1e3);
   matrix.add(1, 1, 0.);
   const std::vector<double> scaling = matching_scaling(matrix);
   expect_matching_scaling_properties(matrix, scaling);
   EXPECT_NEAR(scaling[0] * 1e3 * scaling[1], 1., 1e-12);
}

TEST(MatchingScaling, BadlyScaledKKTMatrix) {
   const SymmetricMatrix matrix = badly_scaled_kkt();
   expect_matching_scaling_properties(matrix, matching_scaling(matrix));
}

// the duplicates are summed and either triangle may be given
TEST(MatchingScaling, DuplicatesAndTriangles) {
   const SymmetricMatrix reference = badly_scaled_kkt();
   SymmetricMatrix split{reference.dimension};
   for (size_t index = 0; index < reference.values.size(); ++index) {
      const auto row = static_cast<size_t>(reference.rows[index]), column = static_cast<size_t>(reference.columns[index]);
      split.add(column, row, 0.25 * reference.values[index]);
      split.add(row, column, 0.75 * reference.values[index]);
   }
   expect_matching_scaling_properties(split, matching_scaling(split));
}

// random sparse symmetric matrices (entries over 20 orders of magnitude), structurally nonsingular (nonzero diagonal, or
// a saddle point structure [H J^T; J 0] with J of full structural rank)
TEST(MatchingScaling, RandomMatrices) {
   std::mt19937 generator(42);
   std::uniform_real_distribution<double> exponent(-10., 10.);
   std::uniform_real_distribution<double> uniform(0., 1.);
   for (size_t trial = 0; trial < 50; ++trial) {
      const size_t n = 3 + trial % 20, m = trial % 2 == 0 ? 0 : 1 + trial % 5;
      SymmetricMatrix matrix{n + m};
      const auto value = [&]() { return (uniform(generator) < 0.5 ? -1. : 1.) * std::pow(10., exponent(generator)); };
      for (size_t i = 0; i < n; ++i) {
         matrix.add(i, i, value());
         for (size_t j = 0; j < i; ++j) {
            if (uniform(generator) < 0.3) {
               matrix.add(i, j, value());
            }
         }
      }
      for (size_t k = 0; k < m; ++k) {
         matrix.add(n + k, k, value()); // structural full rank of J
         for (size_t j = 0; j < n; ++j) {
            if (j != k && uniform(generator) < 0.3) {
               matrix.add(n + k, j, value());
            }
         }
      }
      expect_matching_scaling_properties(matrix, matching_scaling(matrix));
   }
}

// sequence of matrices with the same pattern and different values (the column duals of a matrix are the initial duals
// of the next one): the scaling of each matrix has the properties
TEST(MatchingScaling, SequenceWithSamePattern) {
   std::mt19937 generator(3);
   std::uniform_real_distribution<double> exponent(-12., 12.);
   SymmetricMatrix matrix = badly_scaled_kkt();
   for (size_t index = 0; index < 4; ++index) {
      matrix.add(index + 1, index, 1.); // denser pattern
   }
   SymmetricMatrixScaling scaling_method;
   scaling_method.initialize(matrix.dimension, matrix.rows.data(), matrix.columns.data(), matrix.values.size(), 0);
   for (size_t trial = 0; trial < 20; ++trial) {
      for (double& value: matrix.values) {
         if (value != 0.) {
            value = std::pow(10., exponent(generator));
         }
      }
      std::vector<double> scaling;
      scaling_method.compute(linear_system_scaling_method_from_string("matching"), matrix.values.data(), scaling);
      expect_matching_scaling_properties(matrix, scaling);
   }
}

// larger saddle point matrices [0 J^T; J 0] with a zero diagonal and a square J that contains a random permutation
// (structurally nonsingular): the cheap initial phases cannot match everything, the shortest augmenting paths do
TEST(MatchingScaling, LargeSaddlePointMatrices) {
   std::mt19937 generator(5);
   std::uniform_real_distribution<double> exponent(-10., 10.);
   std::uniform_real_distribution<double> uniform(0., 1.);
   for (size_t trial = 0; trial < 5; ++trial) {
      const size_t n = 100 + 20 * trial;
      std::vector<size_t> permutation(n);
      std::iota(permutation.begin(), permutation.end(), size_t(0));
      std::shuffle(permutation.begin(), permutation.end(), generator);
      SymmetricMatrix matrix{2 * n};
      for (size_t k = 0; k < n; ++k) {
         matrix.add(n + k, permutation[k], std::pow(10., exponent(generator)));
         for (size_t j = 0; j < n; ++j) {
            if (j != permutation[k] && uniform(generator) < 0.03) {
               matrix.add(n + k, j, -std::pow(10., exponent(generator)));
            }
         }
      }
      expect_matching_scaling_properties(matrix, matching_scaling(matrix));
   }
}

// long sequence of matrices with the same pattern whose entries randomly vanish (the warm start then comes from a
// matrix with other nonzeros), including structurally singular matrices
TEST(MatchingScaling, SequenceWithVanishingEntries) {
   std::mt19937 generator(8);
   std::uniform_real_distribution<double> exponent(-12., 12.);
   std::uniform_real_distribution<double> uniform(0., 1.);
   const size_t n = 30;
   SymmetricMatrix matrix{n};
   for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j <= i; ++j) {
         if (i == j || uniform(generator) < 0.15) {
            matrix.add(i, j, 1.);
         }
      }
   }
   SymmetricMatrixScaling scaling_method;
   scaling_method.initialize(matrix.dimension, matrix.rows.data(), matrix.columns.data(), matrix.values.size(), 0);
   for (size_t trial = 0; trial < 200; ++trial) {
      const double vanishing_probability = (trial % 4 == 3) ? 0.6 : 0.1;
      for (size_t index = 0; index < matrix.values.size(); ++index) {
         const bool diagonal = matrix.rows[index] == matrix.columns[index];
         // even trials: nonzero diagonal (structurally nonsingular)
         const bool vanishes = !(diagonal && trial % 2 == 0) && uniform(generator) < vanishing_probability;
         matrix.values[index] = vanishes ? 0. : (uniform(generator) < 0.5 ? -1. : 1.) * std::pow(10., exponent(generator));
      }
      std::vector<double> scaling;
      scaling_method.compute(linear_system_scaling_method_from_string("matching"), matrix.values.data(), scaling);
      expect_scaling_properties(matrix, scaling);
      if (trial % 2 == 0) {
         expect_matching_scaling_properties(matrix, scaling);
      }
   }
}

// matrices without entries: no scaling
TEST(MatchingScaling, EmptyMatrices) {
   EXPECT_TRUE(matching_scaling(SymmetricMatrix{0}).empty());
   SymmetricMatrix zeros{3};
   zeros.add(1, 0, 0.);
   EXPECT_EQ(matching_scaling(zeros), std::vector<double>(3, 1.));
}

// structurally singular matrices (Duff & Pralet 2005): the scaling is finite and positive, the entries of S A S are at
// most 1, and the largest entry of each nonempty row is 1
TEST(MatchingScaling, StructurallySingularMatrices) {
   const auto expect_properties = [](const SymmetricMatrix& matrix) {
      expect_scaling_properties(matrix, matching_scaling(matrix));
   };
   // indices 1 and 2 only couple to index 0
   SymmetricMatrix star{3};
   star.add(1, 0, 1e4);
   star.add(2, 0, 1e-4);
   expect_properties(star);
   // saddle point matrix [H J^T; J 0] with a rank-deficient J (two identical constraint rows), H of rank 2
   SymmetricMatrix saddle{5};
   saddle.add(0, 0, 1e-6);
   saddle.add(1, 1, 3e5);
   saddle.add(3, 0, 1e2);
   saddle.add(3, 2, 2.);
   saddle.add(4, 0, 1e-3);
   saddle.add(4, 2, 7e3);
   expect_properties(saddle);
   // random sparse symmetric matrices with zero diagonals (mostly structurally singular)
   std::mt19937 generator(11);
   std::uniform_real_distribution<double> exponent(-8., 8.);
   std::uniform_real_distribution<double> uniform(0., 1.);
   for (size_t trial = 0; trial < 100; ++trial) {
      const size_t n = 2 + trial % 15;
      SymmetricMatrix matrix{n};
      for (size_t i = 0; i < n; ++i) {
         for (size_t j = 0; j < i; ++j) {
            if (uniform(generator) < 0.2) {
               matrix.add(i, j, std::pow(10., exponent(generator)));
            }
         }
      }
      if (!matrix.values.empty()) {
         expect_properties(matrix);
      }
   }
}

TEST(MatchingScaling, FillGrowthOptions) {
   Options options;
   DefaultOptions::load(options);
   EXPECT_EQ(options.get_string("linear_system_scaling"), "ruiz");
   EXPECT_DOUBLE_EQ(options.get_double("linear_system_scaling_fill_growth_threshold"), 2.);
}
