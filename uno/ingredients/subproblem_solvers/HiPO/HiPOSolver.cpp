// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "HiPOSolver.hpp"
#include "linear_algebra/Vector.hpp"
#include "symbolic/Range.hpp"

namespace uno {
   HiPOSolver::HiPOSolver(): DirectSymmetricIndefiniteLinearSolver<double>() {
      // initialize the global HiPO scheduler once for the whole process
      static const HighsInt initialize_status = FactorHighs_initialise(0); // 0: default number of threads
      if (initialize_status != 0) {
         throw std::runtime_error("HiPO is not available in this HiGHS build");
      }
      this->symbolic = FactorHighs_symbolic_create();
      this->solver = FactorHighs_create();
      // Uno performs its own inertia correction, so HiPO must factorize the matrix as is:
      // no internal regularization, and report the true inertia
      FactorHighs_setLogging(this->solver, 0);              // no printing
      FactorHighs_setRegularisation(this->solver, 0.0, 0.0); // no static regularization
      FactorHighs_setPivoting(this->solver, 1);             // Bunch-Kaufman pivoting for stability
      FactorHighs_setOneIndexing(this->solver, 0);          // 0-based CSC
   }

   HiPOSolver::~HiPOSolver() {
      if (this->symbolic != nullptr) {
         FactorHighs_symbolic_destroy(this->symbolic);
      }
      if (this->solver != nullptr) {
         FactorHighs_destroy(this->solver);
      }
      // note: the global scheduler (FactorHighs_terminate) is shared between instances and is left
      // to be cleaned up at process exit
   }

   void HiPOSolver::initialize_memory() {
      this->dimension = static_cast<HighsInt>(this->linear_system.dimension);
   }

   void HiPOSolver::set_expected_inertia(const Inertia& expected_inertia) {
      // HiPO needs the expected sign of each pivot: the leading positive pivots, then the negative ones
      this->number_positive_pivots = expected_inertia.positive;
      this->number_negative_pivots = expected_inertia.negative;
   }

   void HiPOSolver::do_symbolic_analysis() {
      assert(!this->analysis_performed);

      const HighsInt n = this->dimension;
      const size_t coo_number_nonzeros = this->linear_system.number_nonzeros;
      const std::vector<uno_int>& coo_rows = this->linear_system.matrix_row_indices;
      const std::vector<uno_int>& coo_columns = this->linear_system.matrix_column_indices;

      // sort the COO entries by (column, row): this is the order of a lower-triangular CSC matrix
      std::vector<size_t> order(coo_number_nonzeros);
      std::iota(order.begin(), order.end(), 0);
      std::sort(order.begin(), order.end(), [&](size_t i, size_t j) {
         if (coo_columns[i] != coo_columns[j]) {
            return coo_columns[i] < coo_columns[j];
         }
         return coo_rows[i] < coo_rows[j];
      });

      // build the unique CSC sparsity pattern and the COO -> CSC mapping. Two HiPO requirements are
      // enforced here:
      //  - duplicate (row, column) entries are merged: HiPO does not sum them (their values are summed
      //    during the numerical factorization);
      //  - every column must contain an explicit diagonal entry, even if zero: HiPO's analyse requires
      //    it. Uno's augmented system has no diagonal on the dual (2, 2) block unless dual
      //    regularization is applied, so any missing diagonal is inserted here with a zero value.
      // In a lower-triangular CSC column the diagonal has the smallest row index, so it comes first.
      this->row_indices.clear();
      this->row_indices.reserve(coo_number_nonzeros + static_cast<size_t>(n));
      this->coo_to_csc.assign(coo_number_nonzeros, 0);
      this->column_pointers.assign(static_cast<size_t>(n) + 1, 0);
      HighsInt current_slot = -1;
      size_t sorted_index = 0;
      for (HighsInt column = 0; column < n; ++column) {
         uno_int previous_row = -1;
         bool diagonal_present = false;
         while (sorted_index < coo_number_nonzeros &&
               coo_columns[order[sorted_index]] == static_cast<uno_int>(column)) {
            const size_t entry = order[sorted_index];
            const uno_int row = coo_rows[entry];
            // insert the (column, column) diagonal before the first strictly sub-diagonal entry
            if (!diagonal_present && static_cast<HighsInt>(row) > column) {
               ++current_slot;
               this->row_indices.push_back(column);
               ++this->column_pointers[static_cast<size_t>(column) + 1];
               diagonal_present = true;
               previous_row = static_cast<uno_int>(column);
            }
            if (row != previous_row) {
               // new unique entry
               ++current_slot;
               this->row_indices.push_back(static_cast<HighsInt>(row));
               ++this->column_pointers[static_cast<size_t>(column) + 1];
               previous_row = row;
               if (static_cast<HighsInt>(row) == column) {
                  diagonal_present = true;
               }
            }
            // duplicate entries reuse current_slot
            this->coo_to_csc[entry] = current_slot;
            ++sorted_index;
         }
         // empty column, or column without an explicit diagonal: add a zero diagonal entry
         if (!diagonal_present) {
            ++current_slot;
            this->row_indices.push_back(column);
            ++this->column_pointers[static_cast<size_t>(column) + 1];
         }
      }
      this->number_nonzeros = current_slot + 1;
      this->values.resize(static_cast<size_t>(this->number_nonzeros));
      // prefix sum to turn per-column counts into column starts
      for (size_t column = 0; column < static_cast<size_t>(n); ++column) {
         this->column_pointers[column + 1] += this->column_pointers[column];
      }
      assert(this->column_pointers[static_cast<size_t>(n)] == this->number_nonzeros);

      // expected sign of each pivot, derived from the expected inertia. HiPO is designed for symmetric
      // quasi-definite systems: the leading positive pivots get +1 (primal (1, 1) block), the next
      // negative pivots get -1 (dual (2, 2) block), and any remaining pivot gets 0 (unknown sign).
      this->pivot_signs.assign(static_cast<size_t>(n), 0);
      size_t index = 0;
      for (size_t count = 0; count < this->number_positive_pivots && index < static_cast<size_t>(n); ++count, ++index) {
         this->pivot_signs[index] = 1;
      }
      for (size_t count = 0; count < this->number_negative_pivots && index < static_cast<size_t>(n); ++count, ++index) {
         this->pivot_signs[index] = -1;
      }

      // compute a fill-reducing ordering with METIS, falling back to the identity ordering
      this->permutation.resize(static_cast<size_t>(n));
      const HighsInt reorder_status = FactorHighs_reorderMetis(this->solver, n, this->number_nonzeros,
         this->row_indices.data(), this->column_pointers.data(), this->permutation.data());
      if (reorder_status != 0) {
         std::iota(this->permutation.begin(), this->permutation.end(), 0);
      }

      const HighsInt analyse_status = FactorHighs_analyse(this->solver, this->symbolic, n, this->number_nonzeros,
         this->row_indices.data(), this->column_pointers.data(), this->pivot_signs.data(), this->permutation.data());
      if (analyse_status != 0) {
         throw std::runtime_error("HiPO could not compute the symbolic analysis");
      }
      this->analysis_performed = true;
   }

   void HiPOSolver::do_numerical_factorization(bool /*is_matrix_positive_definite*/) {
      assert(this->analysis_performed);

      // scatter the COO values into the CSC values, summing duplicate entries
      std::fill(this->values.begin(), this->values.end(), 0.);
      for (size_t k = 0; k < this->linear_system.number_nonzeros; ++k) {
         this->values[static_cast<size_t>(this->coo_to_csc[k])] += this->linear_system.matrix_values[k];
      }

      const HighsInt factorise_status = FactorHighs_factorise(this->solver, this->symbolic, this->dimension,
         this->number_nonzeros, this->row_indices.data(), this->column_pointers.data(), this->values.data());
      if (factorise_status != 0) {
         throw std::runtime_error("HiPO could not compute the numerical factorization");
      }
      this->factorization_performed = true;
   }

   void HiPOSolver::solve_indefinite_system(double* result) {
      assert(this->factorization_performed);

      // copy rhs into result (overwritten by HiPO with the solution)
      for (size_t index: Range(static_cast<size_t>(this->dimension))) {
         result[index] = this->linear_system.rhs[index];
      }
      const HighsInt solve_status = FactorHighs_solve(this->solver, result, 1);
      if (solve_status != 0) {
         throw std::runtime_error("HiPO could not solve the linear system");
      }
   }

   void HiPOSolver::solve_indefinite_system(const double* rhs, double* solution, size_t number_of_rhs) {
      assert(this->factorization_performed);

      // HiPO solves in place, with the right-hand sides stored column-major (each column of length
      // dimension); copy the rhs block into the solution block and solve all columns at once.
      const size_t total_size = static_cast<size_t>(this->dimension) * number_of_rhs;
      for (size_t index = 0; index < total_size; ++index) {
         solution[index] = rhs[index];
      }
      const HighsInt solve_status = FactorHighs_solve(this->solver, solution, static_cast<HighsInt>(number_of_rhs));
      if (solve_status != 0) {
         throw std::runtime_error("HiPO could not solve the linear system");
      }
   }

   void HiPOSolver::compute_inertia(HighsInt& positive, HighsInt& negative, HighsInt& zero) const {
      constexpr double zero_tolerance = 1e-16;
      FactorHighs_inertia(this->solver, &positive, &negative, &zero, zero_tolerance);
   }

   Inertia HiPOSolver::get_inertia() const {
      HighsInt positive = 0, negative = 0, zero = 0;
      this->compute_inertia(positive, negative, zero);
      return {static_cast<size_t>(positive), static_cast<size_t>(negative), static_cast<size_t>(zero)};
   }

   size_t HiPOSolver::number_negative_eigenvalues() const {
      HighsInt positive = 0, negative = 0, zero = 0;
      this->compute_inertia(positive, negative, zero);
      return static_cast<size_t>(negative);
   }

   bool HiPOSolver::matrix_is_singular() const {
      HighsInt positive = 0, negative = 0, zero = 0;
      this->compute_inertia(positive, negative, zero);
      return (zero > 0);
   }

   size_t HiPOSolver::rank() const {
      HighsInt positive = 0, negative = 0, zero = 0;
      this->compute_inertia(positive, negative, zero);
      return static_cast<size_t>(this->dimension - zero);
   }

   LinearSystem& HiPOSolver::get_linear_system() {
      return this->linear_system;
   }

   COOLinearSystem& HiPOSolver::get_coo_linear_system() {
      return this->linear_system;
   }
} // namespace
