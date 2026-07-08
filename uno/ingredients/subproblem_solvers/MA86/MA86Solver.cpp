// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <algorithm>
#include <cassert>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "MA86Solver.hpp"
#include "linear_algebra/VectorView.hpp"

#ifdef HSL_RUNTIME_LOADING
// route the calls through the runtime-resolved function pointers
#define MA86_default_control uno::hsl_ma86_default_control
#define MA86_analyse uno::hsl_ma86_analyse
#define MA86_factor uno::hsl_ma86_factor
#define MA86_solve uno::hsl_ma86_solve
#define MA86_finalise uno::hsl_ma86_finalise
#define MC68_default_control uno::hsl_mc68_default_control
#define MC68_order uno::hsl_mc68_order
#else
#define MA86_default_control ma86_default_control_d
#define MA86_analyse ma86_analyse_d
#define MA86_factor ma86_factor_d
#define MA86_solve ma86_solve_d
#define MA86_finalise ma86_finalise_d
#define MC68_default_control mc68_default_control_i
#define MC68_order mc68_order_i
#endif

namespace uno {
#ifndef HSL_RUNTIME_LOADING
   extern "C" {
      void MA86_default_control(ma86_control* control);
      void MA86_analyse(const int n, const int ptr[], const int row[], int order[], void** keep,
         const ma86_control* control, ma86_info* info);
      void MA86_factor(const int n, const int ptr[], const int row[], const double val[], const int order[],
         void** keep, const ma86_control* control, ma86_info* info, const double scale[]);
      void MA86_solve(const int job, const int nrhs, const int ldx, double* x, const int order[], void** keep,
         const ma86_control* control, ma86_info* info, const double scale[]);
      void MA86_finalise(void** keep, const ma86_control* control);
      void MC68_default_control(mc68_control* control);
      void MC68_order(const int ord, const int n, const int ptr[], const int row[], int perm[],
         const mc68_control* control, mc68_info* info);
   }
#endif

   // MC68 ordering algorithms (first argument of mc68_order)
   constexpr int MC68_AMD = 1;

   // ma86_solve job: 0 solves the full system A x = b
   constexpr int MA86_SOLVE_FULL_SYSTEM = 0;

   MA86Solver::MA86Solver(int solver_indexing):
         DirectSymmetricIndefiniteLinearSolver<double>(),
         solver_indexing(solver_indexing),
         linear_system(solver_indexing) {
#ifdef HSL_RUNTIME_LOADING
      if (!ma86_symbols_available()) {
         throw std::runtime_error("Uno: the MA86 solver was requested but the HSL library could not be loaded at "
            "runtime, or it does not export the MA86/MC68 C interface (set the UNO_HSL_LIBRARY environment variable "
            "to point at a libhsl providing ma86_*_d and mc68_*_i)");
      }
#endif
      // set the default values of the controlling parameters
      MA86_default_control(&this->control);
      // build_csc_from_coo emits the CSC in the COO base (solver_indexing); MA86's C interface
      // interprets ptr/row/order as 1-based when f_arrays != 0
      this->control.f_arrays = this->solver_indexing;
      this->control.diagnostics_level = -1; // no printing
      // Uno performs its own inertia correction, so MA86 must keep going even if the matrix is
      // singular (otherwise the inertia cannot be reported) and must not scale the system.
      this->control.action = 1;
      this->control.scaling = 0;
   }

   MA86Solver::~MA86Solver() {
      if (this->keep != nullptr) {
         MA86_finalise(&this->keep, &this->control);
      }
   }

   void MA86Solver::initialize_memory() {
      this->n = static_cast<int>(this->linear_system.dimension);
   }

   // convert the lower-triangular COO matrix assembled by Uno into the CSC format expected by MA86,
   // merging duplicate entries and inserting a (possibly zero) explicit diagonal in every column.
   // The CSC pointer/row arrays are emitted in the COO base (solver_indexing), so MA86/MC68 must be
   // told about it via f_arrays / f_array_in / f_array_out.
   void MA86Solver::build_csc_from_coo() {
      const size_t coo_number_nonzeros = this->linear_system.number_nonzeros;
      const auto& coo_rows = this->linear_system.matrix_row_indices;
      const auto& coo_columns = this->linear_system.matrix_column_indices;
      const uno_int base = static_cast<uno_int>(this->solver_indexing);

      // sort the COO entries by (column, row): this is the order of a lower-triangular CSC matrix
      std::vector<size_t> sorted(coo_number_nonzeros);
      std::iota(sorted.begin(), sorted.end(), 0);
      std::sort(sorted.begin(), sorted.end(), [&](size_t i, size_t j) {
         if (coo_columns[i] != coo_columns[j]) {
            return coo_columns[i] < coo_columns[j];
         }
         return coo_rows[i] < coo_rows[j];
      });

      // build the unique CSC sparsity pattern and the COO -> CSC mapping. Duplicate (row, column)
      // entries are merged (their values are summed during the numerical factorization), and each
      // column is given an explicit diagonal entry (with the smallest row index in a lower-triangular
      // column, hence inserted first), even if zero. Row/column indices stay in the COO base.
      this->row_indices.clear();
      this->row_indices.reserve(coo_number_nonzeros + static_cast<size_t>(this->n));
      this->coo_to_csc.assign(coo_number_nonzeros, 0);
      this->column_pointers.assign(static_cast<size_t>(this->n) + 1, 0);
      int current_slot = -1; // 0-based slot in the row_indices/values arrays (base-independent)
      size_t sorted_index = 0;
      for (uno_int column = base; column < static_cast<uno_int>(this->n) + base; ++column) {
         const size_t column_slot = static_cast<size_t>(column - base); // 0-based column position
         uno_int previous_row = base - 1; // sentinel: below any valid row in this base
         bool diagonal_present = false;
         while (sorted_index < coo_number_nonzeros && coo_columns[sorted[sorted_index]] == column) {
            const size_t entry = sorted[sorted_index];
            const uno_int row = coo_rows[entry];
            // insert the (column, column) diagonal before the first strictly sub-diagonal entry
            if (!diagonal_present && row > column) {
               ++current_slot;
               this->row_indices.push_back(static_cast<int>(column));
               ++this->column_pointers[column_slot + 1];
               diagonal_present = true;
               previous_row = column;
            }
            if (row != previous_row) {
               // new unique entry
               ++current_slot;
               this->row_indices.push_back(static_cast<int>(row));
               ++this->column_pointers[column_slot + 1];
               previous_row = row;
               if (row == column) {
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
            this->row_indices.push_back(static_cast<int>(column));
            ++this->column_pointers[column_slot + 1];
         }
      }
      this->number_nonzeros = current_slot + 1;
      this->values.resize(static_cast<size_t>(this->number_nonzeros));
      // prefix sum to turn per-column counts into column starts (0-based offsets)
      for (size_t column = 0; column < static_cast<size_t>(this->n); ++column) {
         this->column_pointers[column + 1] += this->column_pointers[column];
      }
      assert(this->column_pointers[static_cast<size_t>(this->n)] == this->number_nonzeros);
      // shift the column pointers into the COO base to match the row indices (no-op when base == 0)
      for (size_t column = 0; column <= static_cast<size_t>(this->n); ++column) {
         this->column_pointers[column] += base;
      }
   }

   void MA86Solver::do_symbolic_analysis() {
      assert(!this->analysis_performed);

      this->build_csc_from_coo();

      // compute a fill-reducing ordering with MC68 (AMD); MA86 requires this ordering as input
      this->order.resize(static_cast<size_t>(this->n));
      mc68_control control68{};
      mc68_info info68{};
      MC68_default_control(&control68);
      control68.f_array_in = this->solver_indexing;  // the CSC arrays are in the COO base
      control68.f_array_out = this->solver_indexing; // emit the ordering in the same base
      MC68_order(MC68_AMD, this->n, this->column_pointers.data(), this->row_indices.data(), this->order.data(),
         &control68, &info68);
      if (info68.flag < 0) {
         throw std::runtime_error("MC68 could not compute the ordering for MA86 (flag = " + std::to_string(info68.flag) + ")");
      }

      MA86_analyse(this->n, this->column_pointers.data(), this->row_indices.data(), this->order.data(), &this->keep,
         &this->control, &this->info);
      if (this->info.flag < 0) {
         throw std::runtime_error("MA86 could not compute the symbolic analysis (flag = " + std::to_string(this->info.flag) + ")");
      }
      this->analysis_performed = true;
   }

   void MA86Solver::do_numerical_factorization(bool /*is_matrix_positive_definite*/) {
      assert(this->analysis_performed);

      // scatter the COO values into the CSC values, summing duplicate entries
      std::fill(this->values.begin(), this->values.end(), 0.);
      for (size_t k = 0; k < this->linear_system.number_nonzeros; ++k) {
         this->values[static_cast<size_t>(this->coo_to_csc[k])] += this->linear_system.matrix_values[k];
      }

      MA86_factor(this->n, this->column_pointers.data(), this->row_indices.data(), this->values.data(),
         this->order.data(), &this->keep, &this->control, &this->info, nullptr);
      if (this->info.flag < 0) {
         throw std::runtime_error("MA86 could not compute the numerical factorization (flag = " + std::to_string(this->info.flag) + ")");
      }
      this->factorization_performed = true;
   }

   void MA86Solver::solve_indefinite_system(double* solution) {
      return this->solve_indefinite_system(this->linear_system.rhs.data(), solution, 1);
   }

   void MA86Solver::solve_indefinite_system(const double* rhs, double* solution, size_t number_of_rhs) {
      assert(this->factorization_performed);

      // copy the rhs block into the solution block (overwritten by MA86); both are column-major with
      // leading dimension n
      const size_t total_size = static_cast<size_t>(this->n) * number_of_rhs;
      view(solution, total_size) = view(rhs, total_size);

      MA86_solve(MA86_SOLVE_FULL_SYSTEM, static_cast<int>(number_of_rhs), this->n, solution, this->order.data(),
         &this->keep, &this->control, &this->info, nullptr);
      if (this->info.flag < 0) {
         throw std::runtime_error("MA86 could not solve the linear system (flag = " + std::to_string(this->info.flag) + ")");
      }
   }

   Inertia MA86Solver::get_inertia() const {
      // rank = number_positive_eigenvalues + number_negative_eigenvalues
      // n = rank + number_zero_eigenvalues
      const size_t rank = this->rank();
      const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
      const size_t number_positive_eigenvalues = rank - number_negative_eigenvalues;
      const size_t number_zero_eigenvalues = static_cast<size_t>(this->n) - rank;
      return {number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues};
   }

   size_t MA86Solver::number_negative_eigenvalues() const {
      return static_cast<size_t>(this->info.num_neg);
   }

   bool MA86Solver::matrix_is_singular() const {
      return (this->info.matrix_rank < this->n);
   }

   size_t MA86Solver::rank() const {
      return static_cast<size_t>(this->info.matrix_rank);
   }

   LinearSystem& MA86Solver::get_linear_system() {
      return this->linear_system;
   }

   COOLinearSystem& MA86Solver::get_coo_linear_system() {
      return this->linear_system;
   }
} // namespace