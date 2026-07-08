// Copyright (c) 2026 Alexis Montoison and Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA86SOLVER_H
#define UNO_MA86SOLVER_H

#include <vector>
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "../COOLinearSystem.hpp"
#include "ingredients/subproblem_solvers/HSL/HSLLoader.hpp"
#include "linear_algebra/Indexing.hpp"

namespace uno {
   // Interface to HSL MA86, a supernodal sparse symmetric indefinite linear solver. Unlike
   // MA27/MA57, MA86 exposes a C interface built around opaque control/info structs and a
   // persistent factor handle (void* keep), and it expects the lower triangle of the matrix in
   // CSC format together with an externally-computed fill-reducing ordering (via MC68). Uno
   // assembles the matrix in COO format, so a COO -> CSC conversion (summing duplicate entries)
   // is performed once during the symbolic analysis and reused afterwards.
   class MA86Solver: public DirectSymmetricIndefiniteLinearSolver<double> {
   public:
      // solver_indexing selects the base of the COO matrix Uno assembles (C/0-based by default,
      // or Fortran/1-based); the CSC is built in that base and passed to MA86/MC68 via f_arrays
      explicit MA86Solver(int solver_indexing = Indexing::C_indexing);
      ~MA86Solver() override;

      void initialize_memory() override;

      void do_symbolic_analysis() override;
      void do_numerical_factorization(bool is_matrix_positive_definite) override;
      void solve_indefinite_system(double* solution) override;
      void solve_indefinite_system(const double* rhs, double* solution, size_t number_of_rhs) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

      [[nodiscard]] LinearSystem& get_linear_system() override;
      [[nodiscard]] COOLinearSystem& get_coo_linear_system();

   protected:
      // MA86 control/info and the opaque factor handle
      ma86_control control{};
      ma86_info info{};
      void* keep{nullptr};
      int n{0};

      // base of the COO matrix Uno assembles (Indexing::C_indexing or Indexing::Fortran_indexing)
      const int solver_indexing;
      // Uno assembles the matrix in COO format; build_csc_from_coo converts it to the 0-based CSC below
      COOLinearSystem linear_system;

      // lower triangle in CSC format (in the solver_indexing base) expected by MA86
      int number_nonzeros{0}; // number of unique entries after duplicate summation
      std::vector<int> column_pointers{}; // size n + 1
      std::vector<int> row_indices{};     // size number_nonzeros
      std::vector<double> values{};       // size number_nonzeros
      std::vector<int> order{};           // fill-reducing ordering from MC68 (size n)

      // maps each COO entry to its slot in the CSC arrays (size = COO number of nonzeros)
      std::vector<int> coo_to_csc{};

      bool analysis_performed{false};
      bool factorization_performed{false};

      void build_csc_from_coo();
   };
} // namespace

#endif // UNO_MA86SOLVER_H
