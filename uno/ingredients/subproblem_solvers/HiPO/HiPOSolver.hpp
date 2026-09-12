// Copyright (c) 2026 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_HIPOSOLVER_H
#define UNO_HIPOSOLVER_H

#include <vector>
#include "FactorHighs_c_api.h" // vendored copy of HiGHS's HiPO C API header (HiGHS does not install it)
#include "../DirectSymmetricIndefiniteLinearSolver.hpp"
#include "../COOLinearSystem.hpp"
#include "linear_algebra/Indexing.hpp"

namespace uno {
   // Interface to HiPO, the sparse symmetric indefinite linear solver shipped with HiGHS
   // (highs/ipm/hipo/factorhighs). HiPO expects the lower triangle of the matrix in CSC format,
   // whereas Uno assembles the matrix in COO format: a COO -> CSC conversion (with summation of
   // duplicate entries) is performed once during the symbolic analysis and reused afterwards.
   class HiPOSolver: public DirectSymmetricIndefiniteLinearSolver<double> {
   public:
      HiPOSolver();
      ~HiPOSolver() override;

      void initialize_memory() override;

      void set_expected_inertia(const Inertia& expected_inertia) override;
      void do_symbolic_analysis() override;
      void do_numerical_factorization(bool is_matrix_positive_definite) override;
      void solve_indefinite_system(double* result) override;
      // native multiple right-hand side solve (HiPO solves all columns in one call)
      void solve_indefinite_system(const double* rhs, double* solution, size_t number_of_rhs) override;

      [[nodiscard]] Inertia get_inertia() const override;
      [[nodiscard]] size_t number_negative_eigenvalues() const override;
      [[nodiscard]] bool matrix_is_singular() const override;
      [[nodiscard]] size_t rank() const override;

      [[nodiscard]] LinearSystem& get_linear_system() override;
      [[nodiscard]] COOLinearSystem& get_coo_linear_system();

   protected:
      // HiPO opaque objects: solver and symbolic factorization
      void* solver{nullptr};
      void* symbolic{nullptr};

      // Uno assembles the matrix in COO format (with C indexing)
      COOLinearSystem linear_system{Indexing::C_indexing};

      // lower triangle in CSC format (0-based) expected by HiPO
      HighsInt dimension{0};
      HighsInt number_nonzeros{0}; // number of unique entries after duplicate summation
      std::vector<HighsInt> column_pointers{};  // size dimension + 1
      std::vector<HighsInt> row_indices{};      // size number_nonzeros
      std::vector<double> values{};             // size number_nonzeros
      std::vector<HighsInt> pivot_signs{};      // expected sign of each pivot (size dimension)
      std::vector<HighsInt> permutation{};      // fill-reducing ordering (size dimension)

      // expected inertia (set by set_expected_inertia): the leading number_positive_pivots pivots are
      // expected positive and the next number_negative_pivots negative; used to fill pivot_signs
      size_t number_positive_pivots{0};
      size_t number_negative_pivots{0};

      // maps each COO entry to its slot in the CSC arrays (size = COO number of nonzeros)
      std::vector<HighsInt> coo_to_csc{};

      bool analysis_performed{false};
      bool factorization_performed{false};

      void compute_inertia(HighsInt& positive, HighsInt& negative, HighsInt& zero) const;
   };
} // namespace

#endif // UNO_HIPOSOLVER_H
