// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#ifndef UNO_MA57SOLVER_H
#define UNO_MA57SOLVER_H

#include <array>
#include <vector>
#include "DirectIndefiniteLinearSolver.hpp"
#include "linear_algebra/COOSymmetricMatrix.hpp"

// forward declaration
template <typename ElementType>
class Vector;

struct MA57Factorization {
   int n{};
   int nnz{};
   int lfact{};
   int lifact{};

   MA57Factorization() = default;
};

/*! \class MA57Solver
 * \brief Interface for MA57
 * see https://github.com/YimingYAN/linSolve
 *
 *  Interface to the symmetric indefinite linear solver MA57
 */
class MA57Solver : public DirectIndefiniteLinearSolver<size_t, double> {
public:
   MA57Solver(size_t dimension, size_t number_nonzeros);
   ~MA57Solver() override = default;

   void do_symbolic_factorization(const SymmetricMatrix<size_t, double>& matrix) override;
   void do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) override;
   void solve_indefinite_system(const SymmetricMatrix<size_t, double>& matrix, const Vector<double>& rhs, Vector<double>& result, bool from_scratch) override;

   [[nodiscard]] std::tuple<size_t, size_t, size_t> get_inertia() const override;
   [[nodiscard]] size_t number_negative_eigenvalues() const override;
   // [[nodiscard]] bool matrix_is_positive_definite() const override;
   [[nodiscard]] bool matrix_is_singular() const override;
   [[nodiscard]] size_t rank() const override;

private:
   // internal matrix representation: COO matrix for the (augmented) system + views for individual blocks
   COOSymmetricMatrix<int, double> COO_matrix;
   COOSymmetricMatrix<int, double>::View hessian;
   COOSymmetricMatrix<int, double>::View jacobian;

   // factorization
   MA57Factorization factorization{};
   std::vector<double> fact{0}; // do not initialize, resize at every iteration
   std::vector<int> ifact{0}; // do not initialize, resize at every iteration
   const int lkeep;
   std::vector<int> keep{};
   std::vector<int> iwork{};
   int lwork;
   std::vector<double> work{};

   // for ma57id_ (default values of controlling parameters)
   std::array<double, 5> cntl{};
   std::array<int, 20> icntl{};
   std::array<double, 20> rinfo{};
   std::array<int, 40> info{};

   const int nrhs{1}; // number of right hand side being solved
   const int job{1};
   std::vector<double> residuals;
   const size_t fortran_shift{1};

   bool use_iterative_refinement{false};
   void save_matrix_to_local_format(const SymmetricMatrix<size_t, double>& row_index);
};

#endif // UNO_MA57SOLVER_H
