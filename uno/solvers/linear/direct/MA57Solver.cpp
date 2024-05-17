// Copyright (c) 2018-2024 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <cassert>
#include "MA57Solver.hpp"
#include "linear_algebra/SymmetricMatrix.hpp"
#include "linear_algebra/Vector.hpp"

extern "C" {
// MA57
// default values of controlling parameters
void ma57id_(double cntl[], int icntl[]);
// symbolic factorization
void ma57ad_(const int* n, const int* ne, const int irn[], const int jcn[], const int* lkeep, int keep[], int iwork[], int icntl[], int info[], double
rinfo[]);
// numerical factorization
void ma57bd_(const int* n, int* ne, const double a[], /* out */ double fact[], const int* lfact, /* out */ int ifact[], const int* lifact,
      const int* lkeep, const int keep[], int iwork[], int icntl[], double cntl[], /* out */ int info[], /* out */ double rinfo[]);
// linear system solve without iterative refinement
void ma57cd_(const int* job, const int* n, double fact[], int* lfact, int ifact[], int* lifact, const int* nrhs, double rhs[], const int* lrhs, double
work[], int* lwork, int iwork[], int icntl[], int info[]);
// linear system solve with iterative refinement
void ma57dd_(const int* job, const int* n, int* ne, const double a[], const int irn[], const int jcn[], double fact[], int* lfact, int ifact[], int*
lifact, const double rhs[], double x[], double resid[], double work[], int iwork[], int icntl[],
      double cntl[], int info[], double rinfo[]);
}

MA57Solver::MA57Solver(size_t dimension, size_t number_nonzeros) : DirectIndefiniteLinearSolver<size_t, double>(dimension),
      COO_matrix(dimension, number_nonzeros, false),
      hessian(this->COO_matrix.view(0, 0)),
      jacobian(this->COO_matrix.view(0, 0)),
      lkeep(static_cast<int>(5 * dimension + number_nonzeros + std::max(dimension, number_nonzeros) + 42)),
      keep(static_cast<size_t>(lkeep)),
      iwork(5 * dimension),
      lwork(static_cast<int>(1.2 * static_cast<double>(dimension))),
      work(static_cast<size_t>(this->lwork)), residuals(dimension) {
   // set the default values of the controlling parameters
   ma57id_(this->cntl.data(), this->icntl.data());
   // suppress warning messages
   this->icntl[4] = 0;
   // iterative refinement enabled
   this->icntl[8] = 1;
}

void MA57Solver::do_symbolic_factorization(const SymmetricMatrix<size_t, double>& matrix) {
   assert(matrix.dimension <= this->dimension && "MA57Solver: the dimension of the matrix is larger than the preallocated size");
   assert(matrix.number_nonzeros <= this->COO_matrix.number_nonzeros &&
      "MA57Solver: the number of nonzeros of the matrix is larger than the preallocated size");

   // build the internal matrix representation
   this->save_matrix_to_local_format(matrix);

   const int n = static_cast<int>(matrix.dimension);
   const int nnz = static_cast<int>(matrix.number_nonzeros);

   // symbolic factorization
   ma57ad_(/* const */ &n,
         /* const */ &nnz,
         /* const */ this->COO_matrix.row_indices_pointer(),
         /* const */ this->COO_matrix.column_indices_pointer(),
         /* const */ &this->lkeep,
         /* const */ this->keep.data(),
         /* out */ this->iwork.data(),
         /* const */ this->icntl.data(),
         /* out */ this->info.data(),
         /* out */ this->rinfo.data());

   assert(0 <= info[0] && "MA57: the symbolic factorization failed");
   if (0 < info[0]) {
      WARNING << "MA57 has issued a warning: info(1) = " << info[0] << '\n';
   }
   int lfact = 2 * this->info[8];
   std::vector<double> fact(static_cast<size_t>(lfact));
   int lifact = 2 * this->info[9];
   std::vector<int> ifact(static_cast<size_t>(lifact));

   // store the symbolic factorization
   this->factorization = {n, nnz, std::move(fact), lfact, std::move(ifact), lifact};
}

void MA57Solver::do_numerical_factorization(const SymmetricMatrix<size_t, double>& matrix) {
   assert(matrix.dimension <= this->dimension && "MA57Solver: the dimension of the matrix is larger than the preallocated size");
   assert(this->factorization.nnz == static_cast<int>(matrix.number_nonzeros) && "MA57Solver: the numbers of nonzeros do not match");

   const int n = static_cast<int>(matrix.dimension);
   int nnz = static_cast<int>(matrix.number_nonzeros);

   // numerical factorization
   ma57bd_(&n,
         &nnz,
         /* const */ this->COO_matrix.data_raw_pointer(),
         /* out */ this->factorization.fact.data(),
         /* const */ &this->factorization.lfact,
         /* out */ this->factorization.ifact.data(),
         /* const */ &this->factorization.lifact,
         /* const */ &this->lkeep,
         /* const */ this->keep.data(), this->iwork.data(), this->icntl.data(), this->cntl.data(),
         /* out */ this->info.data(),
         /* out */ this->rinfo.data());
}

void MA57Solver::solve_indefinite_system(const SymmetricMatrix<size_t, double>& matrix, const std::vector<double>& rhs, std::vector<double>& result,
      bool from_scratch) {
   if (from_scratch) {
      this->do_symbolic_factorization(matrix);
   }
   this->do_numerical_factorization(matrix);

   // solve
   const int n = static_cast<int>(matrix.dimension);
   int nnz = static_cast<int>(matrix.number_nonzeros);
   const int lrhs = n; // integer, length of rhs

   // solve the linear system
   if (this->use_iterative_refinement) {
      ma57dd_(&this->job, &n, &nnz, matrix.data_raw_pointer(), this->COO_matrix.row_indices_pointer(), this->COO_matrix.column_indices_pointer(),
            this->factorization.fact.data(), &this->factorization.lfact, this->factorization.ifact.data(), &this->factorization.lifact,
            rhs.data(), result.data(), this->residuals.data(), this->work.data(), this->iwork.data(), this->icntl.data(),
            this->cntl.data(), this->info.data(), this->rinfo.data());
   }
   else {
      // copy rhs into result (overwritten by MA57)
      copy_from(result, rhs);

      ma57cd_(&this->job, &n, this->factorization.fact.data(), &this->factorization.lfact, this->factorization.ifact.data(),
            &this->factorization.lifact, &this->nrhs, result.data(), &lrhs, this->work.data(), &this->lwork, this->iwork.data(),
            this->icntl.data(), this->info.data());
   }
}

std::tuple<size_t, size_t, size_t> MA57Solver::get_inertia() const {
   // rank = number_positive_eigenvalues + number_negative_eigenvalues
   // n = rank + number_zero_eigenvalues
   const size_t rank = this->rank();
   const size_t number_negative_eigenvalues = this->number_negative_eigenvalues();
   const size_t number_positive_eigenvalues = rank - number_negative_eigenvalues;
   const size_t number_zero_eigenvalues = static_cast<size_t>(this->factorization.n) - rank;
   return std::make_tuple(number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues);
}

size_t MA57Solver::number_negative_eigenvalues() const {
   return static_cast<size_t>(this->info[23]);
}

/*
bool MA57Solver::matrix_is_positive_definite() const {
   // positive definite = non-singular and no negative eigenvalues
   return not this->matrix_is_singular() && this->number_negative_eigenvalues() == 0;
}
*/

bool MA57Solver::matrix_is_singular() const {
   return (this->info[0] == 4);
}

size_t MA57Solver::rank() const {
   return static_cast<size_t>(this->info[24]);
}

void MA57Solver::save_matrix_to_local_format(const SymmetricMatrix<size_t, double>& matrix) {
   // build the internal matrix representation
   this->COO_matrix.reset();
   this->COO_matrix.dimension = matrix.dimension;
   matrix.for_each([&](size_t row_index, size_t column_index, double element) {
      this->COO_matrix.insert(element, static_cast<int>(row_index + this->fortran_shift), static_cast<int>(column_index + this->fortran_shift));
   });
}