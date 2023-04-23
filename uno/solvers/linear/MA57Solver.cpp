// Copyright (c) 2018-2023 Charlie Vanaret
// Licensed under the MIT license. See LICENSE file in the project directory for details.

#include <iostream>
#include <cassert>
#include "MA57Solver.hpp"
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

MA57Solver::MA57Solver(size_t max_dimension, size_t max_number_nonzeros) : SymmetricIndefiniteLinearSolver<double>(max_dimension),
   iwork(5 * max_dimension),
   lwork(static_cast<int>(1.2 * static_cast<double>(max_dimension))),
   work(static_cast<size_t>(this->lwork)), residuals(max_dimension) {
   this->row_indices.reserve(max_number_nonzeros);
   this->column_indices.reserve(max_number_nonzeros);
   // set the default values of the controlling parameters
   ma57id_(this->cntl.data(), this->icntl.data());
   // suppress warning messages
   this->icntl[4] = 0;
   // iterative refinement enabled
   this->icntl[8] = 1;
}

void MA57Solver::factorize(const SymmetricMatrix<double>& matrix) {
   // general factorization method: symbolic factorization and numerical factorization
   this->do_symbolic_factorization(matrix);
   this->do_numerical_factorization(matrix);
}

void MA57Solver::do_symbolic_factorization(const SymmetricMatrix<double>& matrix) {
   assert(matrix.dimension <= this->max_dimension && "MA57Solver: the dimension of the matrix is larger than the preallocated size");
   assert(matrix.number_nonzeros <= this->row_indices.capacity() &&
      "MA57Solver: the number of nonzeros of the matrix is larger than the preallocated size");

   // build the internal matrix representation
   this->save_matrix_to_local_format(matrix);

   const int n = static_cast<int>(matrix.dimension);
   const int nnz = static_cast<int>(matrix.number_nonzeros);

   // sparsity pattern
   const int lkeep = 5 * n + nnz + std::max(n, nnz) + 42;
   std::vector<int> keep(static_cast<size_t>(lkeep));

   // symbolic factorization
   ma57ad_(/* const */ &n,
         /* const */ &nnz,
         /* const */ this->row_indices.data(),
         /* const */ this->column_indices.data(),
         /* const */ &lkeep,
         /* const */ keep.data(),
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
   this->factorization = {n, nnz, std::move(fact), lfact, std::move(ifact), lifact, lkeep, std::move(keep)};
}

void MA57Solver::do_numerical_factorization(const SymmetricMatrix<double>& matrix) {
   assert(matrix.dimension <= this->max_dimension && "MA57Solver: the dimension of the matrix is larger than the preallocated size");
   assert(this->factorization.nnz == static_cast<int>(matrix.number_nonzeros) && "MA57Solver: the numbers of nonzeros do not match");

   const int n = static_cast<int>(matrix.dimension);
   // numerical factorization
   ma57bd_(&n,
         &this->factorization.nnz,
         /* const */ matrix.data_raw_pointer(),
         /* out */ this->factorization.fact.data(),
         /* const */ &this->factorization.lfact,
         /* out*/ this->factorization.ifact.data(),
         /* const */ &this->factorization.lifact,
         /* const */ &this->factorization.lkeep,
         /* const */ this->factorization.keep.data(), this->iwork.data(), this->icntl.data(), this->cntl.data(),
         /* out */ this->info.data(),
         /* out */ this->rinfo.data());
}

void MA57Solver::solve_indefinite_system(const SymmetricMatrix<double>& matrix, const std::vector<double>& rhs, std::vector<double>& result) {
   // solve
   const int n = static_cast<int>(matrix.dimension);
   const int lrhs = n; // integer, length of rhs

   // solve the linear system
   if (this->use_iterative_refinement) {
      ma57dd_(&this->job, &n, &this->factorization.nnz, matrix.data_raw_pointer(), this->row_indices.data(), this->column_indices.data(),
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

void MA57Solver::save_matrix_to_local_format(const SymmetricMatrix<double>& matrix) {
   // build the internal matrix representation
   this->row_indices.clear();
   this->column_indices.clear();
   matrix.for_each([&](size_t i, size_t j, double /*entry*/) {
      this->row_indices.push_back(static_cast<int>(i + this->fortran_shift));
      this->column_indices.push_back(static_cast<int>(j + this->fortran_shift));
   });
}