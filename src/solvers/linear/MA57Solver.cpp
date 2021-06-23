#include <iostream>
#include "MA57Solver.hpp"
#include "Vector.hpp"
#include <cassert>

extern "C" {
// MA57
// default values of controlling parameters
void ma57id_(double cntl[], int icntl[]);
// symbolic factorization
void
ma57ad_(int* n, int* ne, const int irn[], const int jcn[], int* lkeep, int keep[], int iwork[], int icntl[], int info[], double rinfo[]);
// numerical factorization
void ma57bd_(int* n, int* ne, const double a[], /* out */ double fact[], const int* lfact, /* out */ int ifact[], const int* lifact,
      const int* lkeep, const int keep[], int iwork[], int icntl[], double cntl[], /* out */ int info[], /* out */ double rinfo[]);
// linear system solve without iterative refinement
void ma57cd_(int* job, int* n, double fact[], int* lfact, int ifact[], int* lifact, int* nrhs, double rhs[], int* lrhs, double work[],
      int* lwork, int iwork[], int icntl[], int info[]);
// linear system solve with iterative refinement
void ma57dd_(int* job, int* n, int* ne, const double a[], const int irn[], const int jcn[], double fact[], int* lfact, int ifact[], int* lifact,
      double rhs[], double x[], double resid[], double work[], int iwork[], int icntl[],
      double cntl[], int info[], double rinfo[]);
}

MA57Solver::MA57Solver() : use_fortran(1), cntl_(5), icntl_(20), rinfo_(20) {
}

void MA57Solver::factorize(COOMatrix& matrix) {
   // general factorization method: symbolic factorization and numerical factorization
   // for more flexibility, call the two methods independently
   this->do_symbolic_factorization(matrix);
   this->do_numerical_factorization(matrix);
}

void MA57Solver::do_symbolic_factorization(const COOMatrix& matrix) {
   assert(matrix.fortran_indexing == this->use_fortran && "MA57Solver: please use the correct Fortran indexing");

   size_t n = matrix.dimension;
   size_t nnz = matrix.number_nonzeros;

   /* set the default values of the controlling parameters */
   ma57id_(this->cntl_.data(), this->icntl_.data());
   // suppress warning messages
   this->icntl_[4] = 0;
   // iterative refinement enabled
   this->icntl_[8] = 1;

   /* sparsity pattern */
   int lkeep = 5 * n + nnz + std::max(n, nnz) + 42;
   std::vector<int> keep(lkeep);
   std::vector<int> iwork(5 * n);

   /* symbolic factorization */
   std::vector<int> info(40);
   ma57ad_(/* const */ (int*) &n,
         /* const */ (int*) &nnz,
         /* const */ matrix.row_indices.data(),
         /* const */ matrix.column_indices.data(),
         /* const */ &lkeep,
         /* const */ keep.data(),
         /* out */ iwork.data(),
         /* const */ this->icntl_.data(),
         /* out */ info.data(),
         /* out */ this->rinfo_.data());

   // TODO check that info[0] == 0

   int lfact = 2 * info[8];
   std::vector<double> fact(lfact);
   int lifact = 2 * info[9];
   std::vector<int> ifact(lifact);
   // build the factorization object
   this->factorization_ = {n, nnz, fact, lfact, ifact, lifact, lkeep, keep, iwork, info};
}

void MA57Solver::do_numerical_factorization(const COOMatrix& matrix) {
   assert(matrix.fortran_indexing == this->use_fortran && "MA57Solver: please use the correct Fortran indexing");
   assert(this->factorization_.n == matrix.dimension && "MA57Solver: the dimensions do not match");
   assert(this->factorization_.nnz == matrix.number_nonzeros && "MA57Solver: the numbers of nonzeros do not match");

   /* numerical factorization */
   ma57bd_((int*) &this->factorization_.n,
         (int*) &this->factorization_.nnz,
         /* const */ matrix.matrix.data(),
         /* out */ this->factorization_.fact.data(),
         /* const */ &this->factorization_.lfact,
         /* out*/ this->factorization_.ifact.data(),
         /* const */ &this->factorization_.lifact,
         /* const */ &this->factorization_.lkeep,
         /* const */ this->factorization_.keep.data(), this->factorization_.iwork.data(), this->icntl_.data(), this->cntl_.data(),
         /* out */ this->factorization_.info.data(),
         /* out */ this->rinfo_.data());
}

std::vector<double> MA57Solver::solve(const COOMatrix& matrix, std::vector<double>& rhs) {
   /* solve */
   int n = this->factorization_.n;
   int nrhs = 1; // number of right hand side being solved
   int lrhs = n; // integer, length of rhs
   int lwork = 1.2 * n * nrhs; // length of w; lw >= n*nrhs
   std::vector<double> work(lwork);

   // solve the linear system
   if (this->use_iterative_refinement) {
      int job = 1;
      std::vector<double> x(n);
      std::vector<double> residuals(n);

      ma57dd_(&job, &n, (int*) &this->factorization_.nnz, matrix.matrix.data(), matrix.row_indices.data(), matrix.column_indices.data(),
         this->factorization_.fact.data(), &this->factorization_.lfact, this->factorization_.ifact.data(), &this->factorization_.lifact,
         rhs.data(), x.data(), residuals.data(), work.data(), this->factorization_.iwork.data(), this->icntl_.data(),
         this->cntl_.data(), this->factorization_.info.data(), this->rinfo_.data());
      return x;
   }
   else {
      int job = 1;
      ma57cd_(&job, &n, this->factorization_.fact.data(), &this->factorization_.lfact, this->factorization_.ifact.data(),
            &this->factorization_.lifact, &nrhs, rhs.data(), &lrhs, work.data(), &lwork, this->factorization_.iwork.data(),
            this->icntl_.data(), this->factorization_.info.data());
      // the solution is copied in rhs
      return rhs;
   }
}

std::tuple<int, int, int> MA57Solver::get_inertia() const {
   int rank = this->rank();
   int number_negative_eigenvalues = this->number_negative_eigenvalues();
   int number_positive_eigenvalues = rank - number_negative_eigenvalues;
   int number_zero_eigenvalues = this->factorization_.n - rank;
   return std::make_tuple(number_positive_eigenvalues, number_negative_eigenvalues, number_zero_eigenvalues);
}

size_t MA57Solver::number_negative_eigenvalues() const {
   return this->factorization_.info[23];
}

bool MA57Solver::matrix_is_singular() const {
   return (this->factorization_.info[0] == 4);
}

int MA57Solver::rank() const {
   return this->factorization_.info[24];
}

/*
const int rank = info[24];
const int n_negative_eigenvalues = info[23];
const int n_positive_eigenvalues = rank - n_negative_eigenvalues;
assert(0 <= n_positive_eigenvalues);

return {n_positive_eigenvalues, n_negative_eigenvalues, n - rank};
*/
